"""Times the neo4j cypher query held in ``cypher.map``.

The query expands from an ``F2`` node matching a *smiles* along ``FRAG``
relationships (up to *num_hops* levels) and then returns the ``Mol`` nodes
reached by a ``FRAG`` edge carrying the given *synthon*.

``README.md`` lists 5 SMILES/synthon pairs that respond quickly and 5 that do
not; both sets are reproduced here as :data:`FAST_PAIRS` and :data:`SLOW_PAIRS`.

Expose the graph first with::

    kubectl port-forward service/graph 7687 -n graph-x

then run, for example::

    uv run cypher-query-tester --set fast
"""

from __future__ import annotations

import argparse
import os
import statistics
import sys
import time
from dataclasses import dataclass, field

from neo4j import Driver, GraphDatabase
from neo4j.exceptions import Neo4jError

# The query exactly as held in cypher.map. Kept as a %-format template (rather
# than converted to cypher parameters) so that what we time is what the
# application actually sends to the server.
CYPHER_TEMPLATE = """
MATCH (a:F2 {smiles: '%(smiles)s'})
CALL apoc.path.expandConfig(a, {
  relationshipFilter: '<FRAG',
  minLevel: 0,
  maxLevel: %(num_hops)d,
  uniqueness: 'RELATIONSHIP_GLOBAL'
}) YIELD path
WITH last(nodes(path)) AS b
MATCH (b)<-[e:FRAG]-(c:Mol)
WHERE e.prop_synthon = '%(synthon)s'
RETURN c.smiles AS smi, c.cmpd_ids AS ids
LIMIT %(limit)d
"""

DEFAULT_URI = "bolt://localhost:7687"
DEFAULT_NUM_HOPS = 2
DEFAULT_LIMIT = 5
DEFAULT_REPEATS = 3


@dataclass(frozen=True)
class Pair:
    """A SMILES/synthon pair, as listed in README.md."""

    smiles: str
    synthon: str


# The 5 pairs README.md says generate fast results.
FAST_PAIRS: tuple[Pair, ...] = (
    Pair("CCNCC", "[Xe]c1cscn1"),
    Pair("CC(C)(CO)NC=O", "[Xe]n1cnc2ccccc21"),
    Pair("CC(C)(CO)NC=O", "[Xe]c1ncc1"),
    Pair("CCC=O", "[Xe]N1CCCC1"),
    Pair("C1CCCCC1", "[Xe]N1CCNCC1"),
)

# The 5 pairs README.md flags as problematic - these take too long to respond,
# so always run them with a --timeout.
SLOW_PAIRS: tuple[Pair, ...] = (
    Pair("c1cscn1", "CN(C)C(=O)CO[Xe]"),
    Pair("c1ccc2ncccc2c1", "[Xe]N1CCCc2ccccc21"),
    Pair("c1ccc2cnc2c1", "CN(C)C(=O)CO[Xe]"),
    Pair("CCC=O", "O=c1cccnn1[Xe]"),
    Pair("c1cscn1", "O=C1CC([Xe])=NN1"),
)

PAIR_SETS: dict[str, tuple[Pair, ...]] = {
    "fast": FAST_PAIRS,
    "slow": SLOW_PAIRS,
    "all": FAST_PAIRS + SLOW_PAIRS,
}


def build_query(
    smiles: str,
    synthon: str,
    num_hops: int = DEFAULT_NUM_HOPS,
    limit: int = DEFAULT_LIMIT,
) -> str:
    """Return the cypher.map query with its 4 variables substituted."""
    return CYPHER_TEMPLATE % {
        "smiles": smiles,
        "synthon": synthon,
        "num_hops": num_hops,
        "limit": limit,
    }


# --- Cost prediction -------------------------------------------------------
#
# The query is lazy: LIMIT 5 lets it abandon the expansion as soon as 5 edges
# carrying the synthon are found. So its cost is
#
#     min(edges scanned before 5 matches, edges in the whole 2-hop neighbourhood)
#
# A pair is therefore only slow when BOTH factors are bad - a synthon too rare
# to trigger the early exit AND a neighbourhood too large to exhaust cheaply.
# Both are measurable with bounded probes that never run the real expansion:
# degrees come from the degree store (O(1)) and the synthon hit rate from a
# LIMITed sample of the edges the query would filter.

# A neighbourhood smaller than this is cheap to exhaust even with no early exit.
# Observed: ~300k edges exhausts in 0.7s, while ~23M does not finish in 300s.
CHEAP_TO_EXHAUST_EDGES = 5_000_000

# Matches in the sample needed before we trust the early exit to fire.
CONFIDENT_HITS = 5

PROBE_NEIGHBOURS = 2_000
PROBE_EDGES = 50_000

DEGREE_PROBE = """
MATCH (a:F2 {smiles: $smiles})
WITH a, apoc.node.degree(a, '<FRAG') AS level1
MATCH (a)<-[:FRAG]-(b)
WITH level1, b LIMIT $neighbours
RETURN level1,
       count(b)                          AS sampled,
       avg(apoc.node.degree(b, '<FRAG')) AS avg_degree
"""

SYNTHON_PROBE = """
MATCH (a:F2 {smiles: $smiles})<-[:FRAG]-(b)
WITH b LIMIT $neighbours
MATCH (b)<-[e:FRAG]-(c:Mol)
WITH e LIMIT $edges
RETURN count(e) AS sampled,
       sum(CASE WHEN e.prop_synthon = $synthon THEN 1 ELSE 0 END) AS hits
"""


@dataclass
class Prediction:
    """A cheap estimate of whether one pair's query will be fast or slow."""

    pair: Pair
    exists: bool
    level1_degree: int = 0
    avg_level2_degree: float = 0.0
    sampled_edges: int = 0
    synthon_hits: int = 0
    probe_seconds: float = 0.0

    @property
    def estimated_edges(self) -> float:
        """Roughly how many FRAG edges a full 2-hop expansion would filter."""
        level2_nodes = self.level1_degree * self.avg_level2_degree
        return (self.level1_degree + level2_nodes) * self.avg_level2_degree

    @property
    def hit_rate(self) -> float:
        return self.synthon_hits / self.sampled_edges if self.sampled_edges else 0.0

    @property
    def edges_until_limit(self) -> float | None:
        """Edges expected to be scanned before LIMIT 5 is satisfied."""
        return DEFAULT_LIMIT / self.hit_rate if self.hit_rate else None

    @property
    def verdict(self) -> str:
        if not self.exists:
            return "FAST"
        if self.synthon_hits >= CONFIDENT_HITS:
            return "FAST"
        if self.estimated_edges < CHEAP_TO_EXHAUST_EDGES:
            return "FAST"
        return "SLOW"

    @property
    def reason(self) -> str:
        if not self.exists:
            return "no F2 node with this smiles - returns 0 rows immediately"
        if self.synthon_hits >= CONFIDENT_HITS:
            return (
                f"synthon is common ({self.synthon_hits} hits in"
                f" {self.sampled_edges:,} edges) so LIMIT {DEFAULT_LIMIT} exits"
                f" after ~{self.edges_until_limit:,.0f} edges"
            )
        if self.estimated_edges < CHEAP_TO_EXHAUST_EDGES:
            return (
                f"synthon is rare, but the neighbourhood is only"
                f" ~{self.estimated_edges:,.0f} edges - cheap to exhaust"
            )
        return (
            f"synthon absent from {self.sampled_edges:,} sampled edges, so all"
            f" ~{self.estimated_edges:,.0f} edges must be scanned"
        )


def predict(
    driver: Driver,
    pair: Pair,
    database: str | None = None,
    timeout: float | None = 180,
) -> Prediction:
    """Estimate a pair's cost without running its query."""
    started = time.perf_counter()
    with driver.session(database=database) as session:
        with session.begin_transaction(timeout=timeout) as tx:
            degrees = tx.run(
                DEGREE_PROBE, smiles=pair.smiles, neighbours=PROBE_NEIGHBOURS
            ).single()
            if degrees is None or not degrees["sampled"]:
                return Prediction(
                    pair=pair,
                    exists=False,
                    probe_seconds=time.perf_counter() - started,
                )
            synthons = tx.run(
                SYNTHON_PROBE,
                smiles=pair.smiles,
                synthon=pair.synthon,
                neighbours=PROBE_NEIGHBOURS,
                edges=PROBE_EDGES,
            ).single()
    return Prediction(
        pair=pair,
        exists=True,
        level1_degree=degrees["level1"],
        avg_level2_degree=degrees["avg_degree"] or 0.0,
        sampled_edges=synthons["sampled"],
        synthon_hits=synthons["hits"],
        probe_seconds=time.perf_counter() - started,
    )


@dataclass
class PairResult:
    """The outcome of running one pair's query ``repeats`` times."""

    pair: Pair
    durations: list[float] = field(default_factory=list)
    records: list[dict] = field(default_factory=list)
    error: str | None = None

    @property
    def average_seconds(self) -> float | None:
        return statistics.mean(self.durations) if self.durations else None


def run_query(
    driver: Driver,
    pair: Pair,
    num_hops: int = DEFAULT_NUM_HOPS,
    limit: int = DEFAULT_LIMIT,
    timeout: float | None = None,
    database: str | None = None,
) -> tuple[list[dict], float]:
    """Run one query, returning its records and how long the server took.

    ``timeout`` (seconds) is applied as a neo4j transaction timeout, so a
    problematic pair is abandoned by the server rather than left running.
    """
    query = build_query(pair.smiles, pair.synthon, num_hops, limit)
    with driver.session(database=database) as session:
        started = time.perf_counter()
        with session.begin_transaction(timeout=timeout) as tx:
            records = [record.data() for record in tx.run(query)]
        return records, time.perf_counter() - started


def time_pairs(
    driver: Driver,
    pairs: tuple[Pair, ...],
    num_hops: int = DEFAULT_NUM_HOPS,
    limit: int = DEFAULT_LIMIT,
    repeats: int = DEFAULT_REPEATS,
    timeout: float | None = None,
    database: str | None = None,
    report=print,
) -> list[PairResult]:
    """Run every pair ``repeats`` times, reporting progress as we go."""
    results: list[PairResult] = []
    for pair in pairs:
        result = PairResult(pair=pair)
        for _ in range(repeats):
            try:
                records, duration = run_query(
                    driver, pair, num_hops, limit, timeout, database
                )
            except Neo4jError as error:
                result.error = f"{type(error).__name__}: {error}"
                break
            result.durations.append(duration)
            result.records = records
        results.append(result)
        report(_format_result(result))
    return results


def _format_result(result: PairResult) -> str:
    pair = result.pair
    prefix = f"{pair.smiles} + {pair.synthon}"
    if result.error:
        return f"FAIL  {prefix}: {result.error}"
    average = result.average_seconds
    return (
        f"OK    {prefix}: {average:.6f} s average over {len(result.durations)}"
        f" queries, {len(result.records)} rows"
    )


def connect(uri: str, user: str | None, password: str | None) -> Driver:
    """Open a driver, using no authentication when no password is given."""
    auth = (user, password) if password is not None else None
    driver = GraphDatabase.driver(uri, auth=auth)
    driver.verify_connectivity()
    return driver


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--set",
        dest="pair_set",
        choices=sorted(PAIR_SETS),
        default="fast",
        help="which README.md pairs to run (default: fast)",
    )
    parser.add_argument("--uri", default=os.environ.get("NEO4J_URI", DEFAULT_URI))
    parser.add_argument("--user", default=os.environ.get("NEO4J_USERNAME", "neo4j"))
    parser.add_argument("--password", default=os.environ.get("NEO4J_PASSWORD"))
    parser.add_argument("--database", default=os.environ.get("NEO4J_DATABASE"))
    parser.add_argument("--num-hops", type=int, default=DEFAULT_NUM_HOPS)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument(
        "--predict",
        action="store_true",
        help="estimate each pair's cost with bounded probes instead of running it",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="server-side transaction timeout in seconds (advised for --set slow)",
    )
    args = parser.parse_args(argv)

    pairs = PAIR_SETS[args.pair_set]

    if args.predict:
        print(f"Predicting {len(pairs)} '{args.pair_set}' pairs against {args.uri}\n")
        with connect(args.uri, args.user, args.password) as driver:
            for pair in pairs:
                prediction = predict(driver, pair, database=args.database)
                print(
                    f"{prediction.verdict:5} {pair.smiles} + {pair.synthon}"
                    f"  ({prediction.probe_seconds:.2f}s probe)"
                )
                print(f"      {prediction.reason}")
        return 0

    print(f"Running {len(pairs)} '{args.pair_set}' pairs against {args.uri}")
    print(f"num_hops={args.num_hops} limit={args.limit} repeats={args.repeats}\n")

    with connect(args.uri, args.user, args.password) as driver:
        results = time_pairs(
            driver,
            pairs,
            num_hops=args.num_hops,
            limit=args.limit,
            repeats=args.repeats,
            timeout=args.timeout,
            database=args.database,
        )

    failed = [result for result in results if result.error]
    timed = [result for result in results if result.average_seconds is not None]
    print(
        f"\n{len(timed)} of {len(results)} pairs completed"
        + (f", slowest {max(r.average_seconds for r in timed):.6f} s" if timed else "")
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

# Cypher query tester (m2ms-2112)

Times the neo4j cypher query in [`cypher.map`](cypher.map) against a set of
SMILES/synthon pairs, and — more usefully — **predicts whether a pair will be
fast or slow before you run it**.

Some pairs return in under a tenth of a second. Others have never been observed
to finish: four of the ten pairs below exceeded a 300 second timeout, and one was
still running when it was cancelled after 13 minutes. The prediction mode tells
the two apart — usually in a fraction of a second — so you don't have to find out
the hard way.

## Quick start

You need [uv](https://docs.astral.sh/uv/) and a port-forward to the graph:

```bash
kubectl port-forward service/graph 7687 -n graph-x
```

Credentials come from the environment (or `--user` / `--password`):

```bash
export NEO4J_USERNAME=neo4j
export NEO4J_PASSWORD=...
```

Then predict the cost of the known pairs:

```bash
uv run cypher-query-tester --set all --predict
```

## Predicting cost

`--predict` runs two bounded probes per pair instead of the real query. Neither
probe can run away: degrees are read from the degree store in O(1), and the
synthon frequency comes from a `LIMIT`ed sample of the edges the query would
filter.

```
$ uv run cypher-query-tester --set all --predict
Predicting 10 'all' pairs against bolt://localhost:7687

FAST  CCNCC + [Xe]c1cscn1  (1.05s probe)
      synthon is common (86 hits in 13,782 edges) so LIMIT 5 exits after ~801 edges
FAST  CC(C)(CO)NC=O + [Xe]c1ncc1  (0.58s probe)
      synthon is rare, but the neighbourhood is only ~320,671 edges - cheap to exhaust
SLOW  c1cscn1 + CN(C)C(=O)CO[Xe]  (0.31s probe)
      synthon absent from 13,433 sampled edges, so all ~137,276,031 edges must be scanned
FAST  c1ccc2cnc2c1 + CN(C)C(=O)CO[Xe]  (0.06s probe)
      no F2 node with this smiles - returns 0 rows immediately
SLOW  CCC=O + O=c1cccnn1[Xe]  (0.77s probe)
      synthon absent from 10,997 sampled edges, so all ~66,600,634 edges must be scanned
```

That run covered all ten pairs in under six seconds, against a worst case of
100+ minutes for a single one of them. All ten verdicts match what was measured
by actually running the queries.

The verdict follows three checks, cheapest first:

1. **No `F2` node for the SMILES?** Fast — it returns zero rows immediately.
2. **Synthon seen 5+ times in the sample?** Fast — `LIMIT 5` will exit early,
   after roughly `5 / hit_rate` edges.
3. **Otherwise, neighbourhood under 5 million edges?** Fast anyway, cheap to
   exhaust. Above that, **slow**.

A pair is only slow when *both* factors are bad: a synthon too rare to trigger
the early exit **and** a neighbourhood too large to exhaust. Either alone is
harmless.

Use it as an admission check — probe first, and defer or refuse anything that
comes back `SLOW` rather than letting it hold a connection for an hour.

### Tuning the probes

Probe cost is usually well under a second, but it is not guaranteed: probing a
large-degree start node on a cold cache has taken 32s and 74s, because the sample
has to fault pages in from disk. That is still far cheaper than the query it
avoids, but don't assume it is instant.

`PROBE_NEIGHBOURS` (2,000) and `PROBE_EDGES` (50,000) in
[`src/cypher_query.py`](src/cypher_query.py) control sample size. Lowering them
makes probes cheaper but less reliable; the thresholds below them were
calibrated at the current values, so re-check the verdicts if you change them.

## Timing the real queries

```bash
uv run cypher-query-tester --set fast                 # the 5 fast pairs
uv run cypher-query-tester --set slow --timeout 300   # ALWAYS use a timeout here
uv run cypher-query-tester --set all --repeats 3
```

| Flag | Default | Notes |
| --- | --- | --- |
| `--set {fast,slow,all}` | `fast` | Which pairs to run |
| `--predict` | off | Estimate instead of running |
| `--timeout SECONDS` | none | Server-side transaction timeout. **Essential for `--set slow`** |
| `--repeats N` | 3 | Runs per pair; the average is reported |
| `--num-hops N` | 2 | `maxLevel` of the expansion |
| `--limit N` | 5 | Result `LIMIT` |
| `--uri` / `--user` / `--password` / `--database` | `bolt://localhost:7687`, `neo4j`, — , — | Also read from `NEO4J_URI`, `NEO4J_USERNAME`, `NEO4J_PASSWORD`, `NEO4J_DATABASE` |

Exit status is 1 if any pair failed, 0 otherwise.

> **Always pass `--timeout` when running the slow set.** Without it a
> problematic pair will hold a server connection indefinitely — one was still
> running after 100 minutes. With it, the server abandons the transaction and
> the run moves on.

## The query

[`cypher.map`](cypher.map) is a python `%`-format template taking four
variables: `smiles`, `synthon`, `num_hops` and `limit`. It expands from an `F2`
node along `FRAG` relationships up to `num_hops` levels, then returns the `Mol`
nodes reached by a `FRAG` edge carrying the given synthon.

The template is duplicated as `CYPHER_TEMPLATE` in `src/cypher_query.py` so that
what gets timed is exactly what an application would send.

## The pairs

These are hardcoded as `FAST_PAIRS` and `SLOW_PAIRS` in `src/cypher_query.py`.
**If you edit one, edit the other** — nothing keeps them in sync automatically.

Expected to be fast:

| SMILES | Synthon |
| --- | --- |
| `CCNCC` | `[Xe]c1cscn1` |
| `CC(C)(CO)NC=O` | `[Xe]n1cnc2ccccc21` |
| `CC(C)(CO)NC=O` | `[Xe]c1ncc1` |
| `CCC=O` | `[Xe]N1CCCC1` |
| `C1CCCCC1` | `[Xe]N1CCNCC1` |

Expected to be problematic:

| SMILES | Synthon | |
| --- | --- | --- |
| `c1cscn1` | `CN(C)C(=O)CO[Xe]` | confirmed slow |
| `c1ccc2ncccc2c1` | `[Xe]N1CCCc2ccccc21` | confirmed slow |
| `c1ccc2cnc2c1` | `CN(C)C(=O)CO[Xe]` | **see below** |
| `CCC=O` | `O=c1cccnn1[Xe]` | confirmed slow |
| `c1cscn1` | `O=C1CC([Xe])=NN1` | confirmed slow |

## Gotchas

**`c1ccc2cnc2c1` does not exist in the graph.** There is no `F2` node with that
SMILES, so the pair returns 0 rows in ~0.1s and cannot demonstrate the problem
it was listed for. The string also describes a four-membered aromatic ring,
which is chemically implausible — it looks like a transcription error and should
be checked at source.

**Timings depend heavily on cache state.** The same query has been measured at
93.6s cold and 0.5s warm, and at 195.4s cold and 1.7s warm. The page cache holds
about 3.5% of the store, so a first touch pays for disk. Record whether a
measurement was cold or warm, or it isn't comparable to anything.

**`CCC=O` appears in both lists.** It is fast with `[Xe]N1CCCC1` and has never
finished with `O=c1cccnn1[Xe]`. The starting molecule alone tells you nothing —
the synthon decides.

## Findings

Two write-ups cover the investigation behind this tool:

- **The Synthon Cost Cliff** — why four pairs never finish, and how the
  prediction rule was derived and validated.
- **Fixing the Cost Cliff** — whether a query rewrite, indexes, or more memory
  would help, with measurements for each.

Headline conclusions:

- The cost is set by synthon selectivity multiplied by neighbourhood size, not
  by the starting molecule's fan-out.
- The workload is I/O-bound. The page cache covers 3.5% of a 4,806 GB store.
- Deduplicating the expansion (`DISTINCT`, or `NODE_GLOBAL` uniqueness) is worth
  about 1.5–1.75× at small scale and nothing at large scale.
- **The query as written returns duplicate rows.** `RELATIONSHIP_GLOBAL`
  deduplicates relationships rather than end nodes, so `LIMIT 5` can return
  three distinct molecules. Adding `DISTINCT` fixes it.

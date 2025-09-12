import mrich
import logging
from mrich import print

import json
import time
import asyncio
from rdkit.Chem import MolFromSmiles
from neo4j import GraphDatabase, AsyncGraphDatabase

from .config import CONFIG
from .tools import load_sig_factory, calc_pharm_fp


def check_config():
    graph_vars = ["GRAPH_LOCATION", "GRAPH_USERNAME", "GRAPH_PASSWORD"]
    missing = []
    for var in graph_vars:
        if var not in CONFIG:
            missing.append(var)
    if missing:
        mrich.error("Configuration missing:", missing)
        raise ValueError(f"Configuration missing: {missing}")


def get_driver():
    check_config()
    return GraphDatabase.driver(
        CONFIG["GRAPH_LOCATION"],
        auth=(CONFIG["GRAPH_USERNAME"], CONFIG["GRAPH_PASSWORD"]),
    )


async def aget_driver():
    check_config()
    return AsyncGraphDatabase.driver(
        CONFIG["GRAPH_LOCATION"],
        auth=(CONFIG["GRAPH_USERNAME"], CONFIG["GRAPH_PASSWORD"]),
    )


async def arun_query(query, **kwargs):
    driver = await aget_driver()
    async with driver:
        async with driver.session() as session:
            result = await session.run(query, **kwargs)
            records = [record async for record in result]
            return records


def run_query(query, **kwargs):
    driver = get_driver()
    with driver:
        with driver.session() as session:
            result = session.run(query, **kwargs)
            records = [record for record in result]
            return records


async def aget_subnodes(
    smiles: str,
    terminal_nodes: bool = CONFIG["FRAGMENT_TERMINAL_SUBNODES"],
    progress=None,
    task=None,
):
    """
    Get subnodes for a given node (retrieve using SMILES)

    :param smiles: SMILES string for node to retrieve subnodes
    :param terminal_subnodes: whether to only return 'terminal' subnodes (can't be broken down further)
    :return: list of unique subnode SMILES
    """

    if terminal_nodes:
        query = """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2)
        WHERE NOT ()-[:FRAG]-(f)-[:FRAG]->()
        RETURN f
        """
    else:
        query = """
        MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2) 
        RETURN f
        """

    records = await arun_query(query, smiles=smiles)
    subnodes = [record["f"]["smiles"] for record in records]

    if progress:
        progress.update(task, advance=1)

    return set(subnodes)


async def aget_synthons(
    smiles: str,
    terminal_nodes: bool = CONFIG["FRAGMENT_TERMINAL_SYNTHONS"],
    progress=None,
    task=None,
):
    """
    Get constituent synthons (compounds added or removed during transformation) for a given node SMILES.
    [Xe] denotes the attachment point.

    :param smiles: SMILES string of node to retrieve synthons
    :param terminal_synthons: whether to return 'terminal' synthons, i.e. can't be broken down more
    :return: list of constituent synthon SMILES strings
    """

    if terminal_nodes:
        query = """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..15]->(b:F2)
        WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
        RETURN e[-1] as edge
        """
    else:
        query = """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..15]->(b:F2)
        RETURN e[-1] as edge
        """

    records = await arun_query(query, smiles=smiles)
    edges = [edge for record in records if (edge := record["edge"])]

    synthons = set()

    for edge in edges:
        for p in ["prop_synthon", "prop_core"]:
            i = edge[p]
            if not i:
                continue

            if i.count("Xe") != 1:
                continue

            synthons.add(i)

    if progress:
        progress.update(task, advance=1)

    return synthons


async def aget_r_groups(
    smiles: str,
    progress=None,
    task=None,
):

    query = """
    MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..15]->(b:F2)
    WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
    AND e[-1].prop_synthon contains '[Xe]'
    AND NOT e[-1].prop_synthon=e[-2].prop_synthon
    RETURN e[-1].prop_synthon as synthon, e[-2].prop_synthon as r_group;
    """

    records = await arun_query(query, smiles=smiles)

    results = []
    for record in records:
        results.append((record["synthon"], record["r_group"]))

    if progress:
        progress.update(task, advance=1)

    return results


def get_pure_expansions(
    smiles: str,
    synthon: str,
    num_hops: int = 2,
    limit: int = 5,
    index: int | None = None,
    cache_dir=None,
    cached_only=False,
):

    if cache_dir:
        cache_file = cache_dir / f"pure_{smiles}_{synthon}_{num_hops}_{limit}.json"
        if cache_file.exists():
            logging.info(f"Using cache {index} {smiles} {synthon}")
            return json.load(open(cache_file, "rt"))
        elif cached_only:
            return None

    query = """
    MATCH (a:F2 {smiles: $smiles})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
    WHERE e.prop_synthon=$synthon
    WITH c.smiles as smi, c.cmpd_ids as ids
    RETURN smi, ids
    """ % {
        "num_hops": num_hops
    }

    if limit:
        query = query + f" LIMIT {limit}"

    logging.info(f"Starting impure expansion {index} {smiles} {synthon}")

    try:
        records = run_query(query, smiles=smiles, synthon=synthon)
    except Exception as e:
        mrich.error(index, e)
        raise Exception(f"{smiles=} {synthon=} {e}")

    results = []
    for record in records:
        results.append((record["ids"], record["smi"]))

    if cache_dir:
        json.dump(results, open(cache_file, "wt"), indent=2)

    logging.info(f"Success {index} {smiles} {synthon} #results: {len(results)}")

    return results


def get_impure_expansions(
    smiles: str,
    synthon: str,
    num_hops: int = 2,
    limit: int = 5,
    index: int | None = None,
    cache_dir=None,
    cached_only=False,
):
    
    if cache_dir:
        cache_file = cache_dir / f"impure_{smiles}_{synthon}_{num_hops}_{limit}.json"
        if cache_file.exists():
            logging.info(f"Using cache {index} {smiles} {synthon}")
            return json.load(open(cache_file, "rt"))
        elif cached_only:
            return None
    
    logging.info(f"Starting impure expansion {index} {smiles} {synthon}")

    sig_factory = load_sig_factory(
        fdef_file=CONFIG["FINGERPRINT_FDEF"],
        max_point_count=CONFIG["FINGERPRINT_MAXPOINTCOUNT"],
        bins=json.loads(CONFIG["FINGERPRINT_BINS"]),
    )

    vector = calc_pharm_fp(MolFromSmiles(synthon), sig_factory, as_str=False)

    threshold = CONFIG["KNITWORK_SIMILARITY_THRESHOLD"]
    metric = CONFIG["KNITWORK_SIMILARITY_METRIC"]

    query = """
    MATCH (a:F2 {smiles: $smiles})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
    WHERE e.prop_pharmfp IS NOT NULL
    WITH usersimilarity.tanimoto_similarity(e.prop_pharmfp, $vector) as sim, c.smiles as smi, e.prop_synthon as syn, c.cmpd_ids as ids
    WHERE sim >= $threshold
    AND NOT e.prop_synthon=$query_synthon
    RETURN smi, syn, sim, ids
    """ % {
        "num_hops":num_hops,
    }

    if limit:
        query = query + f" LIMIT {limit}"

    try:
        records = run_query(query, 
            smiles=smiles, 
            query_synthon=synthon,
            vector=vector,
            threshold=threshold,
            metric=metric,
            # num_hops=num_hops,
        )
    except Exception as e:
        mrich.error(index, e)
        raise Exception(f"{smiles=} {synthon=} {e}")

    results = []
    for record in records:
        results.append(
            (
                record["smi"],  # expansion
                record["syn"],  # synthon
                record["sim"],  # similarity
                record["ids"],  # compound names / IDs
            )
        )

    if cache_dir:
        json.dump(results, open(cache_file, "wt"), indent=2)

    logging.info(f"Success {index} {smiles} {synthon} #results: {len(results)}")

    return results

from neo4j import GraphDatabase, AsyncGraphDatabase
from .config import CONFIG
import mrich
from mrich import print
import asyncio


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


async def run_query(query, **kwargs):
    driver = await aget_driver()
    async with driver:
        async with driver.session() as session:
            result = await session.run(query, **kwargs)
            records = [record async for record in result]
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

    records = await run_query(query, smiles=smiles)
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
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(b:F2)
        WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
        RETURN e[-1] as edge
        """
    else:
        query = """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(b:F2)
        RETURN e[-1] as edge
        """

    records = await run_query(query, smiles=smiles)
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

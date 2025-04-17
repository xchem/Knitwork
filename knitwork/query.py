
from neo4j import GraphDatabase
from .config import CONFIG

graph_vars = ["GRAPH_LOCATION", "GRAPH_USERNAME", "GRAPH_PASSWORD"]
missing = []
for var in graph_vars:
    if var not in CONFIG:
        missing.append(var)
if missing:
    mrich.error("Configuration missing:", missing)
    raise ValueError(f"Configuration missing: {missing}")

DRIVER = GraphDatabase.driver(
    CONFIG["GRAPH_LOCATION"], 
    auth=(CONFIG["GRAPH_USERNAME"], CONFIG["GRAPH_PASSWORD"]),
)

def get_subnodes(smiles: str, terminal_subnodes: bool = CONFIG["FRAGMENT_TERMINAL_SUBNODES"]) -> list:
    """
    Get subnodes for a given node (retrieve using SMILES)

    :param smiles: SMILES string for node to retrieve subnodes
    :param terminal_subnodes: whether to only return 'terminal' subnodes (can't be broken down further)
    :return: list of unique subnode SMILES
    """
    terminal_subnode_query = (
        """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2)
        WHERE NOT ()-[:FRAG]-(f)-[:FRAG]->()
        RETURN f
        """
    )
    subnode_query = (
        "MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2) RETURN f"
    )
    if terminal_subnodes:
        query = terminal_subnode_query
    else:
        query = subnode_query
        
    with DRIVER.session() as session:
        subnodes = [record['f']['smiles'] for record in session.run(query, smiles=smiles)]
    return list(set(subnodes))

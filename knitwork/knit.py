import pandas as pd
from mrich import print
import mrich
from .query import aget_pure_expansions
from .config import CONFIG, print_config
import time
from rich.progress import Progress
import asyncio

# import sys
# import json
from pathlib import Path
from rdkit.Chem import MolFromSmiles, PandasTools


def pure_merge(
    pairs_df: "pd.DataFrame",
    output_dir: str = "knitwork_output",
):

    mrich.h2("knitwork.knit.pure_merge()")
    print_config("GRAPH_LOCATION")
    print_config("KNITWORK")

    output_dir = Path(output_dir)
    if not output_dir.exists():
        mrich.writing(output_dir)
        output_dir.mkdir(parents=True)

    cache_dir = output_dir / "cache"
    mrich.var("cache_dir", cache_dir)
    if not cache_dir.exists():
        mrich.writing(cache_dir)
        cache_dir.mkdir(parents=True)

    substructure_pairs = get_unique_substructure_pairs(pairs_df)

    # asynchronous expansion queries
    with Progress() as progress:
        n_unique = len(substructure_pairs)
        task = progress.add_task("query pure_expansions", total=n_unique)
        results = asyncio.run(
            pure_merge_task(substructure_pairs, cache_dir, progress, task)
        )

    if not results:
        mrich.error("No results")
        return None

    # process results
    data = []
    for ((hit1, hit2), subnode, synthon), result in zip(substructure_pairs, results):
        for names, merge in result:
            data.append(
                dict(
                    ID_A=hit1,
                    ID_B=hit2,
                    subnode_A=subnode,
                    synthon_B=synthon,
                    merge_smiles=merge,
                    catalogue_names=names,
                )
            )

    mrich.var("#merges", len(data))

    df = pd.DataFrame(data)
    df.loc[:, "ROMol"] = df["merge_smiles"].apply(MolFromSmiles)

    # write pickle
    df_path = output_dir / "pure_merges.pkl.gz"
    mrich.writing(df_path)
    df.to_pickle(df_path)

    df.loc[:, "ID"] = df.index

    # write SDF
    sdf_path = output_dir / "pure_merges.sdf"
    mrich.writing(sdf_path)
    PandasTools.WriteSDF(
        df,
        str(sdf_path),
        molColName="ROMol",
        idName="ID",
        properties=df.columns,
    )

    return df


def get_unique_substructure_pairs(df):

    # get unique substructure pairs

    substructure_pairs = set()
    for i, row in df.iterrows():
        for subnode_A in row["subnodes_A"]:
            for synthon_B in row["synthons_B"]:
                substructure_pairs.add((i, subnode_A, synthon_B))

    mrich.var("#unique substructure pairs", len(substructure_pairs))

    return substructure_pairs


semaphore = asyncio.Semaphore(CONFIG["KNITWORK_NUM_CONNECTIONS"])


async def limited_aget_pure_expansions(*args, **kwargs):
    async with semaphore:
        return await aget_pure_expansions(*args, **kwargs)


async def pure_merge_task(
    substructure_pairs, cache_dir=None, progress=None, prog_task=None
):

    results = None

    try:
        async with asyncio.TaskGroup() as tg:
            tasks = [
                tg.create_task(
                    aget_pure_expansions(
                        smiles,
                        synthon,
                        cache_dir=cache_dir,
                        progress=progress,
                        task=prog_task,
                    )
                )
                for _, smiles, synthon in substructure_pairs
            ]

        results = [task.result() for task in tasks]

    except* Exception as eg:
        for i, e in enumerate(eg.exceptions):
            mrich.error(f"[{i}] Task failed: {e}")

    return results

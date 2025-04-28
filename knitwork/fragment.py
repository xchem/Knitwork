import mrich
from mrich import print
from pathlib import Path
from .config import CONFIG, print_config
import asyncio
from .query import aget_subnodes, aget_synthons, aget_r_groups
from rich.progress import Progress

MOL_CACHE = {}
MOL_C = None

def fragment(
    input_sdf: Path,
    output_dir: Path | str = "fragment_output",
    overlap_cutoff: float = CONFIG["FRAGMENT_OVERLAP_CUTOFF"],
    distance_cutoff: float = CONFIG["FRAGMENT_DISTANCE_CUTOFF"],
):

    import pandas as pd
    from rdkit.Chem import PandasTools, MolToSmiles, MolFromSmarts
    from itertools import permutations
    from .tools import pair_overlap, pair_min_distance

    global MOL_C
    MOL_C = MolFromSmarts("[#6]")

    mrich.h2("knitwork.fragment.fragment()")

    input_sdf = Path(input_sdf)
    output_dir = Path(output_dir)

    # print config
    mrich.var("input_sdf", input_sdf)
    mrich.var("output_dir", output_dir)
    print_config("GRAPH_LOCATION")
    print_config("FRAGMENT")

    # check paths
    assert input_sdf.exists()
    if not output_dir.exists():
        mrich.writing(output_dir)
        output_dir.mkdir(parents=True)

    # get mols
    mol_df = PandasTools.LoadSDF(str(input_sdf.resolve()))
    mol_df = mol_df[["ID", "ROMol"]]
    mol_df["smiles"] = mol_df.apply(lambda x: MolToSmiles(x.ROMol), axis=1)
    MOL_CACHE.update({r["smiles"]: r["ROMol"] for i, r in mol_df.iterrows()})
    n_molecules = len(mol_df)
    mrich.var("#molecules", n_molecules)

    # asynchronous mol tasks
    with Progress() as progress:
        smiles_list = mol_df["smiles"].unique()
        n_unique = len(smiles_list)
        mrich.var("#unique smiles", n_unique)
        t1 = progress.add_task("query subnodes", total=n_unique)
        t2 = progress.add_task("query synthons", total=n_unique)
        t3 = progress.add_task("query r_groups", total=n_unique)
        results = asyncio.run(fragment_tasks(smiles_list, progress, (t1, t2, t3)))

    # filter results
    for smiles, v in results.items():
        v["subnodes"] = filter_smiles_list(v["subnodes"], synthons=False)
        v["synthons"] = filter_smiles_list(v["synthons"], synthons=True)

    # update molecule dataframe
    mol_df["subnodes"] = mol_df["smiles"].map(lambda s: results[s]["subnodes"])
    mol_df["synthons"] = mol_df["smiles"].map(lambda s: results[s]["synthons"])
    mol_df["r_groups"] = mol_df["smiles"].map(lambda s: results[s]["r_groups"])

    # write mol_df
    mol_df_path = output_dir / "molecules.pkl.gz"
    mrich.writing(mol_df_path)
    mol_df.to_pickle(mol_df_path)
    mol_sdf_path = output_dir / "molecules.sdf"
    mrich.writing(mol_sdf_path)
    PandasTools.WriteSDF(
        mol_df,
        str(mol_sdf_path),
        molColName="ROMol",
        idName="ID",
        properties=mol_df.columns,
    )

    # get pairs
    pair_df = mol_df.copy()
    pair_df["key"] = 1
    pair_df = pair_df.merge(pair_df, on="key", suffixes=["_A", "_B"])
    pair_df = pair_df.drop(columns="key")
    pair_df = pair_df[pair_df["ID_A"] != pair_df["ID_B"]]
    pair_df = pair_df.set_index(["ID_A", "ID_B"])
    mrich.var("#pairs", len(pair_df))

    # filter by overlap
    pair_df["overlap"] = pair_df.apply(
        lambda x: pair_overlap(x["ROMol_A"], x["ROMol_B"]), axis=1
    )
    overlapping = pair_df["overlap"] > overlap_cutoff
    mrich.var(f"#(overlap > {overlap_cutoff})", len(pair_df[overlapping]), "pairs")
    pair_df = pair_df[~overlapping]

    # filter by distance
    pair_df["distance"] = pair_df.apply(
        lambda x: pair_min_distance(x["ROMol_A"], x["ROMol_B"]), axis=1
    )
    distant = pair_df["distance"] > distance_cutoff
    mrich.var(f"#(distance > {distance_cutoff})", len(pair_df[distant]), "pairs")
    pair_df = pair_df[~distant]

    mrich.var(f"#pairs (post-filter)", len(pair_df), "pairs")

    # write pair_df
    pair_df_path = output_dir / "pairs.pkl.gz"
    mrich.writing(pair_df_path)
    pair_df.to_pickle(pair_df_path)


async def fragment_tasks(smiles_list, progress, tasks):

    t1, t2, t3 = tasks

    tasks = [
        asyncio.gather(
            aget_subnodes(smiles, progress=progress, task=t1),
            aget_synthons(smiles, progress=progress, task=t2),
            aget_r_groups(smiles, progress=progress, task=t3),
        )
        for smiles in smiles_list
    ]

    results = await asyncio.gather(*tasks)

    return {
        smiles: {"subnodes": subnodes, "synthons": synthons, "r_groups":r_groups}
        for smiles, (subnodes, synthons, r_groups) in zip(smiles_list, results)
    }


def filter_smiles_list(smiles_list, synthons: bool):

    from rdkit import Chem

    filtered = [s for s in smiles_list]
    if CONFIG["FRAGMENT_CHECK_SINGLE_MOL"]:
        filtered = [s for s in filtered if "." not in s]

    if CONFIG["FRAGMENT_CHECK_CARBONS"]:
        n_c = CONFIG["FRAGMENT_MIN_CARBONS"]
        filtered = [
            s
            for s in filtered
            if len(
                MOL_CACHE.setdefault(s, Chem.MolFromSmiles(s)).GetSubstructMatches(MOL_C)
            )
            >= n_c
        ]

    if CONFIG["FRAGMENT_CHECK_CARBON_RING"]:

        new_filtered = []
        for s in filtered:
            mol = MOL_CACHE.setdefault(s, Chem.MolFromSmiles(s))

            num_rings = mol.GetRingInfo().NumRings()
            if num_rings != 1:
                new_filtered.append(s)

            num_ring_atoms = sum([atom.IsInRing() for atom in mol.GetAtoms()])
            num_carbons = len(mol.GetSubstructMatches(MOL_C))
            num_atoms = mol.GetNumAtoms()

            # one atom will be a xenon (attachment point)
            if synthons and (
                num_ring_atoms != (num_atoms - 1) or num_carbons != (num_atoms - 1)
            ):
                new_filtered.append(s)

            elif not (num_ring_atoms != num_atoms or num_carbons != num_atoms):
                new_filtered.append(s)

        filtered = new_filtered

    return filtered

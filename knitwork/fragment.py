import mrich
from mrich import print
from pathlib import Path
from .config import CONFIG


def fragment(
    input_sdf: Path,
    output_dir: Path | str = "fragment_output",
    overlap_cutoff: float = CONFIG["FRAGMENT_OVERLAP_CUTOFF"],
    distance_cutoff: float = CONFIG["FRAGMENT_DISTANCE_CUTOFF"],
):

    from .query import get_subnodes
    import pandas as pd
    from rdkit.Chem import PandasTools, MolToSmiles
    from itertools import permutations
    from .tools import pair_overlap, pair_min_distance

    mrich.h2("knitwork.fragment.fragment()")

    input_sdf = Path(input_sdf)
    output_dir = Path(output_dir)

    mrich.var("input_sdf", input_sdf)
    mrich.var("output_dir", output_dir)

    assert input_sdf.exists()
    if not output_dir.exists():
        mrich.writing(output_dir)
        output_dir.mkdir(parents=True)

    # get mols
    mol_df = PandasTools.LoadSDF(str(input_sdf.resolve()))
    mol_df = mol_df[["ID", "ROMol"]]
    mol_df["smiles"] = mol_df.apply(lambda x: MolToSmiles(x.ROMol), axis=1)
    mrich.var("#molecules", len(mol_df))

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

    # get subnodes

    for smiles in mol_df["smiles"].values:
        print(smiles)
        subnodes = get_subnodes(smiles)
        print(subnodes)
        # break

    return pair_df

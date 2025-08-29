import mrich
from typer import Typer
from pathlib import Path
import json

app = Typer()


@app.command()
def fragment(
    input_sdf: str,
    output_dir: str = "fragment_output",
):
    mrich.h1("FRAGMENT")
    from .fragment import fragment
    from rdkit.Chem import PandasTools

    input_sdf = Path(input_sdf)
    mrich.var("input_sdf", input_sdf)
    mol_df = PandasTools.LoadSDF(str(input_sdf.resolve()))

    fragment(mol_df, output_dir)


@app.command()
def pure_merge(
    # pairs_df: str,
    fragment_dir: str = "fragment_output",
    output_dir: str = "knitwork_output",
    cached_only: bool = False,
):
    mrich.h1("PURE MERGE")
    from .knit import pure_merge
    import pandas as pd

    fragment_dir = Path(fragment_dir)
    mrich.var("fragment_dir", fragment_dir)
    mrich.var("output_dir", output_dir)
    assert fragment_dir.exists()
    assert fragment_dir.is_dir()

    pairs_df = fragment_dir / "pairs.pkl.gz"
    mrich.var("pairs_df", pairs_df)
    pairs_df = pd.read_pickle(pairs_df)

    pure_merge(pairs_df=pairs_df, output_dir=output_dir, cached_only=cached_only)


@app.command()
def impure_merge():
    mrich.h1("IMPURE MERGE")
    raise NotImplementedError


@app.command()
def configure(
    var: str,
    value: str,
):

    from .config import VARIABLES, CONFIG, dump_config

    if var not in VARIABLES:
        mrich.error("Variable is not configurable:", var)
        mrich.var("accepted variables", VARIABLES, separator=":")
        raise ValueError(f"'{var}' is not configurable")

    t = VARIABLES[var]

    mrich.var(var, value)

    if t == bool:
        if value == "False":
            value = False
        elif value == "True":
            value = True
        else:
            raise ValueError("value must be 'True' or 'False'")
    elif t != str:
        value = t(value)

    mrich.var(var, value)

    CONFIG[var] = value
    dump_config(CONFIG)


if __name__ == "__main__":
    app()

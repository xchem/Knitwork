import mrich
from typer import Typer
from pathlib import Path
import json

app = Typer()


@app.command()
def fragment(
    input_sdf: str,
    output_dir: str = "fragment_output",
    config_path: str = None,
):
    """Fragment and pair up input molecules so that substructure matching can be run"""

    mrich.h1("FRAGMENT")
    
    init_config(config_path=config_path)

    from .fragment import fragment
    from rdkit.Chem import PandasTools

    input_sdf = Path(input_sdf)
    mrich.var("input_sdf", input_sdf)
    mol_df = PandasTools.LoadSDF(str(input_sdf.resolve()))

    fragment(mol_df, output_dir)


@app.command()
def pure_merge(
    fragment_dir: str = "fragment_output",
    output_dir: str = "knitwork_output",
    cached_only: bool = False,
    limit: int = 5,
    config_path: str = None,
):
    """Enumerate 'pure' knitwork merges"""

    mrich.h1("PURE MERGE")

    init_config(config_path=config_path)

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

    pure_merge(
        pairs_df=pairs_df, output_dir=output_dir, cached_only=cached_only, limit=limit
    )


@app.command()
def impure_merge(
    fragment_dir: str = "fragment_output",
    output_dir: str = "knitwork_output",
    cached_only: bool = False,
    limit: int = 5,
    config_path: str = None,
):
    """Enumerate 'impure' knitwork merges"""

    mrich.h1("IMPURE MERGE")
    
    init_config(config_path=config_path)
    
    from .knit import impure_merge
    import pandas as pd

    fragment_dir = Path(fragment_dir)
    mrich.var("fragment_dir", fragment_dir)
    mrich.var("output_dir", output_dir)
    assert fragment_dir.exists()
    assert fragment_dir.is_dir()

    pairs_df = fragment_dir / "pairs.pkl.gz"
    mrich.var("pairs_df", pairs_df)
    pairs_df = pd.read_pickle(pairs_df)

    impure_merge(
        pairs_df=pairs_df, output_dir=output_dir, cached_only=cached_only, limit=limit
    )


@app.command()
def configure(
    var: str,
    value: str,
    config_path: str = None,
):
    """Set a configuration variable"""
    
    init_config(config_path=config_path)

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
    dump_config(CONFIG, config_path=config_path)


@app.command()
def combine_inputs(
    inputs: list[str],
    output: str,
):
    """Combine SDF inputs into a single file"""

    mrich.var("#inputs", len(inputs))
    mrich.var("inputs", inputs)
    mrich.var("output", output)

    from rdkit import Chem

    mols = []
    for sdf_file in inputs:
        suppl = Chem.SDMolSupplier(sdf_file)
        mols.extend([m for m in suppl if m is not None])

    mrich.writing(output)
    writer = Chem.SDWriter(output)
    for mol in mols:
        writer.write(mol)
    writer.close()

def init_config(
    config_path: str | Path | None = None,
) -> None:
    from .config import setup_config
    setup_config(config_path=config_path)

if __name__ == "__main__":
    app()

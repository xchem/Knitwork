import mrich
from typer import Typer
from pathlib import Path

app = Typer()


@app.command()
def fragment(
    input_sdf: str,
    output_dir: str = "fragment_output",
):
    mrich.h1("FRAGMENT")
    from .fragment import fragment

    fragment(input_sdf, output_dir)


@app.command()
def knit():
    mrich.h1("KNIT")
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

    if value == "False":
        value = False
    elif value == "True":
        value = True

    mrich.var(var, value)

    CONFIG[var] = value
    dump_config(CONFIG)


if __name__ == "__main__":
    app()

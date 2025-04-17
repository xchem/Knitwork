import json
import mrich
from pathlib import Path

CONFIG = None
CONFIG_PATH = (Path(__file__).parent.parent / "config.json").resolve()

VARIABLES = [
    "GRAPH_LOCATION",
    "GRAPH_USERNAME",
    "GRAPH_PASSWORD",
    "FRAGMENT_OVERLAP_CUTOFF",
    "FRAGMENT_DISTANCE_CUTOFF",
    "FRAGMENT_TERMINAL_SYNTHONS",
    "FRAGMENT_TERMINAL_SUBNODES",
]

DEFAULTS = {
    "FRAGMENT_TERMINAL_SYNTHONS": True,
    "FRAGMENT_TERMINAL_SUBNODES": True,
    "FRAGMENT_OVERLAP_CUTOFF": 0.56,
    "FRAGMENT_DISTANCE_CUTOFF": 5.0,
}


def load_config():

    if CONFIG_PATH.exists():
        # mrich.reading(CONFIG_PATH)
        return json.load(open(CONFIG_PATH, "rt"))
    else:
        config = DEFAULTS.copy()
        dump_config(config)
        return config


def dump_config(config):
    mrich.writing(CONFIG_PATH)
    json.dump(config, open(CONFIG_PATH, "wt"), indent=2)


def setup_config():
    global CONFIG
    CONFIG = load_config()


# if __name__ == "__main__":
setup_config()

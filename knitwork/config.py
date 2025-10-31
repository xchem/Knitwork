import json
import mrich
from pathlib import Path

CONFIG = None
CONFIG_PATH = (Path(__file__).parent.parent / "config.json").resolve()

VARIABLES = {
    "GRAPH_LOCATION": str,
    "GRAPH_USERNAME": str,
    "GRAPH_PASSWORD": str,
    "FRAGMENT_OVERLAP_CUTOFF": float,
    "FRAGMENT_DISTANCE_CUTOFF": float,
    "FRAGMENT_TERMINAL_SYNTHONS": bool,
    "FRAGMENT_TERMINAL_SUBNODES": bool,
    "FRAGMENT_CHECK_SINGLE_MOL": bool,
    "FRAGMENT_CHECK_CARBONS": bool,
    "FRAGMENT_CHECK_CARBON_RING": bool,
    "KNITWORK_NUM_CONNECTIONS": int,
    "KNITWORK_SIMILARITY_THRESHOLD": float,
    "KNITWORK_SIMILARITY_METRIC": str,
    "FINGERPRINT_FDEF": str,
    "FINGERPRINT_MAXPOINTCOUNT": int,
    "FINGERPRINT_BINS": str,
}

DEFAULTS = {
    "FRAGMENT_TERMINAL_SYNTHONS": True,
    "FRAGMENT_TERMINAL_SUBNODES": True,
    "FRAGMENT_OVERLAP_CUTOFF": 0.56,
    "FRAGMENT_DISTANCE_CUTOFF": 5.0,
    "FRAGMENT_CHECK_SINGLE_MOL": True,
    "FRAGMENT_CHECK_CARBONS": True,
    "FRAGMENT_CHECK_CARBON_RING": True,
    "FRAGMENT_MIN_CARBONS": 3,
    "KNITWORK_NUM_CONNECTIONS": 4,
    "KNITWORK_SIMILARITY_THRESHOLD": 0.9,
    "KNITWORK_SIMILARITY_METRIC": "usersimilarity.tanimoto_similarity",
    "FINGERPRINT_FDEF": "FeatureswAliphaticXenon.fdef",
    "FINGERPRINT_MAXPOINTCOUNT": 2,
    "FINGERPRINT_BINS": "[[0, 2], [2, 5], [5, 8]]",
}


def load_config(
    config_path: str | Path | None = None,
) -> dict:
    """Load configuration from JSON, or use defaults"""
    
    config_path = config_path or CONFIG_PATH

    if config_path.exists():
        return json.load(open(config_path, "rt"))
    else:
        config = DEFAULTS.copy()
        dump_config(config, config_path=config_path)
        return config


def dump_config(
    config: dict, 
    config_path: str | Path | None = None,
) -> None:
    """Dump configuration as JSON"""

    config_path = config_path or CONFIG_PATH
    mrich.writing(config_path)
    json.dump(config, open(config_path, "wt"), indent=2)


def setup_config(
    config_path: str | Path | None = None,
) -> None:
    """Set global CONFIG variable"""
    global CONFIG
    CONFIG = load_config(config_path=config_path)


def print_config(prefix: str = None) -> None:
    """Print configuration"""

    if prefix:
        mrich.var(
            f"CONFIG ({prefix})",
            {k: v for k, v in CONFIG.items() if k.startswith(prefix)},
        )
    else:
        mrich.var("CONFIG", CONFIG)


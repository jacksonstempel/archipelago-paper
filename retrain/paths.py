"""The data root: one folder that holds every input, intermediate file and result (layout in data/README.md).

It is given by a command-line flag (--data) or, without one, by the environment variable ARCHIPELAGO_DATA. There is no
default: a missing root is an error, so nothing is ever read from or written to an unintended place.
"""
import os
from pathlib import Path

ENVIRONMENT = "ARCHIPELAGO_DATA"


def data_root(path=None):
    """The data root: `path` if given, else $ARCHIPELAGO_DATA. Fails if neither is set or the folder does not exist."""
    path = path or os.environ.get(ENVIRONMENT)
    if not path:
        raise SystemExit(f"No data root: pass --data DIR or set {ENVIRONMENT} (layout in data/README.md)")
    root = Path(path).expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"Data root {root} does not exist")
    return root


def add_argument(parser):
    """Add the --data flag to an argparse parser."""
    parser.add_argument("--data", metavar="DIR", help=f"data root (default: ${ENVIRONMENT})")

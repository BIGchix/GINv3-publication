"""Repository-default paths for GIN network and curation files."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_GIN_VERSION = "GINv3.11a"
DEFAULT_NODE_LIST_DIR = REPO_ROOT / "nodeList"
DEFAULT_SYMBOL_BLACKLIST = REPO_ROOT / "config" / "gin_symbol_blacklist.txt"
DEFAULT_SYMBOL_ALIASES = REPO_ROOT / "config" / "symbol_aliases.tsv"


def _resolve_default_sifi() -> Path:
    for candidate in (
        REPO_ROOT / "GINv3.11a.sifi",
        REPO_ROOT / "data" / "GINv3.11a.sifi",
    ):
        if candidate.is_file():
            return candidate
    return REPO_ROOT / "GINv3.11a.sifi"


DEFAULT_SIFI = _resolve_default_sifi()

"""GIN bipartite multi-omics integration utilities."""

from gin_tool.user_stories import USER_STORIES, UserStory

__all__ = [
    "USER_STORIES",
    "UserStory",
    "active_tf",
    "analyze",
    "benchmark",
    "context_subgraph",
    "enrichment",
    "flows",
    "io",
    "mapping",
    "module_flow",
    "modules",
    "propagate",
    "regulatory_filter",
    "schematic",
    "signaling_subgraph",
    "spn_export",
    "spn_rehydrate",
    "subgraph",
]


def __getattr__(name: str):
    if name in __all__ and name not in {"USER_STORIES", "UserStory"}:
        import importlib

        return importlib.import_module(f"gin_tool.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

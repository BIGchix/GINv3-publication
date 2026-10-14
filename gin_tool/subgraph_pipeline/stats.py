"""Node / edge count helpers for subgraph pipeline QC."""

from __future__ import annotations

from typing import Any

from gin_tool.subgraph_pipeline.manifest import NodeCountStats


def node_count_stats(walk: Any, *, n_seeds_no_loop: int = 0) -> NodeCountStats:
    return NodeCountStats(
        n_real_nodes=len(getattr(walk, "real_nodes", set()) or set()),
        n_virtual_nodes=len(getattr(walk, "virtual_nodes", set()) or set()),
        n_gin_edges=len(getattr(walk, "gin_edges", set()) or set()),
        n_compressed_edges=len(getattr(walk, "compressed_edges", set()) or set()),
        n_seeds_no_loop=n_seeds_no_loop,
    )

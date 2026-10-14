"""Adaptive caps aligned with scBionet subgraph presets."""

from __future__ import annotations

from typing import Any

import numpy as np


def resolve_max_real_nodes(n_union_genes: int) -> tuple[int, dict[str, Any]]:
    """``clamp(round(3.5 * n_union), 600, 1100)`` — same as scBionet."""
    raw = int(round(3.5 * max(int(n_union_genes), 0)))
    resolved = int(np.clip(raw, 600, 1100))
    return resolved, {
        "n_union_genes": int(n_union_genes),
        "raw": raw,
        "max_real_nodes_resolved": resolved,
        "formula": "clamp(round(3.5*n_union), 600, 1100)",
    }


def fast_preset_config_overrides(n_union_genes: int) -> dict[str, Any]:
    """scBionet ``low_noncan`` walk settings (bipartite pipeline keeps virtual nodes)."""
    max_nodes, budget_meta = resolve_max_real_nodes(n_union_genes)
    return {
        "max_real_hops": 4,
        "max_real_nodes": max_nodes,
        "require_seed_self_loops": False,
        "enable_connectivity_l2": False,
        "regulatory_relaxed": False,
        "budget_meta": budget_meta,
    }

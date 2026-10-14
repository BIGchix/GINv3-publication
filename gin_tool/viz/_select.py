"""Select a tractable vertex subset for network layout."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from gin_tool.mapping import is_virtual_gin_node


def select_viz_node_indices_named(
    sub_sym: sparse.csr_matrix,
    scores: np.ndarray,
    seed_indices: list[int],
    names: list[str],
    *,
    max_nodes: int,
) -> set[int]:
    """Like ``select_viz_node_indices`` but drop low-score virtual nodes first when trimming."""
    n = int(sub_sym.shape[0])
    if n <= max_nodes:
        return set(range(n))

    selected: set[int] = set(int(i) for i in seed_indices if 0 <= i < n)
    for s in list(selected):
        lo, hi = sub_sym.indptr[s], sub_sym.indptr[s + 1]
        for j in sub_sym.indices[lo:hi]:
            selected.add(int(j))

    order = np.argsort(-scores)
    for i in order:
        if len(selected) >= max_nodes:
            break
        selected.add(int(i))

    if len(selected) <= max_nodes:
        return selected

    ranked = sorted(
        selected,
        key=lambda i: (
            i in seed_indices,
            float(scores[i]),
            0 if is_virtual_gin_node(names[i]) else 1,
        ),
        reverse=True,
    )
    return set(ranked[:max_nodes])

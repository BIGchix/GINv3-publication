"""Relate Louvain modules (clusters) to directed information-flow paths."""

from __future__ import annotations

import json
from collections import Counter
from typing import Callable, Sequence

import numpy as np
from scipy import sparse

from gin_tool.mapping import is_virtual_gin_node


def path_module_incidence_matrix(
    paths: Sequence[Sequence[int]],
    modules: Sequence[frozenset[int]],
) -> np.ndarray:
    """
    Shape ``(num_paths, num_modules)``: count of vertices of path ``p`` lying in module ``m``.
    """
    p_n = len(paths)
    m_n = len(modules)
    mat = np.zeros((p_n, m_n), dtype=np.int32)
    for pi, path in enumerate(paths):
        s = set(path)
        for mi, mod in enumerate(modules):
            mat[pi, mi] = len(s & mod)
    return mat


def node_primary_module(node: int, modules: Sequence[frozenset[int]]) -> int | None:
    for mi, mod in enumerate(modules):
        if node in mod:
            return mi
    return None


def interface_real_nodes(
    paths: Sequence[Sequence[int]],
    modules: Sequence[frozenset[int]],
    adj_sym: sparse.csr_matrix,
    name_fn: Callable[[int], str],
) -> set[int]:
    """
    Real (non-virtual) nodes that either:

    - appear on **≥2** top paths (flow–flow interface), or
    - lie on a **sym** edge whose endpoints belong to **different** modules
      (module–module interface, often mechanistic bridges).
    """
    cross_path: set[int] = set()
    cnt: Counter[int] = Counter()
    for p in paths:
        for v in p:
            cnt[v] += 1
    for v, k in cnt.items():
        if k >= 2:
            cross_path.add(v)

    node_mod: dict[int, int] = {}
    for mi, mod in enumerate(modules):
        for v in mod:
            node_mod[v] = mi

    cross_mod: set[int] = set()
    n = adj_sym.shape[0]
    for u in range(n):
        lo, hi = adj_sym.indptr[u], adj_sym.indptr[u + 1]
        for idx in range(lo, hi):
            v = int(adj_sym.indices[idx])
            if v <= u:
                continue
            mu, mv = node_mod.get(u), node_mod.get(v)
            if mu is not None and mv is not None and mu != mv:
                cross_mod.add(u)
                cross_mod.add(v)

    raw = cross_path | cross_mod
    return {v for v in raw if not is_virtual_gin_node(name_fn(v))}


def write_module_flow_json(
    path: str,
    incidence: np.ndarray,
    interfaces: set[int],
    names: Sequence[str],
) -> None:
    payload = {
        "path_module_incidence": incidence.tolist(),
        "interface_real_nodes": sorted(names[i] for i in interfaces),
        "notes": (
            "Rows: ranked information-flow paths; columns: Louvain modules. "
            "High counts mean a path spends many steps inside that module. "
            "Interface nodes are candidates for schematic detail labels."
        ),
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

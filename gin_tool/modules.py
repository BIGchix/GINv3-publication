"""Functional modules: community structure and inter-module graph distances."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Sequence

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path


def louvain_communities_subgraph(adj_sym: sparse.csr_matrix) -> list[frozenset[int]]:
    """
    Partition vertices of an **undirected** subgraph using Louvain modularity.

    ``adj_sym`` must be symmetric (e.g. k-hop support from ``extract_khop_region``).
    """
    G = nx.from_scipy_sparse_array(adj_sym, create_using=nx.Graph)
    comms = nx.community.louvain_communities(G, seed=42)
    return [frozenset(c) for c in comms]


def module_pairwise_min_geodesic(
    adj_sym: sparse.csr_matrix,
    modules: Sequence[frozenset[int]],
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each unordered pair of modules, compute the **minimum shortest-path length**
    between any node in module A and any node in module B (graph geodesic on the
    undirected support).

    Returns ``(dist_matrix, unreachable_mask)`` of shape (M, M). Diagonal is 0;
    ``unreachable_mask`` is True where no finite path exists between the two sets.
    """
    n = adj_sym.shape[0]
    m = len(modules)
    if m == 0:
        return np.zeros((0, 0)), np.zeros((0, 0), dtype=bool)
    if n > 6000:
        raise ValueError(
            f"Subgraph has n={n} nodes; all-pairs shortest paths are too heavy. "
            "Increase locality (smaller k / fewer seeds) or implement sampling."
        )
    sp = shortest_path(adj_sym, directed=False, unweighted=True)
    dist = np.full((m, m), np.inf, dtype=np.float64)
    unreachable = np.zeros((m, m), dtype=bool)
    for i in range(m):
        dist[i, i] = 0.0
        for j in range(i + 1, m):
            dmin = np.inf
            for u in modules[i]:
                row = sp[u]
                for v in modules[j]:
                    d = float(row[v])
                    if d < dmin:
                        dmin = d
            dist[i, j] = dist[j, i] = dmin
            unreachable[i, j] = unreachable[j, i] = not np.isfinite(dmin)
    return dist, unreachable


def write_module_distance_tsv(
    path: str,
    modules: Sequence[frozenset[int]],
    dist: np.ndarray,
    unreachable: np.ndarray,
) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("module_i\tmodule_j\tmin_geodesic\tunreachable\n")
        m = len(modules)
        for i in range(m):
            for j in range(i + 1, m):
                d = dist[i, j]
                unr = unreachable[i, j]
                f.write(f"{i}\t{j}\t{d if np.isfinite(d) else 'NA'}\t{int(unr)}\n")


def modules_to_json(path: str, modules: Sequence[frozenset[int]], names: Sequence[str]) -> None:
    out = []
    for mi, nodes in enumerate(modules):
        out.append(
            {
                "id": mi,
                "size": len(nodes),
                "nodes": sorted(names[i] for i in nodes),
            }
        )
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

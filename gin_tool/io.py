"""Load GIN .sifi edge lists and build sparse adjacency (subset or full)."""

from __future__ import annotations

import gzip
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import numpy as np
from scipy import sparse


REL_POS = 1
REL_NEG = -1


def _is_gzip_path(path: str | Path) -> bool:
    name = Path(path).name.lower()
    return name.endswith(".gz")


@contextmanager
def open_sifi_text(
    path: str | Path,
    encoding: str = "utf-8",
    errors: str = "strict",
) -> Iterator[TextIO]:
    """Open a ``.sifi`` or ``.sifi.gz`` edge list as text (transparent gzip)."""
    p = Path(path)
    if _is_gzip_path(p):
        with gzip.open(p, "rt", encoding=encoding, errors=errors) as fh:
            yield fh  # type: ignore[misc]
    else:
        with p.open("r", encoding=encoding, errors=errors) as fh:
            yield fh


def collect_unique_node_names(
    path: str,
    encoding: str = "utf-8",
    max_lines: int | None = None,
) -> set[str]:
    """
    Single pass over the .sifi file to collect node names (no edge buffer).

    Use this for ``gin_node_is_member`` mapping QC before loading all edges.
    Accepts plain ``.sifi`` or ``.sifi.gz``.
    """
    names: set[str] = set()
    n = 0
    with open_sifi_text(path, encoding=encoding) as f:
        for line in f:
            if max_lines is not None and n >= max_lines:
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            names.add(parts[0])
            names.add(parts[1])
            n += 1
    return names


def iter_sifi_edges(
    path: str,
    encoding: str = "utf-8",
    max_lines: int | None = None,
) -> Iterator[tuple[str, str, int]]:
    """
    Stream (source, target, sign) from a GIN tab-separated edge list.

    Third column must be ``positive`` or ``negative`` (case-sensitive as in GIN).
    Accepts plain ``.sifi`` or ``.sifi.gz``.
    """
    n = 0
    with open_sifi_text(path, encoding=encoding) as f:
        for line in f:
            if max_lines is not None and n >= max_lines:
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            s, t, rel = parts[0], parts[1], parts[2]
            if rel == "positive":
                w = REL_POS
            elif rel == "negative":
                w = REL_NEG
            else:
                continue
            yield s, t, w
            n += 1


@dataclass
class IndexedEdges:
    """Integer-indexed directed edges with optional relation weights."""

    num_nodes: int
    sources: np.ndarray  # int32
    targets: np.ndarray  # int32
    weights: np.ndarray  # int8 or float
    index_to_name: list[str]
    name_to_index: dict[str, int]


def build_indexed_edges(
    path: str,
    *,
    max_lines: int | None = None,
    encoding: str = "utf-8",
) -> IndexedEdges:
    """
    Two-pass: collect node vocabulary, then fill edge arrays.

    For the full ~16M-edge GIN file this is memory- and time-heavy; use
    ``max_lines`` for development or subgraph pipelines that pre-filter the TSV.
    """
    # First pass: unique nodes
    names: list[str] = []
    seen: set[str] = set()
    buf_s: list[str] = []
    buf_t: list[str] = []
    buf_w: list[int] = []
    for s, t, w in iter_sifi_edges(path, encoding=encoding, max_lines=max_lines):
        buf_s.append(s)
        buf_t.append(t)
        buf_w.append(w)
        for n in (s, t):
            if n not in seen:
                seen.add(n)
                names.append(n)
    name_to_index = {n: i for i, n in enumerate(names)}
    m = len(buf_s)
    sources = np.empty(m, dtype=np.int32)
    targets = np.empty(m, dtype=np.int32)
    weights = np.empty(m, dtype=np.int8)
    for i in range(m):
        sources[i] = name_to_index[buf_s[i]]
        targets[i] = name_to_index[buf_t[i]]
        weights[i] = buf_w[i]
    return IndexedEdges(
        num_nodes=len(names),
        sources=sources,
        targets=targets,
        weights=weights,
        index_to_name=names,
        name_to_index=name_to_index,
    )


def directed_csr_from_indexed(ed: IndexedEdges) -> sparse.csr_matrix:
    """Directed weighted adjacency (rows=source, cols=target)."""
    return sparse.coo_matrix(
        (ed.weights.astype(np.float64), (ed.sources, ed.targets)),
        shape=(ed.num_nodes, ed.num_nodes),
    ).tocsr()


def symmetrize_undirected_unweighted(adj: sparse.csr_matrix) -> sparse.csr_matrix:
    """k-hop neighborhood expansion on an undirected unweighted support graph."""
    sup = (adj != 0).astype(np.float64)
    sym = sup + sup.T
    sym.data[:] = 1.0
    sym.eliminate_zeros()
    return sym.tocsr()


def khop_nodes_per_seed(
    adj: sparse.csr_matrix,
    seed_indices: Sequence[int],
    k_by_seed: Sequence[int],
) -> set[int]:
    """Union of k-hop balls where each seed may use a different ``k``."""
    if len(seed_indices) != len(k_by_seed):
        raise ValueError("seed_indices and k_by_seed must have the same length")
    nodes: set[int] = set()
    for idx, k in zip(seed_indices, k_by_seed):
        nodes |= khop_nodes_csr(adj, [int(idx)], int(k))
    return nodes


def khop_nodes_csr(adj: sparse.csr_matrix, seeds: Sequence[int], k: int) -> set[int]:
    """Return all nodes within k hops (unweighted) on symmetric CSR adjacency."""
    frontier: set[int] = set(int(x) for x in seeds)
    visited: set[int] = set(frontier)
    for _ in range(k):
        next_front: set[int] = set()
        for u in frontier:
            lo, hi = adj.indptr[u], adj.indptr[u + 1]
            neigh = adj.indices[lo:hi]
            for v in neigh:
                if v not in visited:
                    visited.add(v)
                    next_front.add(v)
        frontier = next_front
        if not frontier:
            break
    return visited


def khop_frontier_sizes(
    adj: sparse.csr_matrix, seeds: Sequence[int], k: int
) -> list[int]:
    """
    Cumulative visited count after each hop expansion (includes seeds at step 0).

    ``sizes[0]`` = number of seeds; ``sizes[i]`` after i hop rounds (i in 1..k).
    """
    frontier: set[int] = set(int(x) for x in seeds)
    visited: set[int] = set(frontier)
    sizes: list[int] = [len(visited)]
    for _ in range(k):
        next_front: set[int] = set()
        for u in frontier:
            lo, hi = adj.indptr[u], adj.indptr[u + 1]
            for v in adj.indices[lo:hi]:
                v = int(v)
                if v not in visited:
                    visited.add(v)
                    next_front.add(v)
        frontier = next_front
        sizes.append(len(visited))
        if not frontier:
            while len(sizes) < k + 1:
                sizes.append(len(visited))
            break
    return sizes


def induced_subgraph_csr(
    adj: sparse.csr_matrix,
    node_subset: set[int],
) -> tuple[sparse.csr_matrix, dict[int, int], dict[int, int]]:
    """
    Extract CSR on relabeled 0..n-1 indices for nodes in ``node_subset``.

    Returns (sub_adj, old_to_new, new_to_old).
    """
    old_ids = sorted(node_subset)
    old_to_new = {o: j for j, o in enumerate(old_ids)}
    n = len(old_ids)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for new_u, old_u in enumerate(old_ids):
        lo, hi = adj.indptr[old_u], adj.indptr[old_u + 1]
        for old_v, val in zip(adj.indices[lo:hi], adj.data[lo:hi]):
            if old_v in old_to_new:
                rows.append(new_u)
                cols.append(old_to_new[old_v])
                data.append(float(val))
    sub = sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    new_to_old = {j: o for o, j in old_to_new.items()}
    return sub, old_to_new, new_to_old

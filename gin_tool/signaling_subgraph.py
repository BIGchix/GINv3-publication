"""GINtoSPN-style signaling walk on strict bipartite GIN with real-only output."""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse

from gin_tool.regulatory_filter import allow_directed_edge
from gin_tool.io import IndexedEdges, directed_csr_from_indexed, symmetrize_undirected_unweighted
from gin_tool.mapping import is_virtual_gin_node
from gin_tool.regulatory_filter import RegulatoryPrefilterStats, filter_indexed_edges_regulatory


def _out_neighbors_csr(adj: sparse.csr_matrix, u: int) -> np.ndarray:
    lo, hi = adj.indptr[u], adj.indptr[u + 1]
    return adj.indices[lo:hi]



@dataclass
class SignalingWalkResult:
    path_pairs: list[tuple[str, str]] = field(default_factory=list)
    real_nodes: set[str] = field(default_factory=set)
    compressed_edges: set[tuple[str, str]] = field(default_factory=set)
    pair_min_real_hops: dict[tuple[str, str], int] = field(default_factory=dict)
    virtual_nodes_visited: int = 0
    status: str = "ok"


def _edge_allowed(
    u_name: str,
    v_name: str,
    *,
    filter_regulatory: bool,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
) -> bool:
    if not filter_regulatory:
        return True
    return allow_directed_edge(
        u_name,
        v_name,
        active_tf_names=active_tf_names,
        seed_symbols=seed_symbols,
        seed_gin_names=seed_gin_names,
    )


def _bfs_from_seed(
    s_idx: int,
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    input_idx: set[int],
    *,
    max_real_hops: int,
    max_real_nodes: int,
    filter_regulatory: bool,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
) -> tuple[dict[int, int], dict[int, set[int]], int]:
    """BFS in real-hop layers; returns dist, predecessors on real nodes, virtual visit count."""
    real_dist: dict[int, int] = {s_idx: 0}
    real_parents: dict[int, set[int]] = {s_idx: set()}
    virtual_visits = 0
    frontier_real = [s_idx]

    for hop in range(max_real_hops):
        if not frontier_real:
            break

        v_from: dict[int, set[int]] = defaultdict(set)
        for r in frontier_real:
            u_name = index_to_name[r]
            for v in _out_neighbors_csr(adj, r):
                v = int(v)
                v_name = index_to_name[v]
                if not is_virtual_gin_node(v_name):
                    continue
                if not _edge_allowed(
                    u_name,
                    v_name,
                    filter_regulatory=filter_regulatory,
                    active_tf_names=active_tf_names,
                    seed_symbols=seed_symbols,
                    seed_gin_names=seed_gin_names,
                ):
                    continue
                v_from[v].add(r)

        if not v_from:
            break

        next_preds: dict[int, set[int]] = defaultdict(set)
        new_hop = hop + 1
        for v, pred_reals in v_from.items():
            virtual_visits += 1
            v_name = index_to_name[v]
            for r2 in _out_neighbors_csr(adj, v):
                r2 = int(r2)
                r2_name = index_to_name[r2]
                if is_virtual_gin_node(r2_name):
                    continue
                if not _edge_allowed(
                    v_name,
                    r2_name,
                    filter_regulatory=filter_regulatory,
                    active_tf_names=active_tf_names,
                    seed_symbols=seed_symbols,
                    seed_gin_names=seed_gin_names,
                ):
                    continue
                old = real_dist.get(r2)
                if old is not None and old < new_hop:
                    continue
                next_preds[r2].update(pred_reals)

        if not next_preds:
            break

        frontier_real = []
        for r2, preds in next_preds.items():
            old = real_dist.get(r2)
            if old is not None and old < new_hop:
                continue
            real_dist[r2] = new_hop
            if old == new_hop:
                real_parents[r2].update(preds)
            else:
                real_parents[r2] = set(preds)
            frontier_real.append(r2)
            if len(real_dist) >= max_real_nodes:
                return real_dist, real_parents, virtual_visits

    return real_dist, real_parents, virtual_visits


def _backtrack_real_nodes(
    s_idx: int,
    e_idx: int,
    hop: int,
    real_parents: dict[int, set[int]],
) -> tuple[set[int], set[tuple[int, int]]]:
    """Collect real indices and compressed edges on all shortest real-hop paths."""
    nodes: set[int] = set()
    edges: set[tuple[int, int]] = set()

    def walk(cur: int, remaining: int) -> None:
        nodes.add(cur)
        if remaining == 0:
            return
        for pred in real_parents.get(cur, ()):
            edges.add((pred, cur))
            walk(pred, remaining - 1)

    walk(e_idx, hop)
    nodes.add(s_idx)
    return nodes, edges


def _prune_real_nodes_to_seeds(
    real_indices: set[int],
    compressed_edges: set[tuple[int, int]],
    seed_indices: set[int],
) -> set[int]:
    """Keep real nodes with a reverse path to any seed in the compressed real graph."""
    if not real_indices:
        return set()
    preds: dict[int, set[int]] = defaultdict(set)
    for u, v in compressed_edges:
        if u in real_indices and v in real_indices:
            preds[v].add(u)

    keep: set[int] = set()
    queue: deque[int] = deque()
    for s in seed_indices:
        if s in real_indices:
            keep.add(s)
            queue.append(s)

    while queue:
        n = queue.popleft()
        for p in preds.get(n, ()):
            if p not in keep:
                keep.add(p)
                queue.append(p)
    return keep


def signaling_real_subgraph(
    input_gs: set[str],
    adj_dir: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    *,
    max_real_hops: int = 12,
    max_real_nodes: int = 100_000,
    filter_regulatory: bool = False,
    active_tf_names: set[str] | None = None,
    seed_symbols: set[str] | None = None,
    seed_gin_names: set[str] | None = None,
) -> SignalingWalkResult:
    """
    Signaling walk with rule③ (continue through non-input reals).

    Traverses R-V-R in the bipartite graph but outputs only real nodes and
    compressed real→real edges (each R-V-R segment is one hop).
    """
    active_tf_names = active_tf_names or set()
    seed_symbols = seed_symbols or set()
    seed_gin_names = seed_gin_names or set(input_gs)

    seed_indices = [name_to_index[n] for n in input_gs if n in name_to_index]
    if not seed_indices:
        return SignalingWalkResult(status="no_seeds")

    input_idx = set(seed_indices)
    global_pair_hops: dict[tuple[int, int], int] = {}
    per_seed: dict[int, tuple[dict[int, int], dict[int, set[int]]]] = {}
    virtual_visits = 0

    for s_idx in seed_indices:
        dist, parents, vv = _bfs_from_seed(
            s_idx,
            adj_dir,
            index_to_name,
            input_idx,
            max_real_hops=max_real_hops,
            max_real_nodes=max_real_nodes,
            filter_regulatory=filter_regulatory,
            active_tf_names=active_tf_names,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_names,
        )
        virtual_visits += vv
        per_seed[s_idx] = (dist, parents)
        for e_idx, hop in dist.items():
            if e_idx == s_idx or e_idx not in input_idx:
                continue
            key = (s_idx, e_idx)
            prev = global_pair_hops.get(key)
            if prev is None or hop < prev:
                global_pair_hops[key] = hop

    real_idx_set: set[int] = set()
    compressed: set[tuple[int, int]] = set()
    for (s_idx, e_idx), hop in global_pair_hops.items():
        dist, parents = per_seed[s_idx]
        if dist.get(e_idx) != hop:
            continue
        nodes, edges = _backtrack_real_nodes(s_idx, e_idx, hop, parents)
        real_idx_set.update(nodes)
        compressed.update(edges)

    seed_idx_set = set(seed_indices)
    pruned_idx = _prune_real_nodes_to_seeds(real_idx_set, compressed, seed_idx_set)

    real_nodes = {index_to_name[i] for i in pruned_idx}
    for name in list(real_nodes):
        if name.endswith(".state"):
            real_nodes.add(name.replace(".state", "") + ".RNA")

    path_pairs: set[tuple[str, str]] = set()
    pair_min: dict[tuple[str, str], int] = {}
    for (s_idx, e_idx), hop in global_pair_hops.items():
        a, b = index_to_name[s_idx], index_to_name[e_idx]
        if a == b:
            continue
        pair = (a, b)
        path_pairs.add(pair)
        pair_min[pair] = hop

    compressed_named = {
        (index_to_name[u], index_to_name[v])
        for u, v in compressed
        if u in pruned_idx and v in pruned_idx
    }

    return SignalingWalkResult(
        path_pairs=sorted(path_pairs),
        real_nodes=real_nodes,
        compressed_edges=compressed_named,
        pair_min_real_hops=pair_min,
        virtual_nodes_visited=virtual_visits,
        status="ok",
    )


def signaling_real_subgraph_from_edges(
    ed: IndexedEdges,
    input_gs: set[str],
    *,
    max_real_hops: int = 12,
    max_real_nodes: int = 100_000,
    prefilter_regulatory: bool = True,
    filter_regulatory: bool = False,
    active_tf_names: set[str] | None = None,
    seed_symbols: set[str] | None = None,
    seed_gin_names: set[str] | None = None,
    tf_pool: set[str] | None = None,
) -> tuple[SignalingWalkResult, sparse.csr_matrix, RegulatoryPrefilterStats | None]:
    """
    Signaling walk on ``ed``, optionally with regulatory pre-filter (default on).

    Pre-filter drops protein/TF→DNA edges except from active TFs (seeds-only
    input; active TFs are **not** added to the walk frontier). When
    ``prefilter_regulatory=True``, ``filter_regulatory`` defaults to off because
    the graph is already masked.
    """
    active_tf_names = active_tf_names or set()
    seed_symbols = seed_symbols or set()
    seed_gin_names = seed_gin_names or set(input_gs)
    tf_pool = tf_pool or set()

    prefilter_stats: RegulatoryPrefilterStats | None = None
    if prefilter_regulatory:
        work_ed, prefilter_stats = filter_indexed_edges_regulatory(
            ed,
            active_tf_names=active_tf_names,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_names,
            tf_pool=tf_pool,
        )
        use_runtime_filter = filter_regulatory
    else:
        work_ed = ed
        use_runtime_filter = filter_regulatory

    adj_dir = directed_csr_from_indexed(work_ed)
    result = signaling_real_subgraph(
        input_gs,
        adj_dir,
        work_ed.index_to_name,
        work_ed.name_to_index,
        max_real_hops=max_real_hops,
        max_real_nodes=max_real_nodes,
        filter_regulatory=use_runtime_filter,
        active_tf_names=active_tf_names,
        seed_symbols=seed_symbols,
        seed_gin_names=seed_gin_names,
    )
    return result, adj_dir, prefilter_stats


def build_sym_subgraph_from_walk(
    result: SignalingWalkResult,
) -> tuple[sparse.csr_matrix, dict[str, int], list[str]]:
    """
    Undirected unweighted subgraph from compressed real→real walk edges.

    Returns ``(sub_sym, name_to_local, ordered_names)``.
    """
    nodes = sorted(result.real_nodes)
    name_to_local = {n: i for i, n in enumerate(nodes)}
    rows: list[int] = []
    cols: list[int] = []
    for a, b in result.compressed_edges:
        if a not in name_to_local or b not in name_to_local:
            continue
        i, j = name_to_local[a], name_to_local[b]
        rows.extend([i, j])
        cols.extend([j, i])
    n = len(nodes)
    if not rows:
        return sparse.csr_matrix((n, n)), name_to_local, nodes
    data = np.ones(len(rows), dtype=np.float64)
    sym = sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    sym.eliminate_zeros()
    return sym, name_to_local, nodes

"""GIN signaling walk preserving virtual transition nodes (no R-R compression)."""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse

from gin_tool.mapping import is_virtual_gin_node
from gin_tool.regulatory_filter import allow_directed_edge
from gin_tool.signaling_subgraph import (
    _backtrack_real_nodes,
    _bfs_from_seed,
    _edge_allowed,
    _out_neighbors_csr,
    _prune_real_nodes_to_seeds,
)
from gin_tool.spn_rehydrate import find_virtual_intermediates


@dataclass
class SignalingBipartiteWalkResult:
    path_pairs: list[tuple[str, str]] = field(default_factory=list)
    real_nodes: set[str] = field(default_factory=set)
    virtual_nodes: set[str] = field(default_factory=set)
    gin_edges: set[tuple[str, str, float]] = field(default_factory=set)
    compressed_edges: set[tuple[str, str]] = field(default_factory=set)
    pair_min_real_hops: dict[tuple[str, str], int] = field(default_factory=dict)
    virtual_nodes_visited: int = 0
    hit_max_real_nodes: bool = False
    status: str = "ok"


def _edge_weight(adj: sparse.csr_matrix, u: int, v: int) -> float:
    return float(adj[u, v])


def add_bipartite_segment(
    walk: SignalingBipartiteWalkResult,
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    source_idx: int,
    target_idx: int,
) -> None:
    """Expand one real-hop (R-V-R) into gin_edges + virtual_nodes."""
    intermediates = find_virtual_intermediates(adj, index_to_name, source_idx, target_idx)
    src = index_to_name[source_idx]
    tgt = index_to_name[target_idx]
    walk.compressed_edges.add((src, tgt))

    if intermediates:
        via, _sign = intermediates[0]
        walk.virtual_nodes.add(via)
        vi = name_to_index[via]
        walk.gin_edges.add((src, via, _edge_weight(adj, source_idx, vi)))
        walk.gin_edges.add((via, tgt, _edge_weight(adj, vi, target_idx)))
        for alt, _ in intermediates[1:]:
            walk.virtual_nodes.add(alt)
    elif abs(_edge_weight(adj, source_idx, target_idx)) > 0:
        walk.gin_edges.add((src, tgt, _edge_weight(adj, source_idx, target_idx)))


def add_real_hop_path(
    walk: SignalingBipartiteWalkResult,
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    s_idx: int,
    e_idx: int,
    hop: int,
    real_parents: dict[int, set[int]],
) -> None:
    nodes, edges = _backtrack_real_nodes(s_idx, e_idx, hop, real_parents)
    walk.real_nodes.update(index_to_name[i] for i in nodes)
    for u, v in edges:
        add_bipartite_segment(walk, adj, index_to_name, name_to_index, u, v)


def rebuild_gin_edges_from_compressed(
    walk: SignalingBipartiteWalkResult,
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
) -> None:
    walk.gin_edges.clear()
    walk.virtual_nodes.clear()
    for src, tgt in sorted(walk.compressed_edges):
        if src not in name_to_index or tgt not in name_to_index:
            continue
        add_bipartite_segment(
            walk,
            adj,
            index_to_name,
            name_to_index,
            name_to_index[src],
            name_to_index[tgt],
        )


def signaling_bipartite_subgraph(
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
) -> SignalingBipartiteWalkResult:
    """
    Signaling walk (rule③) outputting real + virtual nodes and full gin_edges.

    Real-hop BFS matches ``signaling_real_subgraph``; each R-V-R segment is
    rehydrated instead of compressed to a single logical edge.
    """
    active_tf_names = active_tf_names or set()
    seed_symbols = seed_symbols or set()
    seed_gin_names = seed_gin_names or set(input_gs)

    seed_indices = [name_to_index[n] for n in input_gs if n in name_to_index]
    if not seed_indices:
        return SignalingBipartiteWalkResult(status="no_seeds")

    input_idx = set(seed_indices)
    global_pair_hops: dict[tuple[int, int], int] = {}
    per_seed: dict[int, tuple[dict[int, int], dict[int, set[int]]]] = {}
    virtual_visits = 0

    walk = SignalingBipartiteWalkResult()
    hit_cap = False

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
        if len(dist) >= max_real_nodes:
            hit_cap = True
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

    for (s_idx, e_idx), hop in global_pair_hops.items():
        if s_idx not in pruned_idx or e_idx not in pruned_idx:
            continue
        dist, parents = per_seed[s_idx]
        if dist.get(e_idx) != hop:
            continue
        add_real_hop_path(walk, adj_dir, index_to_name, name_to_index, s_idx, e_idx, hop, parents)

    walk.real_nodes = {index_to_name[i] for i in pruned_idx}
    walk.compressed_edges = {
        (index_to_name[u], index_to_name[v])
        for u, v in compressed
        if u in pruned_idx and v in pruned_idx
    }
    rebuild_gin_edges_from_compressed(walk, adj_dir, index_to_name, name_to_index)

    path_pairs: set[tuple[str, str]] = set()
    pair_min: dict[tuple[str, str], int] = {}
    for (s_idx, e_idx), hop in global_pair_hops.items():
        a, b = index_to_name[s_idx], index_to_name[e_idx]
        if a == b or a not in walk.real_nodes or b not in walk.real_nodes:
            continue
        path_pairs.add((a, b))
        pair_min[(a, b)] = hop

    walk.path_pairs = sorted(path_pairs)
    walk.pair_min_real_hops = pair_min
    walk.virtual_nodes_visited = virtual_visits
    walk.hit_max_real_nodes = hit_cap
    walk.status = "ok"
    return walk


def add_bipartite_path_edges(
    walk: SignalingBipartiteWalkResult,
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    start: str,
    goal: str,
    max_hops: int,
) -> bool:
    """BFS on full GIN; merge all nodes/edges on shortest path."""
    if start not in name_to_index or goal not in name_to_index:
        return False
    s, t = name_to_index[start], name_to_index[goal]
    dist = {s: 0}
    parent: dict[int, int] = {}
    q = deque([s])
    while q:
        u = q.popleft()
        if u == t:
            break
        if dist[u] >= max_hops:
            continue
        lo, hi = adj.indptr[u], adj.indptr[u + 1]
        for v in adj.indices[lo:hi]:
            v = int(v)
            if v not in dist:
                dist[v] = dist[u] + 1
                parent[v] = u
                q.append(v)
    if t not in dist:
        return False

    chain: list[int] = []
    cur = t
    while cur in parent:
        chain.append(cur)
        cur = parent[cur]
    chain.append(cur)
    chain.reverse()

    for idx in chain:
        name = index_to_name[idx]
        if is_virtual_gin_node(name):
            walk.virtual_nodes.add(name)
        else:
            walk.real_nodes.add(name)

    for i in range(len(chain) - 1):
        u, v = chain[i], chain[i + 1]
        a, b = index_to_name[u], index_to_name[v]
        w = _edge_weight(adj, u, v)
        walk.gin_edges.add((a, b, w))
        if not is_virtual_gin_node(a) and not is_virtual_gin_node(b):
            walk.compressed_edges.add((a, b))

    return True

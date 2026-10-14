"""GIN bipartite subgraph connectivity repair (L2 bridge + L1 largest component)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from gin_tool.subgraph_pipeline.manifest import ConnectivityQC, NodeCountStats
from gin_tool.subgraph_pipeline.signaling_bipartite import (
    SignalingBipartiteWalkResult,
    add_bipartite_path_edges,
    rebuild_gin_edges_from_compressed,
)
from gin_tool.subgraph_pipeline.stats import node_count_stats


@dataclass
class ConnectivityResult:
    walk: SignalingBipartiteWalkResult
    n_components_before: int
    n_components_after: int
    largest_component_frac: float
    nodes_removed: int
    bridge_edges_added: int
    seeds_in_real_nodes: int
    seeds_total: int
    seed_gin_nodes: list[str]


def _sym_graph_from_walk(walk: SignalingBipartiteWalkResult) -> tuple[sparse.csr_matrix, dict[str, int], list[str]]:
    nodes = sorted(walk.real_nodes)
    name_to_idx = {n: i for i, n in enumerate(nodes)}
    rows, cols = [], []
    for a, b in walk.compressed_edges:
        if a not in name_to_idx or b not in name_to_idx:
            continue
        i, j = name_to_idx[a], name_to_idx[b]
        rows.extend([i, j])
        cols.extend([j, i])
    n = len(nodes)
    if not rows:
        return sparse.csr_matrix((n, n)), name_to_idx, nodes
    data = np.ones(len(rows), dtype=np.float32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n)), name_to_idx, nodes


def _component_labels(walk: SignalingBipartiteWalkResult) -> tuple[int, np.ndarray, list[str]]:
    g, _, nodes = _sym_graph_from_walk(walk)
    if len(nodes) <= 1:
        return 1, np.zeros(len(nodes), dtype=int), nodes
    n_comp, labels = connected_components(g, directed=False)
    return n_comp, labels, nodes


def largest_connected_subgraph(
    walk: SignalingBipartiteWalkResult,
) -> tuple[SignalingBipartiteWalkResult, int, float, int]:
    g, _, nodes = _sym_graph_from_walk(walk)
    n_comp_before, _ = connected_components(g, directed=False) if len(nodes) > 1 else (1, np.zeros(1))
    if len(nodes) <= 1:
        return walk, n_comp_before, 1.0, 0

    n_comp, labels = connected_components(g, directed=False)
    if n_comp <= 1:
        return walk, n_comp_before, 1.0, 0

    counts = np.bincount(labels)
    largest_label = int(np.argmax(counts))
    keep_names = {nodes[i] for i, lab in enumerate(labels) if lab == largest_label}
    removed = len(walk.real_nodes) - len(keep_names)

    new_compressed = {(a, b) for a, b in walk.compressed_edges if a in keep_names and b in keep_names}
    new_pairs = [(a, b) for a, b in walk.path_pairs if a in keep_names and b in keep_names]
    new_hops = {
        k: v for k, v in walk.pair_min_real_hops.items() if k[0] in keep_names and k[1] in keep_names
    }

    pruned = SignalingBipartiteWalkResult(
        path_pairs=new_pairs,
        real_nodes=keep_names,
        virtual_nodes=set(),
        gin_edges=set(),
        compressed_edges=new_compressed,
        pair_min_real_hops=new_hops,
        virtual_nodes_visited=walk.virtual_nodes_visited,
        status=walk.status,
    )
    frac = len(keep_names) / max(len(walk.real_nodes), 1)
    return pruned, n_comp_before, frac, removed


def bridge_seed_components_gin(
    walk: SignalingBipartiteWalkResult,
    seed_gin_nodes: list[str],
    adj_dir: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    *,
    bridge_hops: int = 12,
) -> tuple[SignalingBipartiteWalkResult, int]:
    if len(seed_gin_nodes) < 2:
        return walk, 0

    n_comp, labels, nodes = _component_labels(walk)
    if n_comp <= 1:
        return walk, 0

    name_to_comp = {nodes[i]: int(labels[i]) for i in range(len(nodes))}
    seed_in_walk = [s for s in seed_gin_nodes if s in walk.real_nodes]
    if len(seed_in_walk) < 2:
        return walk, 0

    comps: dict[int, list[str]] = {}
    for s in seed_in_walk:
        comps.setdefault(name_to_comp.get(s, -1), []).append(s)
    if len(comps) <= 1:
        return walk, 0

    added = 0
    comp_ids = list(comps.keys())
    anchor = comps[comp_ids[0]][0]

    for cid in comp_ids[1:]:
        target = comps[cid][0]
        before = len(walk.compressed_edges)
        if add_bipartite_path_edges(
            walk, adj_dir, index_to_name, name_to_index, anchor, target, bridge_hops
        ):
            added += len(walk.compressed_edges) - before

    return walk, added


def count_seeds_in_gin_real(seed_gin_nodes: list[str], gin_real: set[str]) -> int:
    return sum(1 for s in seed_gin_nodes if s in gin_real)


def repair_connectivity(
    walk: SignalingBipartiteWalkResult,
    seed_gin_nodes: list[str],
    adj_dir: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    *,
    enable_l2: bool = True,
    bridge_hops: int = 12,
    n_seeds_no_loop: int = 0,
) -> tuple[ConnectivityResult, ConnectivityQC]:
    n_comp_before, _, _ = _component_labels(walk)

    if enable_l2:
        walk, bridge_added = bridge_seed_components_gin(
            walk, seed_gin_nodes, adj_dir, index_to_name, name_to_index, bridge_hops=bridge_hops
        )
    else:
        bridge_added = 0

    pruned, _, lcc_frac, removed = largest_connected_subgraph(walk)
    rebuild_gin_edges_from_compressed(pruned, adj_dir, index_to_name, name_to_index)
    n_comp_after, _, _ = _component_labels(pruned)
    seeds_in = count_seeds_in_gin_real(seed_gin_nodes, pruned.real_nodes)

    result = ConnectivityResult(
        walk=pruned,
        n_components_before=n_comp_before,
        n_components_after=n_comp_after,
        largest_component_frac=lcc_frac,
        nodes_removed=removed,
        bridge_edges_added=bridge_added,
        seeds_in_real_nodes=seeds_in,
        seeds_total=len(seed_gin_nodes),
        seed_gin_nodes=list(seed_gin_nodes),
    )
    qc = ConnectivityQC(
        n_components_before=n_comp_before,
        n_components_after=n_comp_after,
        largest_component_frac=lcc_frac,
        nodes_removed_l1=removed,
        bridge_edges_added=bridge_added,
        seeds_in_real_nodes=seeds_in,
        seeds_total=len(seed_gin_nodes),
        node_stats=node_count_stats(pruned, n_seeds_no_loop=n_seeds_no_loop),
    )
    return result, qc

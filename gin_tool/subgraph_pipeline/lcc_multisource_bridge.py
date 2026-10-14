"""LCC-anchored island bridge via shortest-path DAG (context subgraph post-step)."""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field

import networkx as nx
from scipy import sparse

from gin_tool.io import IndexedEdges, directed_csr_from_indexed


@dataclass
class LccBridgeResult:
    edges: list[tuple[str, str, float]]
    n_cc_before: int
    n_cc_after: int
    main_size_before: int
    main_size_after: int
    bridge_edges_added: int
    bridge_nodes_added: int
    satellite_nodes_scanned: int
    satellite_nodes_reached: int
    satellite_components_before: int
    max_hops: int
    hop_kind: str = "graph-hop"
    mode: str = "shortest_path_dag"
    details: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "n_cc_before": self.n_cc_before,
            "n_cc_after": self.n_cc_after,
            "main_size_before": self.main_size_before,
            "bridge_edges_added": self.bridge_edges_added,
            "bridge_nodes_added": self.bridge_nodes_added,
            "satellite_nodes_scanned": self.satellite_nodes_scanned,
            "satellite_nodes_reached": self.satellite_nodes_reached,
            "satellite_components_before": self.satellite_components_before,
            "max_hops": self.max_hops,
            "hop_kind": self.hop_kind,
            "mode": self.mode,
            "hop_note": (
                "Each directed GIN edge = 1 graph-hop. "
                "graph-hop=2 ≈ one real-hop (R→V→R) on strict bipartite GIN; "
                "this is NOT real-hop BFS."
            ),
            **self.details,
        }


def _undirected_components(edges: list[tuple[str, str, float]]) -> list[set[str]]:
    g = nx.Graph()
    for a, b, _ in edges:
        g.add_edge(a, b)
    return sorted(nx.connected_components(g), key=len, reverse=True)


def _edge_weight(adj: sparse.csr_matrix, u: int, v: int) -> float:
    return float(adj[u, v])


def _multisource_bfs(
    adj: sparse.csr_matrix,
    sources: list[int],
    max_hops: int,
) -> tuple[dict[int, int], dict[int, set[int]]]:
    dist: dict[int, int] = {}
    pred: dict[int, set[int]] = defaultdict(set)
    q: deque[int] = deque()
    for s in sources:
        if s not in dist:
            dist[s] = 0
            q.append(s)
    while q:
        u = q.popleft()
        du = dist[u]
        if du >= max_hops:
            continue
        lo, hi = adj.indptr[u], adj.indptr[u + 1]
        for v in adj.indices[lo:hi]:
            v = int(v)
            nd = du + 1
            if v not in dist:
                dist[v] = nd
                pred[v].add(u)
                q.append(v)
            elif dist[v] == nd:
                pred[v].add(u)
    return dist, pred


def _shortest_path_dag_edges(
    target: int,
    dist: dict[int, int],
    pred: dict[int, set[int]],
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    max_hops: int,
) -> set[tuple[str, str, float]]:
    """All edges on any shortest path from main to target (equal-hop ties kept)."""
    if target not in dist or dist[target] > max_hops:
        return set()

    on_path: set[int] = {target}
    frontier = [target]
    while frontier:
        nxt: list[int] = []
        for v in frontier:
            for u in pred.get(v, ()):
                if u not in on_path:
                    on_path.add(u)
                    if dist[u] > 0:
                        nxt.append(u)
        frontier = nxt

    out: set[tuple[str, str, float]] = set()
    for v in on_path:
        for u in pred.get(v, ()):
            if u in on_path:
                out.add((index_to_name[u], index_to_name[v], _edge_weight(adj, u, v)))
    return out


def lcc_multisource_bridge(
    edges: list[tuple[str, str, float]],
    adj_dir: sparse.csr_matrix,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    *,
    max_hops: int = 2,
) -> LccBridgeResult:
    """
    Largest-CC island repair: multi-source BFS from main, add shortest-path DAG
    edges to each reachable satellite node within ``max_hops`` graph-hops.
    """
    comps = _undirected_components(edges)
    if len(comps) <= 1:
        return LccBridgeResult(
            edges=list(edges),
            n_cc_before=1,
            n_cc_after=1,
            main_size_before=len(comps[0]),
            main_size_after=len(comps[0]),
            bridge_edges_added=0,
            bridge_nodes_added=0,
            satellite_nodes_scanned=0,
            satellite_nodes_reached=0,
            satellite_components_before=0,
            max_hops=max_hops,
        )

    main = comps[0]
    satellite_nodes: set[str] = set()
    for comp in comps[1:]:
        satellite_nodes.update(comp)

    sources = [name_to_index[n] for n in main if n in name_to_index]
    satellite_indices = [name_to_index[n] for n in satellite_nodes if n in name_to_index]
    if not sources or not satellite_indices:
        raise ValueError("main or satellite has no nodes indexed in GIN")

    dist, pred = _multisource_bfs(adj_dir, sources, max_hops)
    bridge_edges: set[tuple[str, str, float]] = set()
    reached = 0
    for t_idx in satellite_indices:
        path_edges = _shortest_path_dag_edges(
            t_idx, dist, pred, adj_dir, index_to_name, max_hops
        )
        if path_edges:
            reached += 1
            bridge_edges.update(path_edges)

    existing_undir = {tuple(sorted((a, b))) for a, b, _ in edges}
    merged = list(edges)
    added_count = 0
    added_nodes: set[str] = set()
    for a, b, w in sorted(bridge_edges):
        key = tuple(sorted((a, b)))
        if key not in existing_undir:
            merged.append((a, b, w))
            existing_undir.add(key)
            added_count += 1
            added_nodes.update((a, b))

    comps_after = _undirected_components(merged)
    return LccBridgeResult(
        edges=merged,
        n_cc_before=len(comps),
        n_cc_after=len(comps_after),
        main_size_before=len(main),
        main_size_after=len(comps_after[0]),
        bridge_edges_added=added_count,
        bridge_nodes_added=len(added_nodes - main),
        satellite_nodes_scanned=len(satellite_indices),
        satellite_nodes_reached=reached,
        satellite_components_before=len(comps) - 1,
        max_hops=max_hops,
        details={
            "island_sizes_before": [len(c) for c in comps[1:]],
            "island_sizes_after": [len(c) for c in comps_after[1:]],
        },
    )


def apply_lcc_bridge_to_edges(
    edges: list[tuple[str, str, float]],
    ed: IndexedEdges,
    *,
    max_hops: int = 2,
) -> LccBridgeResult:
    adj = directed_csr_from_indexed(ed)
    return lcc_multisource_bridge(
        edges,
        adj,
        ed.index_to_name,
        ed.name_to_index,
        max_hops=max_hops,
    )

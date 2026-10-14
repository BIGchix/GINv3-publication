"""Context-aware subgraph: active TF regulatory edges + mechanistic BFS + reverse prune."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from gin_tool.io import (
    IndexedEdges,
    build_indexed_edges,
    directed_csr_from_indexed,
    induced_subgraph_csr,
    symmetrize_undirected_unweighted,
)
from gin_tool.regulatory_filter import (
    RegulatoryPrefilterStats,
    allow_directed_edge,
    filter_indexed_edges_regulatory,
)
from gin_tool.subgraph import seeds_from_names


def _out_neighbors(adj: sparse.csr_matrix, u: int) -> np.ndarray:
    lo, hi = adj.indptr[u], adj.indptr[u + 1]
    return adj.indices[lo:hi]


def _is_promoter_seed_name(name: str) -> bool:
    return ".promoter" in name.lower()


def _resolve_start_indices(
    seed_gin_names: list[str],
    active_tf_gin_names: list[str],
    name_to_index: dict[str, int],
) -> list[int]:
    starts: list[int] = []
    seen: set[int] = set()
    for name in list(seed_gin_names) + list(active_tf_gin_names):
        if name not in name_to_index:
            continue
        idx = name_to_index[name]
        if idx not in seen:
            seen.add(idx)
            starts.append(idx)
    return starts


def directed_expand(
    adj_dir: sparse.csr_matrix,
    index_to_name: list[str],
    start_indices: list[int],
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    max_depth: int,
    max_nodes: int | None,
) -> set[int]:
    """Forward BFS on directed CSR with regulatory edge mask."""
    visited: set[int] = set(start_indices)
    frontier: deque[tuple[int, int]] = deque((i, 0) for i in start_indices)

    while frontier:
        u, depth = frontier.popleft()
        if depth >= max_depth:
            continue
        u_name = index_to_name[u]
        for v in _out_neighbors(adj_dir, u):
            v = int(v)
            v_name = index_to_name[v]
            if not allow_directed_edge(
                u_name,
                v_name,
                active_tf_names=active_tf_names,
                seed_symbols=seed_symbols,
                seed_gin_names=seed_gin_names,
            ):
                continue
            if v not in visited:
                if max_nodes is not None and len(visited) >= max_nodes:
                    return visited
                visited.add(v)
                frontier.append((v, depth + 1))
    return visited


def directed_expand_reverse(
    adj_rev: sparse.csr_matrix,
    index_to_name: list[str],
    start_indices: list[int],
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    max_depth: int,
    max_nodes: int | None,
    mask_edges: bool,
) -> set[int]:
    """
    BFS along **incoming** GIN edges (walk on ``adj_dir.T``).

    Used for initial seed promoters so graph-hop=1 reaches regulating TFs.
    Edge policy is still evaluated on the **forward** arc ``pred → current``.
    """
    visited: set[int] = set(start_indices)
    frontier: deque[tuple[int, int]] = deque((i, 0) for i in start_indices)

    while frontier:
        u, depth = frontier.popleft()
        if depth >= max_depth:
            continue
        u_name = index_to_name[u]
        for pred in _out_neighbors(adj_rev, u):
            pred = int(pred)
            pred_name = index_to_name[pred]
            if mask_edges and not allow_directed_edge(
                pred_name,
                u_name,
                active_tf_names=active_tf_names,
                seed_symbols=seed_symbols,
                seed_gin_names=seed_gin_names,
            ):
                continue
            if pred not in visited:
                if max_nodes is not None and len(visited) >= max_nodes:
                    return visited
                visited.add(pred)
                frontier.append((pred, depth + 1))
    return visited


def reverse_reachable_to_anchors(
    adj_dir: sparse.csr_matrix,
    universe: set[int],
    anchor_indices: set[int],
) -> set[int]:
    """
    Keep nodes in ``universe`` that are **ancestors** of seed anchors
    (including the anchors themselves):

      path node → … → anchor

    Drops forward-only sinks that cannot flow back into any seed (e.g. seed
    protein → chr* dead ends). Upstream nodes from promoter reverse-walk
    (TF → promoter;state → seed.promoter.active) are retained because they
    reach a promoter anchor.

    Historical note: a brief intermediate version kept ancestors ∪ descendants;
    product choice is ancestor-only.
    """
    adj_rev = adj_dir.T.tocsr()
    keep: set[int] = {a for a in anchor_indices if a in universe}
    frontier: set[int] = set(keep)
    while frontier:
        nxt: set[int] = set()
        for v in frontier:
            for u in _out_neighbors(adj_rev, v):
                u = int(u)
                if u in universe and u not in keep:
                    keep.add(u)
                    nxt.add(u)
        frontier = nxt
    return keep


@dataclass
class ContextSubgraphStats:
    nodes_before_prune: int
    nodes_after_prune: int
    anchor_count: int
    active_tf_count: int
    mech_max_depth: int
    max_nodes_cap: int | None
    seed_only_anchors: bool = True
    prefilter: RegulatoryPrefilterStats | None = None
    n_promoter_reverse_starts: int = 0
    n_forward_starts: int = 0
    promoter_reverse_walk: bool = True
    prune_mode: str = "ancestors"

    def to_dict(self) -> dict:
        out = {
            "nodes_before_prune": self.nodes_before_prune,
            "nodes_after_prune": self.nodes_after_prune,
            "anchor_count": self.anchor_count,
            "active_tf_count": self.active_tf_count,
            "mech_max_depth": self.mech_max_depth,
            "max_nodes_cap": self.max_nodes_cap,
            "seed_only_anchors": self.seed_only_anchors,
            "promoter_reverse_walk": self.promoter_reverse_walk,
            "n_promoter_reverse_starts": self.n_promoter_reverse_starts,
            "n_forward_starts": self.n_forward_starts,
            "prune_mode": self.prune_mode,
        }
        if self.prefilter is not None:
            out["prefilter"] = self.prefilter.to_dict()
        return out


def extract_context_region(
    ed: IndexedEdges,
    seed_gin_names: list[str],
    active_tf_gin_names: list[str],
    seed_symbols: set[str],
    *,
    mech_max_depth: int = 2,
    max_nodes: int | None = 100_000,
    active_tf_symbols: set[str] | None = None,
    tf_pool: set[str] | None = None,
    include_active_tf_anchors: bool = False,
    promoter_reverse_walk: bool = True,
    n_jobs: int = 1,
) -> tuple[
    sparse.csr_matrix,
    sparse.csr_matrix,
    dict[int, int],
    dict[int, int],
    ContextSubgraphStats,
]:
    """
    Build context subgraph on ``ed`` (already loaded).

    Seed **promoters** optionally expand along reverse directed edges (to reach
    regulating TFs at small hop). All other seed anchors expand forward. When a
    non-promoter seed later reaches a promoter, subsequent expansion from that
    visit remains forward (promoters are not reclassified as reverse seeds).

    ``n_jobs``: parallel per-seed BFS when ``>1`` (thread pool; shared CSR).
    With parallel expansion, ``max_nodes`` is applied after merging seed balls
    (deterministic trim) rather than mid-BFS.
    """
    from gin_tool.subgraph_pipeline.parallel_util import map_threads, resolve_n_jobs

    name_to_index = ed.name_to_index
    index_to_name = ed.index_to_name
    seed_gin_set = set(seed_gin_names)
    active_tf_set = set(active_tf_gin_names)
    if active_tf_symbols:
        active_tf_set |= set(active_tf_symbols)
    tf_pool = tf_pool or set()
    jobs = resolve_n_jobs(n_jobs)

    prefilter_stats: RegulatoryPrefilterStats | None = None
    if include_active_tf_anchors:
        work_ed = ed
        starts = _resolve_start_indices(seed_gin_names, active_tf_gin_names, name_to_index)
        mask_during_bfs = True
    else:
        work_ed, prefilter_stats = filter_indexed_edges_regulatory(
            ed,
            active_tf_names=active_tf_set,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_set,
            tf_pool=tf_pool,
        )
        starts = _resolve_start_indices(seed_gin_names, [], name_to_index)
        mask_during_bfs = False

    if not starts:
        raise ValueError("No seed nodes found in graph index")

    adj_dir = directed_csr_from_indexed(work_ed)
    adj_rev = adj_dir.T.tocsr()
    sym = symmetrize_undirected_unweighted(adj_dir)
    anchor_indices = set(starts)

    expand_kwargs = (
        dict(
            active_tf_names=active_tf_set,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_set,
        )
        if mask_during_bfs
        else dict(
            active_tf_names=set(),
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_set,
        )
    )

    promoter_starts: list[int] = []
    forward_starts: list[int] = []
    if promoter_reverse_walk:
        for i in starts:
            if _is_promoter_seed_name(index_to_name[i]):
                promoter_starts.append(i)
            else:
                forward_starts.append(i)
    else:
        forward_starts = list(starts)

    def _cap(nodes: set[int]) -> set[int]:
        if max_nodes is None or len(nodes) <= max_nodes:
            return nodes
        keep = set(starts)
        for i in sorted(nodes - keep):
            if len(keep) >= max_nodes:
                break
            keep.add(i)
        return keep

    def _fwd_one(s: int) -> set[int]:
        if mask_during_bfs:
            return directed_expand(
                adj_dir,
                index_to_name,
                [s],
                active_tf_names=expand_kwargs["active_tf_names"],
                seed_symbols=seed_symbols,
                seed_gin_names=seed_gin_set,
                max_depth=mech_max_depth,
                max_nodes=None,
            )
        return _directed_expand_unmasked(
            adj_dir,
            index_to_name,
            [s],
            max_depth=mech_max_depth,
            max_nodes=None,
        )

    def _rev_one(s: int) -> set[int]:
        return directed_expand_reverse(
            adj_rev,
            index_to_name,
            [s],
            active_tf_names=expand_kwargs["active_tf_names"],
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_set,
            max_depth=mech_max_depth,
            max_nodes=None,
            mask_edges=mask_during_bfs,
        )

    nodes: set[int] = set(starts)
    if forward_starts:
        if jobs <= 1 or len(forward_starts) == 1:
            if mask_during_bfs:
                nodes |= directed_expand(
                    adj_dir,
                    index_to_name,
                    forward_starts,
                    active_tf_names=expand_kwargs["active_tf_names"],
                    seed_symbols=seed_symbols,
                    seed_gin_names=seed_gin_set,
                    max_depth=mech_max_depth,
                    max_nodes=max_nodes,
                )
            else:
                nodes |= _directed_expand_unmasked(
                    adj_dir,
                    index_to_name,
                    forward_starts,
                    max_depth=mech_max_depth,
                    max_nodes=max_nodes,
                )
        else:
            parts = map_threads(_fwd_one, forward_starts, n_jobs=jobs)
            for p in parts:
                nodes |= p
            nodes = _cap(nodes)
    if promoter_starts:
        if jobs <= 1 or len(promoter_starts) == 1:
            nodes |= directed_expand_reverse(
                adj_rev,
                index_to_name,
                promoter_starts,
                active_tf_names=expand_kwargs["active_tf_names"],
                seed_symbols=seed_symbols,
                seed_gin_names=seed_gin_set,
                max_depth=mech_max_depth,
                max_nodes=max_nodes,
                mask_edges=mask_during_bfs,
            )
        else:
            parts = map_threads(_rev_one, promoter_starts, n_jobs=jobs)
            for p in parts:
                nodes |= p
            nodes = _cap(nodes)

    before = len(nodes)
    nodes = reverse_reachable_to_anchors(adj_dir, nodes, anchor_indices)
    after = len(nodes)

    sub_sym, old_to_new, new_to_old = induced_subgraph_csr(sym, nodes)
    sub_dir, _, _ = induced_subgraph_csr(adj_dir, nodes)
    stats = ContextSubgraphStats(
        nodes_before_prune=before,
        nodes_after_prune=after,
        anchor_count=len(anchor_indices),
        active_tf_count=len(active_tf_set),
        mech_max_depth=mech_max_depth,
        max_nodes_cap=max_nodes,
        seed_only_anchors=not include_active_tf_anchors,
        prefilter=prefilter_stats,
        n_promoter_reverse_starts=len(promoter_starts),
        n_forward_starts=len(forward_starts),
        promoter_reverse_walk=promoter_reverse_walk,
        prune_mode="ancestors",
    )
    return sub_sym, sub_dir, old_to_new, new_to_old, stats


def _directed_expand_unmasked(
    adj_dir: sparse.csr_matrix,
    index_to_name: list[str],
    start_indices: list[int],
    *,
    max_depth: int,
    max_nodes: int | None,
) -> set[int]:
    """Forward BFS without per-edge mask (graph already pre-filtered)."""
    visited: set[int] = set(start_indices)
    frontier: deque[tuple[int, int]] = deque((i, 0) for i in start_indices)

    while frontier:
        u, depth = frontier.popleft()
        if depth >= max_depth:
            continue
        for v in _out_neighbors(adj_dir, u):
            v = int(v)
            if v not in visited:
                if max_nodes is not None and len(visited) >= max_nodes:
                    return visited
                visited.add(v)
                frontier.append((v, depth + 1))
    return visited


def extract_context_region_from_sifi(
    sifi_path: str,
    seed_gin_names: list[str],
    active_tf_gin_names: list[str],
    seed_symbols: set[str],
    *,
    max_lines: int | None = None,
    mech_max_depth: int = 2,
    max_nodes: int | None = 100_000,
    active_tf_symbols: set[str] | None = None,
    tf_pool: set[str] | None = None,
    include_active_tf_anchors: bool = False,
    promoter_reverse_walk: bool = True,
    n_jobs: int = 1,
) -> tuple[
    sparse.csr_matrix,
    sparse.csr_matrix,
    dict[int, int],
    dict[int, int],
    IndexedEdges,
    ContextSubgraphStats,
]:
    """Load ``.sifi`` then ``extract_context_region``."""
    ed = build_indexed_edges(sifi_path, max_lines=max_lines)
    seeds_from_names(ed.name_to_index, seed_gin_names)  # validate
    sub_sym, sub_dir, o2n, n2o, stats = extract_context_region(
        ed,
        list(seed_gin_names),
        list(active_tf_gin_names),
        seed_symbols,
        mech_max_depth=mech_max_depth,
        max_nodes=max_nodes,
        active_tf_symbols=active_tf_symbols,
        tf_pool=tf_pool,
        include_active_tf_anchors=include_active_tf_anchors,
        promoter_reverse_walk=promoter_reverse_walk,
        n_jobs=n_jobs,
    )
    return sub_sym, sub_dir, o2n, n2o, ed, stats

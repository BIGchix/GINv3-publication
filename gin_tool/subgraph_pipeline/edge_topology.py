"""Classify compressed real edges vs scBionet six-layer forward chain."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Iterable

from gin_tool.subgraph_pipeline.manifest import EdgeTopologyQC, NodeCountStats
from gin_tool.subgraph_pipeline.signaling_bipartite import SignalingBipartiteWalkResult
from gin_tool.subgraph_pipeline.stats import node_count_stats

LAYER_RNA = 0
LAYER_PROTEIN = 1
LAYER_COMPLEX = 2
LAYER_METAB = 3
LAYER_TF = 4
LAYER_PROMOTER = 5
LAYER_OTHER = 6

LAYER_NAMES = {
    LAYER_RNA: "rna",
    LAYER_PROTEIN: "protein/subunit",
    LAYER_COMPLEX: "protein_complex",
    LAYER_METAB: "complex_metab",
    LAYER_TF: "tf",
    LAYER_PROMOTER: "promoter",
    LAYER_OTHER: "other",
}

CANONICAL_FORWARD_PAIRS = {
    (LAYER_RNA, LAYER_PROTEIN),
    (LAYER_PROTEIN, LAYER_COMPLEX),
    (LAYER_COMPLEX, LAYER_METAB),
    (LAYER_METAB, LAYER_TF),
    (LAYER_TF, LAYER_PROMOTER),
    (LAYER_PROMOTER, LAYER_RNA),
}


@dataclass
class EdgeTopologyReport:
    n_edges: int = 0
    canonical_forward: int = 0
    cross_layer_forward: int = 0
    reverse: int = 0
    lateral: int = 0
    involving_other: int = 0
    by_layer_pair: dict[str, int] = field(default_factory=dict)
    samples: dict[str, list[list[str]]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def assign_scbionet_layer(node: str, *, tf_pool: set[str] | None = None) -> int:
    tf_pool = tf_pool or set()
    lower = node.lower()
    if lower.endswith(".rna"):
        return LAYER_RNA
    if ".promoter" in lower and ";" not in node.split(".promoter", 1)[0]:
        return LAYER_PROMOTER
    if "metabol" in lower or lower.endswith(".met"):
        return LAYER_METAB
    sym = node.split(".")[0] if "." in node else node
    if sym in tf_pool:
        return LAYER_TF
    if lower.endswith(".protein"):
        return LAYER_PROTEIN
    if "." not in node:
        return LAYER_TF if sym in tf_pool else LAYER_PROTEIN
    if "complex" in lower:
        return LAYER_COMPLEX
    return LAYER_OTHER


def classify_compressed_edges(
    edges: Iterable[tuple[str, str]],
    *,
    tf_pool: set[str] | None = None,
    sample_n: int = 15,
) -> EdgeTopologyReport:
    report = EdgeTopologyReport()
    for a, b in edges:
        la = assign_scbionet_layer(a, tf_pool=tf_pool)
        lb = assign_scbionet_layer(b, tf_pool=tf_pool)
        key = f"{LAYER_NAMES.get(la, la)}->{LAYER_NAMES.get(lb, lb)}"
        report.by_layer_pair[key] = report.by_layer_pair.get(key, 0) + 1
        report.n_edges += 1

        if LAYER_OTHER in (la, lb):
            report.involving_other += 1
            bucket = "involving_other"
        elif la == lb:
            report.lateral += 1
            bucket = "lateral"
        elif lb < la:
            report.reverse += 1
            bucket = "reverse"
        elif (la, lb) in CANONICAL_FORWARD_PAIRS:
            report.canonical_forward += 1
            bucket = "canonical_forward"
        elif lb > la:
            report.cross_layer_forward += 1
            bucket = "cross_layer_forward"
        else:
            report.involving_other += 1
            bucket = "involving_other"

        if bucket not in report.samples:
            report.samples[bucket] = []
        if len(report.samples[bucket]) < sample_n:
            report.samples[bucket].append([a, b])

    return report


def is_canonical_forward_edge(a: str, b: str, *, tf_pool: set[str] | None = None) -> bool:
    la = assign_scbionet_layer(a, tf_pool=tf_pool)
    lb = assign_scbionet_layer(b, tf_pool=tf_pool)
    return (la, lb) in CANONICAL_FORWARD_PAIRS


def prune_walk_to_canonical_forward(
    walk: SignalingBipartiteWalkResult,
    *,
    protect_nodes: Iterable[str],
    tf_pool: set[str] | None = None,
) -> tuple[SignalingBipartiteWalkResult, dict[str, Any]]:
    protect = {n for n in protect_nodes if n}
    edges_in = set(walk.compressed_edges)
    kept_edges = {(a, b) for a, b in edges_in if is_canonical_forward_edge(a, b, tf_pool=tf_pool)}
    dropped_edges = edges_in - kept_edges

    endpoint_nodes = {n for e in kept_edges for n in e}
    nodes_in = set(walk.real_nodes)
    kept_nodes = endpoint_nodes | (protect & nodes_in)
    dropped_nodes = nodes_in - kept_nodes
    kept_edges = {(a, b) for a, b in kept_edges if a in kept_nodes and b in kept_nodes}

    pruned = SignalingBipartiteWalkResult(
        path_pairs=[(a, b) for a, b in walk.path_pairs if a in kept_nodes and b in kept_nodes],
        real_nodes=kept_nodes,
        virtual_nodes=set(),
        gin_edges=set(),
        compressed_edges=kept_edges,
        pair_min_real_hops={
            k: v for k, v in walk.pair_min_real_hops.items() if k[0] in kept_nodes and k[1] in kept_nodes
        },
        virtual_nodes_visited=walk.virtual_nodes_visited,
        status=walk.status,
    )
    stats = {
        "edges_before": len(edges_in),
        "edges_after": len(kept_edges),
        "edges_dropped": len(dropped_edges),
        "nodes_before": len(nodes_in),
        "nodes_after": len(kept_nodes),
        "nodes_dropped": len(dropped_nodes),
        "protect_kept": len(protect & kept_nodes),
    }
    return pruned, stats


def soften_noncanonical_edges(
    walk: SignalingBipartiteWalkResult,
    *,
    protect_nodes: Iterable[str],
    tf_pool: set[str] | None = None,
    max_noncan_fraction: float = 0.25,
) -> tuple[SignalingBipartiteWalkResult, dict[str, Any]]:
    protect = {n for n in protect_nodes if n}
    edges_in = set(walk.compressed_edges)
    nodes_in = set(walk.real_nodes)
    canon = {(a, b) for a, b in edges_in if is_canonical_forward_edge(a, b, tf_pool=tf_pool)}
    noncan = list(edges_in - canon)

    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> bool:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        parent[rb] = ra
        return True

    for a, b in canon:
        union(a, b)
    for n in protect & nodes_in:
        find(n)

    before_comp = len({find(n) for n in protect & nodes_in})

    def merges_protect(a: str, b: str) -> bool:
        find(a)
        find(b)
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        left = any(find(n) == ra for n in protect & nodes_in)
        right = any(find(n) == rb for n in protect & nodes_in)
        return left and right

    noncan.sort(key=lambda e: (0 if merges_protect(e[0], e[1]) else 1, e[0], e[1]))

    kept_noncan: set[tuple[str, str]] = set()
    for a, b in noncan:
        total_if = len(canon) + len(kept_noncan) + 1
        frac = (len(kept_noncan) + 1) / max(total_if, 1)
        useful = merges_protect(a, b)
        if frac > max_noncan_fraction and not useful:
            continue
        if frac > max(max_noncan_fraction, 0.35):
            break
        if useful or frac <= max_noncan_fraction:
            kept_noncan.add((a, b))
            union(a, b)

    kept_edges = set(canon) | kept_noncan
    endpoint_nodes = {n for e in kept_edges for n in e}
    kept_nodes = set(endpoint_nodes) | (protect & nodes_in)

    softened = SignalingBipartiteWalkResult(
        path_pairs=[(a, b) for a, b in walk.path_pairs if a in kept_nodes and b in kept_nodes],
        real_nodes=kept_nodes,
        virtual_nodes=set(),
        gin_edges=set(),
        compressed_edges=kept_edges,
        pair_min_real_hops={
            k: v for k, v in walk.pair_min_real_hops.items() if k[0] in kept_nodes and k[1] in kept_nodes
        },
        virtual_nodes_visited=walk.virtual_nodes_visited,
        status=walk.status,
    )
    after_comp = len({find(n) for n in protect & kept_nodes}) if (protect & kept_nodes) else 0
    stats = {
        "edges_before": len(edges_in),
        "edges_after": len(kept_edges),
        "canon_kept": len(canon),
        "noncan_before": len(noncan),
        "noncan_kept": len(kept_noncan),
        "nodes_before": len(nodes_in),
        "nodes_after": len(kept_nodes),
        "protect_components_before": before_comp,
        "protect_components_after": after_comp,
    }
    return softened, stats


def evaluate_cross_layer_modeling(report: EdgeTopologyReport) -> dict[str, Any]:
    n = max(report.n_edges, 1)
    non_canonical = report.cross_layer_forward + report.reverse + report.lateral
    return {
        "n_edges": report.n_edges,
        "fraction_non_canonical": round(non_canonical / n, 4),
        "fraction_canonical_forward": round(report.canonical_forward / n, 4),
        "needs_modeling_extension": non_canonical > 0,
    }


def apply_edge_topology(
    walk: SignalingBipartiteWalkResult,
    *,
    protect_nodes: Iterable[str],
    tf_pool: set[str] | None = None,
    canonical_edges_only: bool = False,
    soft_noncan: bool = False,
    soft_noncan_max_fraction: float = 0.25,
    n_seeds_no_loop: int = 0,
    adj_dir=None,
    index_to_name: list[str] | None = None,
    name_to_index: dict[str, int] | None = None,
) -> tuple[SignalingBipartiteWalkResult, EdgeTopologyQC, list[str]]:
    from gin_tool.subgraph_pipeline.signaling_bipartite import rebuild_gin_edges_from_compressed

    warnings: list[str] = []
    final_walk = walk

    if canonical_edges_only and soft_noncan:
        warnings.append("both canonical_edges_only and soft_noncan set; preferring soft_noncan")
        canonical_edges_only = False

    if soft_noncan:
        final_walk, stats = soften_noncanonical_edges(
            final_walk,
            protect_nodes=protect_nodes,
            tf_pool=tf_pool,
            max_noncan_fraction=soft_noncan_max_fraction,
        )
        warnings.append(
            f"soft_noncan: edges {stats['edges_before']}→{stats['edges_after']}"
        )
    elif canonical_edges_only:
        final_walk, stats = prune_walk_to_canonical_forward(
            final_walk,
            protect_nodes=protect_nodes,
            tf_pool=tf_pool,
        )
        warnings.append(
            f"canonical_edges_only: edges {stats['edges_before']}→{stats['edges_after']}"
        )

    if adj_dir is not None and index_to_name is not None and name_to_index is not None:
        rebuild_gin_edges_from_compressed(final_walk, adj_dir, index_to_name, name_to_index)

    report = classify_compressed_edges(final_walk.compressed_edges, tf_pool=tf_pool)
    modeling = evaluate_cross_layer_modeling(report)

    qc = EdgeTopologyQC(
        n_edges=report.n_edges,
        canonical_forward=report.canonical_forward,
        cross_layer_forward=report.cross_layer_forward,
        reverse=report.reverse,
        lateral=report.lateral,
        involving_other=report.involving_other,
        by_layer_pair=report.by_layer_pair,
        modeling=modeling,
        node_stats=node_count_stats(final_walk, n_seeds_no_loop=n_seeds_no_loop),
    )
    if report.cross_layer_forward or report.reverse or report.lateral:
        warnings.append(
            f"non-canonical edges: cross={report.cross_layer_forward} "
            f"reverse={report.reverse} lateral={report.lateral}"
        )
    return final_walk, qc, warnings

"""Train-mode regulatory edge filtering (from scBionet, unchanged semantics)."""

from __future__ import annotations

import numpy as np


def keep_regulatory_edge_train(
    u_name: str,
    v_name: str,
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    tf_pool: set[str],
    relaxed: bool = True,
) -> bool:
    from gin_tool.active_tf import gin_node_to_gene_symbol
    from gin_tool.regulatory_filter import is_protein_tf_to_dna_edge

    if not is_protein_tf_to_dna_edge(u_name, v_name, tf_pool):
        return True

    u_sym = gin_node_to_gene_symbol(u_name) or u_name
    if u_sym not in active_tf_names and u_name not in active_tf_names:
        return False

    lower = v_name.lower()
    if ".promoter" in lower:
        gene = v_name.split(".promoter", 1)[0]
        if gene and ";" not in gene:
            return gene in seed_symbols
        return v_name in seed_gin_names
    return True if relaxed else False


def filter_indexed_edges_regulatory_train(
    ed,
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    tf_pool: set[str],
    relaxed: bool = True,
):
    from gin_tool.io import IndexedEdges
    from gin_tool.regulatory_filter import RegulatoryPrefilterStats

    n = len(ed.sources)
    if n == 0:
        stats = RegulatoryPrefilterStats(0, 0, 0, seed_only_anchors=not relaxed)
        return ed, stats

    keep = np.ones(n, dtype=bool)
    index_to_name = ed.index_to_name
    for i in range(n):
        u_name = index_to_name[int(ed.sources[i])]
        v_name = index_to_name[int(ed.targets[i])]
        if not keep_regulatory_edge_train(
            u_name,
            v_name,
            active_tf_names=active_tf_names,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_names,
            tf_pool=tf_pool,
            relaxed=relaxed,
        ):
            keep[i] = False

    kept_n = int(keep.sum())
    stats = RegulatoryPrefilterStats(
        edges_before=n,
        edges_after=kept_n,
        edges_removed=n - kept_n,
        seed_only_anchors=not relaxed,
    )
    if kept_n == n:
        return ed, stats

    return (
        IndexedEdges(
            num_nodes=ed.num_nodes,
            sources=ed.sources[keep],
            targets=ed.targets[keep],
            weights=ed.weights[keep],
            index_to_name=ed.index_to_name,
            name_to_index=ed.name_to_index,
        ),
        stats,
    )

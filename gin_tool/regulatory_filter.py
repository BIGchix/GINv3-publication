"""Pre-filter protein/TF → DNA regulatory edges (active TF only, seeds-only anchors)."""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from gin_tool.active_tf import gin_node_to_gene_symbol, promoter_target_to_gene
from gin_tool.io import IndexedEdges
from gin_tool.mapping import is_virtual_gin_node

_CHR_PREFIX = re.compile(r"^chr\d", re.IGNORECASE)


def _node_symbol(name: str) -> str:
    return gin_node_to_gene_symbol(name) or name


def allow_directed_edge(
    u_name: str,
    v_name: str,
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
) -> bool:
    """
    Edge filter during expansion or pre-filter refinement.

    Regulatory edges (TF → promoter) only from active TFs to seed-gene promoters.
    All other directed edges are allowed (mechanistic / cross-layer).
    """
    gene = promoter_target_to_gene(v_name)
    if gene is not None or ".promoter" in v_name.lower():
        target_gene = gene or gin_node_to_gene_symbol(v_name)
        u_sym = _node_symbol(u_name)
        if u_sym in active_tf_names or u_name in active_tf_names:
            return target_gene in seed_symbols if target_gene else v_name in seed_gin_names
        if ".promoter" in v_name.lower():
            return False
    return True


def is_dna_element(name: str) -> bool:
    """Promoter transitions, chromatin loci, and related virtual DNA states."""
    if promoter_target_to_gene(name) is not None:
        return True
    lower = name.lower()
    if ".promoter" in lower:
        return True
    if is_virtual_gin_node(name) and (";state" in lower or "promoter" in lower):
        return True
    if name.startswith("chr") or _CHR_PREFIX.match(name):
        return True
    return False


def is_protein_or_tf_source(name: str, tf_pool: set[str]) -> bool:
    """Real-node protein or TF that can emit a regulatory edge toward DNA."""
    if is_virtual_gin_node(name):
        return False
    lower = name.lower()
    if lower.endswith(".protein"):
        return True
    sym = gin_node_to_gene_symbol(name)
    if sym and sym in tf_pool:
        return True
    if sym is not None and "." not in name:
        return True
    return False


def is_protein_tf_to_dna_edge(u_name: str, v_name: str, tf_pool: set[str]) -> bool:
    return is_protein_or_tf_source(u_name, tf_pool) and is_dna_element(v_name)


def keep_regulatory_edge(
    u_name: str,
    v_name: str,
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    tf_pool: set[str],
) -> bool:
    """
    Return whether to keep ``u_name → v_name`` when pre-filtering the graph.

    Non protein/TF→DNA edges are kept. Protein/TF→DNA requires an active TF
    source; promoter targets still require seed-gene overlap (Context semantics).
    """
    if not is_protein_tf_to_dna_edge(u_name, v_name, tf_pool):
        return True
    u_sym = gin_node_to_gene_symbol(u_name) or u_name
    if u_sym not in active_tf_names and u_name not in active_tf_names:
        return False
    gene = promoter_target_to_gene(v_name)
    if gene is not None or ".promoter" in v_name.lower():
        return allow_directed_edge(
            u_name,
            v_name,
            active_tf_names=active_tf_names,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_names,
        )
    return True


@dataclass
class RegulatoryPrefilterStats:
    edges_before: int
    edges_after: int
    edges_removed: int
    seed_only_anchors: bool

    def to_dict(self) -> dict:
        return {
            "edges_before": self.edges_before,
            "edges_after": self.edges_after,
            "edges_removed": self.edges_removed,
            "seed_only_anchors": self.seed_only_anchors,
        }


def filter_indexed_edges_regulatory(
    ed: IndexedEdges,
    *,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    tf_pool: set[str],
) -> tuple[IndexedEdges, RegulatoryPrefilterStats]:
    """Drop protein/TF→DNA edges except active TF→allowed DNA targets."""
    n = len(ed.sources)
    if n == 0:
        stats = RegulatoryPrefilterStats(0, 0, 0, seed_only_anchors=True)
        return ed, stats

    keep = np.ones(n, dtype=bool)
    index_to_name = ed.index_to_name
    for i in range(n):
        u_name = index_to_name[int(ed.sources[i])]
        v_name = index_to_name[int(ed.targets[i])]
        if not keep_regulatory_edge(
            u_name,
            v_name,
            active_tf_names=active_tf_names,
            seed_symbols=seed_symbols,
            seed_gin_names=seed_gin_names,
            tf_pool=tf_pool,
        ):
            keep[i] = False

    kept_n = int(keep.sum())
    stats = RegulatoryPrefilterStats(
        edges_before=n,
        edges_after=kept_n,
        edges_removed=n - kept_n,
        seed_only_anchors=True,
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

"""Bounded met-enzyme recovery for context subgraph hyperedges (v4 step 3c).

For each ``;``-joined virtual already present in the subgraph, pull missing
``*.protein`` enzyme places and their GIN arcs to/from that virtual. Does not
invent islands: only enzymes that touch an in-subgraph hyperedge in the master
GIN are added.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

from gin_tool.io import IndexedEdges, directed_csr_from_indexed
from gin_tool.mapping import is_virtual_gin_node
from gin_tool.subgraph_pipeline.role_seed_enrichment import proteins_from_hyperedge_name


@dataclass
class MetEnzymeRecoveryQC:
    n_hyperedges_scanned: int = 0
    n_enzymes_referenced: int = 0
    n_enzymes_added: int = 0
    n_edges_added: int = 0
    enzyme_sample: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_hyperedges_scanned": self.n_hyperedges_scanned,
            "n_enzymes_referenced": self.n_enzymes_referenced,
            "n_enzymes_added": self.n_enzymes_added,
            "n_edges_added": self.n_edges_added,
            "enzyme_sample": self.enzyme_sample,
        }


def recover_met_enzymes_on_edges(
    edges: Sequence[tuple[str, str, float]],
    ed: IndexedEdges,
    *,
    max_new_enzymes: int = 500,
    n_jobs: int = 1,
) -> tuple[list[tuple[str, str, float]], MetEnzymeRecoveryQC]:
    """Add missing enzyme nodes/edges referenced by subgraph hyperedges."""
    from gin_tool.subgraph_pipeline.parallel_util import map_threads, resolve_n_jobs

    edge_set: dict[tuple[str, str], float] = {}
    nodes: set[str] = set()
    for u, v, w in edges:
        edge_set.setdefault((u, v), float(w))
        nodes.add(u)
        nodes.add(v)

    qc = MetEnzymeRecoveryQC()
    hyperedges = sorted(
        n for n in nodes if is_virtual_gin_node(n) or ";" in n
    )
    qc.n_hyperedges_scanned = len(hyperedges)

    referenced: set[str] = set()
    for hv in hyperedges:
        referenced |= proteins_from_hyperedge_name(hv)
    qc.n_enzymes_referenced = len(referenced)

    adj = directed_csr_from_indexed(ed)
    adj_rev = adj.T.tocsr()
    name_to = ed.name_to_index
    nodes_f = frozenset(nodes)

    def _named_outs(i: int) -> list[tuple[str, float]]:
        lo = adj.indptr[i]
        hi = adj.indptr[i + 1]
        return [
            (ed.index_to_name[int(j)], float(w))
            for j, w in zip(adj.indices[lo:hi], adj.data[lo:hi])
        ]

    def _named_ins(i: int) -> list[tuple[str, float]]:
        lo = adj_rev.indptr[i]
        hi = adj_rev.indptr[i + 1]
        return [
            (ed.index_to_name[int(j)], float(w))
            for j, w in zip(adj_rev.indices[lo:hi], adj_rev.data[lo:hi])
        ]

    candidates = sorted(
        e for e in referenced if e not in nodes_f and e in name_to
    )

    def _one(enz: str) -> tuple[str, list[tuple[str, str, float]]] | None:
        ei = name_to[enz]
        cand_edges: list[tuple[str, str, float]] = []
        for nbr, w in _named_outs(ei):
            if nbr in nodes_f and (is_virtual_gin_node(nbr) or ";" in nbr):
                cand_edges.append((enz, nbr, w))
        for nbr, w in _named_ins(ei):
            if nbr in nodes_f and (is_virtual_gin_node(nbr) or ";" in nbr):
                cand_edges.append((nbr, enz, w))
        if not cand_edges:
            return None
        return enz, cand_edges

    jobs = resolve_n_jobs(n_jobs)
    results = map_threads(_one, candidates, n_jobs=jobs) if candidates else []

    added_enz: list[str] = []
    for res in results:
        if res is None:
            continue
        enz, cand_edges = res
        if len(added_enz) >= max_new_enzymes:
            break
        added_enz.append(enz)
        nodes.add(enz)
        for u, v, w in cand_edges:
            if (u, v) not in edge_set:
                edge_set[(u, v)] = w
                qc.n_edges_added += 1

    qc.n_enzymes_added = len(added_enz)
    qc.enzyme_sample = added_enz[:40]
    out_edges = [(u, v, w) for (u, v), w in edge_set.items()]
    return out_edges, qc

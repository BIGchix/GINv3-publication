"""Rehydrate compressed real→real signaling edges back to GIN R→V→R paths."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from scipy import sparse

from gin_tool.io import IndexedEdges, directed_csr_from_indexed
from gin_tool.mapping import is_virtual_gin_node


def _out_neighbors(adj: sparse.csr_matrix, u: int) -> np.ndarray:
    lo, hi = adj.indptr[u], adj.indptr[u + 1]
    return adj.indices[lo:hi]


def _edge_weight(adj: sparse.csr_matrix, u: int, v: int) -> float:
    return float(adj[u, v])


@dataclass
class RehydratedEdge:
    """One logical connection between two real GIN nodes."""

    source: str
    target: str
    sign: float
    via_virtual: str | None = None
    synthetic_transition: str | None = None
    ambiguous_virtuals: list[str] = field(default_factory=list)


@dataclass
class RehydratedSubgraph:
    real_nodes: set[str]
    virtual_nodes: set[str]
    gin_edges: list[tuple[str, str, float]]
    logical_edges: list[RehydratedEdge]
    rehydrated_count: int = 0
    synthetic_count: int = 0
    direct_count: int = 0
    ambiguous_count: int = 0

    def to_dict(self) -> dict:
        return {
            "real_nodes": len(self.real_nodes),
            "virtual_nodes": len(self.virtual_nodes),
            "gin_edges": len(self.gin_edges),
            "logical_edges": len(self.logical_edges),
            "rehydrated_count": self.rehydrated_count,
            "synthetic_count": self.synthetic_count,
            "direct_count": self.direct_count,
            "ambiguous_count": self.ambiguous_count,
        }


def find_virtual_intermediates(
    adj: sparse.csr_matrix,
    index_to_name: list[str],
    source_idx: int,
    target_idx: int,
) -> list[tuple[str, float]]:
    """
    Return virtual nodes V with source→V→target, sorted by name for stability.

    Sign is product of edge signs along the two-hop path.
    """
    out: list[tuple[str, float]] = []
    s_lo, s_hi = adj.indptr[source_idx], adj.indptr[source_idx + 1]
    for v_idx in adj.indices[s_lo:s_hi]:
        v_idx = int(v_idx)
        v_name = index_to_name[v_idx]
        if not is_virtual_gin_node(v_name):
            continue
        w1 = _edge_weight(adj, source_idx, v_idx)
        t_lo, t_hi = adj.indptr[v_idx], adj.indptr[v_idx + 1]
        for t_idx in adj.indices[t_lo:t_hi]:
            if int(t_idx) == target_idx:
                w2 = _edge_weight(adj, v_idx, target_idx)
                sign = w1 * w2
                out.append((v_name, sign))
                break
    out.sort(key=lambda x: x[0])
    return out


def rehydrate_compressed_edges(
    compressed_edges: set[tuple[str, str]] | list[tuple[str, str]],
    ed: IndexedEdges,
    *,
    real_nodes: set[str] | None = None,
) -> RehydratedSubgraph:
    """
    For each compressed real→real edge, recover GIN virtual intermediates when possible.

    When no virtual path exists, ``synthetic_transition`` is left ``None`` for the
    exporter to infer (Transcription / Translation / Binding / Catalysis).
    """
    adj = directed_csr_from_indexed(ed)
    name_to_index = ed.name_to_index
    index_to_name = ed.index_to_name

    if real_nodes is None:
        real_nodes = set()
        for a, b in compressed_edges:
            real_nodes.add(a)
            real_nodes.add(b)
    else:
        real_nodes = set(real_nodes)

    logical: list[RehydratedEdge] = []
    virtual_nodes: set[str] = set()
    gin_edge_set: set[tuple[str, str, float]] = set()

    rehydrated = synthetic = direct = ambiguous = 0

    for source, target in sorted(compressed_edges):
        if source not in name_to_index or target not in name_to_index:
            continue
        si, ti = name_to_index[source], name_to_index[target]

        if is_virtual_gin_node(source) or is_virtual_gin_node(target):
            continue

        intermediates = find_virtual_intermediates(adj, index_to_name, si, ti)
        if intermediates:
            via, sign = intermediates[0]
            virtual_nodes.add(via)
            gin_edge_set.add((source, via, _edge_weight(adj, si, name_to_index[via])))
            gin_edge_set.add((via, target, _edge_weight(adj, name_to_index[via], ti)))
            amb = [v for v, _ in intermediates[1:]]
            logical.append(
                RehydratedEdge(
                    source=source,
                    target=target,
                    sign=sign,
                    via_virtual=via,
                    ambiguous_virtuals=amb,
                )
            )
            rehydrated += 1
            if amb:
                ambiguous += 1
        elif abs(_edge_weight(adj, si, ti)) > 0:
            w = _edge_weight(adj, si, ti)
            gin_edge_set.add((source, target, w))
            logical.append(
                RehydratedEdge(source=source, target=target, sign=w, via_virtual=None)
            )
            direct += 1
        else:
            logical.append(
                RehydratedEdge(source=source, target=target, sign=1.0, via_virtual=None)
            )
            synthetic += 1

    return RehydratedSubgraph(
        real_nodes=real_nodes,
        virtual_nodes=virtual_nodes,
        gin_edges=sorted(gin_edge_set),
        logical_edges=logical,
        rehydrated_count=rehydrated,
        synthetic_count=synthetic,
        direct_count=direct,
        ambiguous_count=ambiguous,
    )


def rehydrate_from_signaling_result(
    compressed_edges: set[tuple[str, str]],
    real_nodes: set[str],
    ed: IndexedEdges,
) -> RehydratedSubgraph:
    """Convenience wrapper for ``SignalingWalkResult`` outputs."""
    return rehydrate_compressed_edges(compressed_edges, ed, real_nodes=real_nodes)

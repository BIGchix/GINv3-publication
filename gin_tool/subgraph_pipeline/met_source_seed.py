"""Metabolite source seeding for met↔met digraphs (in-degree 0 / condensation sources)."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from gin_tool.spn_export import _is_metabolite


@dataclass(frozen=True)
class MetSourceSeedResult:
    sources: tuple[str, ...]
    """met↔met nodes with in-degree 0 (have ≥1 out-edge)."""
    n_met_nodes: int
    n_met_edges: int
    has_cycle: bool
    n_seeded_new: int
    n_already_present: int
    default_level: float


def _is_met_place(name: str) -> bool:
    return _is_metabolite(str(name))


def met_met_adjacency(
    edges: Iterable[tuple[str, str, float]],
) -> tuple[dict[str, set[str]], dict[str, set[str]], set[str]]:
    """Return (succ, pred, nodes) restricted to metabolite↔metabolite arcs."""
    succ: dict[str, set[str]] = defaultdict(set)
    pred: dict[str, set[str]] = defaultdict(set)
    nodes: set[str] = set()
    for src, dst, _w in edges:
        if not (_is_met_place(src) and _is_met_place(dst)):
            continue
        succ[src].add(dst)
        pred[dst].add(src)
        nodes.add(src)
        nodes.add(dst)
    return succ, pred, nodes


def metabolite_in_degree_zero_sources(
    edges: Iterable[tuple[str, str, float]],
) -> list[str]:
    """Upstream metabolites: in-degree 0 in the met↔met digraph and out-degree ≥1."""
    succ, pred, nodes = met_met_adjacency(edges)
    return sorted(m for m in nodes if len(pred.get(m, ())) == 0 and len(succ.get(m, ())) > 0)


def metabolite_source_places_in_net(
    edges: Iterable[tuple[str, str, float]],
    net_places: set[str] | frozenset[str],
) -> tuple[str, ...]:
    """In-degree-0 met↔met sources that exist as Petri places in the subgraph net."""
    place_set = set(net_places)
    return tuple(
        p for p in metabolite_in_degree_zero_sources(edges) if p in place_set
    )


def met_graph_has_cycle(edges: Iterable[tuple[str, str, float]]) -> bool:
    succ, _pred, nodes = met_met_adjacency(edges)
    WHITE, GRAY, BLACK = 0, 1, 2
    color: dict[str, int] = {}

    def dfs(u: str) -> bool:
        color[u] = GRAY
        for v in succ.get(u, ()):
            c = color.get(v, WHITE)
            if c == GRAY:
                return True
            if c == WHITE and dfs(v):
                return True
        color[u] = BLACK
        return False

    for u in nodes:
        if color.get(u, WHITE) == WHITE and dfs(u):
            return True
    return False


def load_subgraph_edges_tsv(path: str | Path) -> list[tuple[str, str, float]]:
    rows: list[tuple[str, str, float]] = []
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            if parts[0].lower() in {"src", "source", "u"}:
                continue
            w = float(parts[2]) if len(parts) > 2 and parts[2] not in ("", ".") else 1.0
            rows.append((parts[0], parts[1], w))
    return rows


def merge_met_source_levels(
    base_levels: Mapping[str, float],
    *,
    edges: Iterable[tuple[str, str, float]],
    default_level: float = 10.0,
    only_places: set[str] | None = None,
) -> tuple[dict[str, float], MetSourceSeedResult]:
    """Add in-degree-0 met sources missing from ``base_levels``.

    Existing metabolite seeds keep their levels. New sources get ``default_level``.
    If ``only_places`` is set, skip sources not in that place set (e.g. SPN places).
    """
    sources = metabolite_in_degree_zero_sources(edges)
    if only_places is not None:
        sources = [s for s in sources if s in only_places]
    succ, _pred, nodes = met_met_adjacency(edges)
    n_edges = sum(len(v) for v in succ.values())
    out = {str(k): float(v) for k, v in base_levels.items() if float(v) > 0}
    n_new = 0
    n_have = 0
    for s in sources:
        if s in out:
            n_have += 1
            continue
        out[s] = float(default_level)
        n_new += 1
    stats = MetSourceSeedResult(
        sources=tuple(sources),
        n_met_nodes=len(nodes),
        n_met_edges=n_edges,
        has_cycle=met_graph_has_cycle(edges),
        n_seeded_new=n_new,
        n_already_present=n_have,
        default_level=float(default_level),
    )
    return out, stats

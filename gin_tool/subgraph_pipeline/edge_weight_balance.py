"""Per-place edge-weight balance: Σ|w_in| = Σ|w_out| = min by down-scaling only.

Used as a topology prior when stoichiometric coefficients are missing.
Never inflate weights (avoids token minting when w>1 under Ruths sampling).

Protein-layer balance excludes TF/chromatin regulatory arcs
(``*.promoter;state`` / promoter.active hubs): TFs recycle after DNA binding,
so those edges must not be down-weighted by out-degree balancing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Literal, Sequence

from gin_tool.mapping import is_virtual_gin_node
from gin_tool.spn_export import _is_histone_mark, _is_metabolite, _is_rna, classify_place
from gin_tool.spn_export import SpnPlaceKind

BalanceLayer = Literal["metabolite", "protein", "rna"]


def place_balance_layer(name: str) -> BalanceLayer | None:
    """Map a real place to a balance layer, or None if excluded."""
    if name.startswith("__SPN_sink__") or is_virtual_gin_node(name):
        return None
    # Explicit exclusions: chromatin activity / histone marks
    low = name.lower()
    if ".promoter.active" in low or ".promoter;state" in low:
        return None
    if _is_histone_mark(name):
        return None
    if _is_metabolite(name):
        return "metabolite"
    if _is_rna(name):
        return "rna"
    try:
        kind = classify_place(name)
    except ValueError:
        kind = None
    if kind == SpnPlaceKind.METABOLITE:
        return "metabolite"
    if kind == SpnPlaceKind.RNA:
        return "rna"
    if kind in (SpnPlaceKind.PROTEIN, SpnPlaceKind.PROTEIN_COMPLEX) or name.endswith(
        ".protein"
    ):
        return "protein"
    return None


def is_chromatin_regulatory_virtual(name: str) -> bool:
    """GIN virtual hubs for TF/histone → promoter;state / active (non-consuming biology)."""
    low = name.lower()
    return (
        ".promoter;state" in low
        or low.endswith("promoter;state")
        or ".promoter.active" in low
        or low.endswith("promoter.active")
    )


def _arc_is_chromatin_regulatory(a: str, b: str) -> bool:
    """True if either endpoint is a promoter;state / active virtual hub."""
    return is_chromatin_regulatory_virtual(a) or is_chromatin_regulatory_virtual(b)


@dataclass
class PlaceBalanceRow:
    place: str
    layer: str
    status: str
    n_in: int = 0
    n_out: int = 0
    n_in_excl_chromatin: int = 0
    n_out_excl_chromatin: int = 0
    S_in0: float = 0.0
    S_out0: float = 0.0
    T: float = 0.0
    scale_in: float = 1.0
    scale_out: float = 1.0


@dataclass
class BalanceStats:
    layers: list[str] = field(default_factory=list)
    n_places_considered: int = 0
    n_balanced: int = 0
    n_skipped_missing_side: int = 0
    n_already_equal: int = 0
    n_edges_scaled: int = 0
    n_chromatin_arcs_protected: int = 0
    rows: list[PlaceBalanceRow] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "layers": list(self.layers),
            "n_places_considered": self.n_places_considered,
            "n_balanced": self.n_balanced,
            "n_skipped_missing_side": self.n_skipped_missing_side,
            "n_already_equal": self.n_already_equal,
            "n_edges_scaled": self.n_edges_scaled,
            "n_chromatin_arcs_protected": self.n_chromatin_arcs_protected,
        }


def balance_gin_edge_weights(
    edges: Sequence[tuple[str, str, float]],
    *,
    layers: Iterable[str],
) -> tuple[list[tuple[str, str, float]], BalanceStats]:
    """
    Return a new edge list with per-place in/out weight sums equalized.

    For each selected place P (bipartite place↔virtual arcs):
      S_in  = Σ|w| on virtual→P
      S_out = Σ|w| on P→virtual
      T = min(S_in, S_out) when both > 0
      scale the heavier side by T/S_heavy (never > 1).

    Protein layer: arcs touching ``*.promoter;state`` / promoter.active are
    excluded from Σ and never scaled (TF recycle / regulatory semantics).
    """
    layer_set = {str(x).strip().lower() for x in layers if str(x).strip()}
    layer_set.discard("")
    stats = BalanceStats(layers=sorted(layer_set))
    if not layer_set:
        return [tuple(e) for e in edges], stats

    # index arcs by place endpoint
    in_idx: dict[str, list[int]] = {}
    out_idx: dict[str, list[int]] = {}
    out_edges: list[tuple[str, str, float]] = [tuple(e) for e in edges]  # type: ignore[misc]
    chrom_arc: list[bool] = []

    for i, (a, b, _w) in enumerate(out_edges):
        chrom_arc.append(_arc_is_chromatin_regulatory(a, b))
        a_virt = is_virtual_gin_node(a)
        b_virt = is_virtual_gin_node(b)
        if (not a_virt) and b_virt:
            out_idx.setdefault(a, []).append(i)
        elif a_virt and (not b_virt):
            in_idx.setdefault(b, []).append(i)

    places = sorted(set(in_idx) | set(out_idx))
    scales_in: dict[str, float] = {}
    scales_out: dict[str, float] = {}

    for place in places:
        layer = place_balance_layer(place)
        if layer is None or layer not in layer_set:
            continue
        stats.n_places_considered += 1
        ii = in_idx.get(place, [])
        oo = out_idx.get(place, [])
        # Protein: drop chromatin-regulatory arcs from the balance budget.
        # Metabolite/RNA: keep all incident arcs (promoter hubs are not met/RNA places).
        if layer == "protein":
            ii_use = [i for i in ii if not chrom_arc[i]]
            oo_use = [i for i in oo if not chrom_arc[i]]
        else:
            ii_use, oo_use = ii, oo
        s_in = sum(abs(out_edges[i][2]) for i in ii_use)
        s_out = sum(abs(out_edges[i][2]) for i in oo_use)
        row = PlaceBalanceRow(
            place=place,
            layer=layer,
            status="",
            n_in=len(ii),
            n_out=len(oo),
            n_in_excl_chromatin=len(ii_use),
            n_out_excl_chromatin=len(oo_use),
            S_in0=s_in,
            S_out0=s_out,
        )
        if s_in <= 0.0 or s_out <= 0.0:
            row.status = "unbalanced_missing_side"
            stats.n_skipped_missing_side += 1
            stats.rows.append(row)
            continue
        t = min(s_in, s_out)
        row.T = t
        sc_in = 1.0 if s_in <= t + 1e-15 else t / s_in
        sc_out = 1.0 if s_out <= t + 1e-15 else t / s_out
        # never inflate
        sc_in = min(1.0, sc_in)
        sc_out = min(1.0, sc_out)
        row.scale_in = sc_in
        row.scale_out = sc_out
        if abs(sc_in - 1.0) < 1e-15 and abs(sc_out - 1.0) < 1e-15:
            row.status = "already_equal"
            stats.n_already_equal += 1
        else:
            row.status = "balanced_min"
            stats.n_balanced += 1
            scales_in[place] = sc_in
            scales_out[place] = sc_out
        stats.rows.append(row)

    n_scaled = 0
    n_protected = 0
    new_edges: list[tuple[str, str, float]] = []
    for i, (a, b, w) in enumerate(out_edges):
        ww = float(w)
        protect = chrom_arc[i]
        # Never scale TF/protein ↔ promoter;state (or active) hubs.
        if protect:
            n_protected += 1
            new_edges.append((a, b, ww))
            continue
        if a in scales_out:
            ww *= scales_out[a]
            n_scaled += 1
        if b in scales_in:
            ww *= scales_in[b]
            n_scaled += 1
        new_edges.append((a, b, ww))
    stats.n_edges_scaled = n_scaled
    stats.n_chromatin_arcs_protected = n_protected
    return new_edges, stats


def parse_edge_balance_layers(spec: str | None) -> list[str]:
    """Parse CLI ``metabolite,protein,rna`` / ``off`` / empty."""
    if spec is None:
        return []
    s = str(spec).strip().lower()
    if not s or s in {"off", "none", "0"}:
        return []
    out: list[str] = []
    for part in s.replace(";", ",").split(","):
        p = part.strip()
        if not p:
            continue
        if p in {"met", "mets", "metabolite", "metabolites"}:
            out.append("metabolite")
        elif p in {"protein", "proteins", "prot"}:
            out.append("protein")
        elif p in {"rna", "rnas", "mrna"}:
            out.append("rna")
        else:
            raise ValueError(f"unknown edge-balance layer: {part!r}")
    # stable unique
    seen: set[str] = set()
    uniq: list[str] = []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq

"""First-order degradation (half-life) outlets for RNA, protein, histone, metabolites."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import AbstractSet, Any, Mapping

from gin_tool.spn_export import (
    SpnPetriNet,
    SpnPlaceKind,
    SpnTransitionKind,
    classify_place,
    is_mirna_place,
)
from gin_tool.subgraph_pipeline.spn_layer_rules import (
    SINK_CHROMATIN,
    SINK_HISTONE,
    SINK_METABOLITE,
    SINK_PROTEIN,
    SINK_RNA,
)

# Re-export sinks for callers / tests.
__all__ = [
    "SINK_RNA",
    "SINK_PROTEIN",
    "SINK_METABOLITE",
    "SINK_HISTONE",
    "SINK_CHROMATIN",
    "DEFAULT_EXTRA_DEGRADE_PLACES",
    "DegradationLayerStats",
    "add_molecule_degradation_layer",
    "place_kind_for_degradation",
]

# Hub metabolites that otherwise pile up under maintain (no intrinsic outlet).
# SID315674483 = phosphoric acid (PubChem CID 1004).
# SID315675471 = water (PubChem CID 962).
DEFAULT_EXTRA_DEGRADE_PLACES: dict[str, float] = {
    "SID315674483": 5.0,
    "SID315675471": 5.0,
}


@dataclass
class DegradationLayerStats:
    n_rna: int = 0
    n_mirna: int = 0
    n_protein: int = 0
    n_histone: int = 0
    n_promoter_active: int = 0
    n_metabolite: int = 0
    n_extra: int = 0
    n_protein_skipped_ups: int = 0
    rna_half_life_blocks: float = 5.0
    mirna_half_life_blocks: float | None = None
    protein_half_life_blocks: float = 20.0
    histone_half_life_blocks: float | None = None
    promoter_active_half_life_blocks: float | None = None
    metabolite_half_life_blocks: float | None = None
    extra_places: dict[str, float] = field(default_factory=dict)
    sink_places: list[str] = field(default_factory=list)
    skipped_protein_places: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_rna_degradation": self.n_rna,
            "n_mirna_degradation": self.n_mirna,
            "n_protein_degradation": self.n_protein,
            "n_protein_skipped_ups": self.n_protein_skipped_ups,
            "n_histone_degradation": self.n_histone,
            "n_promoter_active_degradation": self.n_promoter_active,
            "n_metabolite_degradation": self.n_metabolite,
            "n_extra_degradation": self.n_extra,
            "rna_half_life_blocks": self.rna_half_life_blocks,
            "mirna_half_life_blocks": self.mirna_half_life_blocks,
            "protein_half_life_blocks": self.protein_half_life_blocks,
            "histone_half_life_blocks": self.histone_half_life_blocks,
            "promoter_active_half_life_blocks": self.promoter_active_half_life_blocks,
            "metabolite_half_life_blocks": self.metabolite_half_life_blocks,
            "extra_places": dict(self.extra_places),
            "sink_places": self.sink_places,
            "skipped_protein_places": list(self.skipped_protein_places),
        }


def _inhibition_weight_for_half_life(half_life_blocks: float) -> float:
    """
    Map discrete half-life (blocks) to Ruths self-inhibition weight.

    One degradation firing per block; simulator removes
    Uniform[1, floor(marking * weight)]. Approximate E[loss] ≈ marking * weight/2,
    so (1 - weight/2)^hl ≈ 0.5  →  weight ≈ 2*(1 - 2^(-1/hl)).
    """
    hl = max(float(half_life_blocks), 1.0)
    return 2.0 * (1.0 - 2.0 ** (-1.0 / hl))


def first_order_turnover_weight(half_life_blocks: float) -> float:
    """Same discrete map as degradation weight; used for reaction out-route budgets (not sinks)."""
    return _inhibition_weight_for_half_life(half_life_blocks)


def _is_promoter_active_place(place: str) -> bool:
    return place.endswith(".promoter.active") or ".promoter.active" in place


def add_molecule_degradation_layer(
    net: SpnPetriNet,
    *,
    rna_half_life_blocks: float = 5.0,
    mirna_half_life_blocks: float | None = None,
    protein_half_life_blocks: float = 20.0,
    include_protein_complex: bool = False,
    extra_places: Mapping[str, float] | None = None,
    histone_half_life_blocks: float | None = None,
    promoter_active_half_life_blocks: float | None = None,
    metabolite_half_life_blocks: float | None = None,
    skip_protein_places: AbstractSet[str] | None = None,
) -> DegradationLayerStats:
    """
    Add degradation transitions (destroy-only; no sink places).

    Topology (per molecule place P)::

        P -->> {P};degradation

    Simulation maps each degradation transition to a **self-inhibition** on P
    (token outlet). Destroyed mass is ledger-only (``sink_remove``); no
    ``__SPN_sink__.*`` bookkeeping places are created.

    miRNA places use ``mirna_half_life_blocks`` when set, else ``rna_half_life_blocks``.
    This outlet matters because translational repression does **not** consume miRNA.

    ``extra_places`` maps place name → half-life (blocks). Defaults to
    ``DEFAULT_EXTRA_DEGRADE_PLACES`` when ``extra_places`` is ``None``.
    Pass an empty dict to disable extras.

    When ``metabolite_half_life_blocks`` is set, every METABOLITE place gets a
    default half-life (extras still override). Same pattern for histone /
    promoter.active when their half-lives are set.
    """
    mir_hl = (
        float(rna_half_life_blocks)
        if mirna_half_life_blocks is None
        else float(mirna_half_life_blocks)
    )
    extras = (
        dict(DEFAULT_EXTRA_DEGRADE_PLACES) if extra_places is None else dict(extra_places)
    )
    skip_prot = set(skip_protein_places or ())
    stats = DegradationLayerStats(
        rna_half_life_blocks=rna_half_life_blocks,
        mirna_half_life_blocks=mir_hl,
        protein_half_life_blocks=protein_half_life_blocks,
        histone_half_life_blocks=histone_half_life_blocks,
        promoter_active_half_life_blocks=promoter_active_half_life_blocks,
        metabolite_half_life_blocks=metabolite_half_life_blocks,
        extra_places=extras,
        sink_places=[],
        skipped_protein_places=sorted(skip_prot),
    )

    if not hasattr(net, "degradation_weights"):
        net.degradation_weights = {}  # type: ignore[attr-defined]

    attached: set[str] = set()

    def _attach(place: str, hl: float) -> None:
        if place in attached:
            return
        tid = f"{place};degradation"
        net.transitions[tid] = SpnTransitionKind.DEGRADATION
        w = _inhibition_weight_for_half_life(hl)
        net.degradation_weights[place] = w  # type: ignore[attr-defined]
        if (place, tid) not in net.arcs:
            net.arcs.append((place, tid))
        attached.add(place)

    for place, kind in list(net.places.items()):
        if place.startswith("__SPN_sink__"):
            continue
        if place in extras:
            continue
        if kind == SpnPlaceKind.RNA:
            if is_mirna_place(place):
                _attach(place, mir_hl)
                stats.n_mirna += 1
            else:
                _attach(place, rna_half_life_blocks)
                stats.n_rna += 1
        elif kind == SpnPlaceKind.PROTEIN or (
            kind == SpnPlaceKind.PROTEIN_COMPLEX and include_protein_complex
        ):
            if place in skip_prot:
                stats.n_protein_skipped_ups += 1
                continue
            _attach(place, protein_half_life_blocks)
            stats.n_protein += 1
        elif kind == SpnPlaceKind.HISTONE_MARK and histone_half_life_blocks is not None:
            _attach(place, float(histone_half_life_blocks))
            stats.n_histone += 1
        elif _is_promoter_active_place(place) and promoter_active_half_life_blocks is not None:
            _attach(place, float(promoter_active_half_life_blocks))
            stats.n_promoter_active += 1
        elif kind == SpnPlaceKind.METABOLITE and metabolite_half_life_blocks is not None:
            _attach(place, float(metabolite_half_life_blocks))
            stats.n_metabolite += 1

    for place, hl in extras.items():
        if place not in net.places or place.startswith("__SPN_sink__"):
            continue
        _attach(place, float(hl))
        stats.n_extra += 1

    return stats


def place_kind_for_degradation(name: str) -> SpnPlaceKind:
    return classify_place(name)

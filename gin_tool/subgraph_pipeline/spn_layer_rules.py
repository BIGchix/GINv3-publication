"""Reusable SPN layer-rule semantics (v2): signal/regulatory vs metabolic transfer.

Rules are kind/name-driven for any subgraph — not pathway-specific patches.
"""

from __future__ import annotations

import math
from enum import Enum
from typing import Any, Mapping, Sequence

from gin_tool.spn_export import (
    SpnPlaceKind,
    _is_histone_mark,
    _is_metabolite,
    _is_protein_complex,
    _is_rna,
    classify_place,
)


class FireMode(str, Enum):
    """Per-interaction firing semantics under layer_rules=v2."""

    LEGACY_CONSUME = "legacy_consume"
    """Historical gin_consume_v1 (activation transfers; inhibition drains both)."""

    TRANSFER_CONSERVING = "transfer_conserving"
    """Activation: move δ from source to target with δ ≤ source (no minting)."""

    REGULATORY = "regulatory"
    """Bounded occupancy on chromatin/activity places; source not consumed.

    Target must be promoter.active / ;state / histone mark. Updates use logistic
    saturation in ``[0, M]`` — never unbounded ``target += δ`` minting.
    """

    CATALYST = "catalyst"
    """Protein enzyme/kinase: activation = presence gate; inhibition = drain target only.

    Metabolite mass moves on met↔met TRANSFER edges; protein copy number via
    transcription / translation / degradation / complex transfer — never
    ``target += f(enzyme)`` on activation.
    """
    PRODUCE_TRANSCRIPTION = "produce_transcription"
    """ΔRNA = α·f(m_active); do not consume promoter.active."""

    PRODUCE_TRANSLATION = "produce_translation"
    """Δprotein = α·min(m_RNA, Mmax); do not consume RNA."""

    DEGRADE = "degrade"
    """Self-drain; destroyed tokens should be credited to a sink place."""


LAYER_RULES_V2 = "v2"
LAYER_RULES_OFF = "off"

# Default capacity for chromatin activity places (token units ≈ high_tokens scale).
DEFAULT_ACTIVITY_SIGMOID_M_MAX = 20.0
DEFAULT_ACTIVITY_SIGMOID_K = 1.0


def _is_promoter_active(name: str) -> bool:
    return name.endswith(".promoter.active") or ".promoter.active" in name


def is_sigmoid_activity_place(name: str) -> bool:
    """Places whose level is a saturating activity/occupancy signal.

    Includes promoter.active / promoter;state and histone marks (H3K*/H4K*/…).
    """
    if _is_histone_mark(name):
        return True
    low = name.lower()
    return (
        ".promoter.active" in low
        or low.endswith("promoter.active")
        or ".promoter;state" in low
        or low.endswith("promoter;state")
    )


def _is_protein_complex_place(name: str) -> bool:
    """GIN complexes are typically ``A_B`` (underscore, no ``.`` / ``;``)."""
    if _is_protein_complex(name):
        return True
    try:
        return classify_place(name) == SpnPlaceKind.PROTEIN_COMPLEX
    except ValueError:
        return False


def _is_protein_monomer(name: str) -> bool:
    """Free protein / monomer place (not A_B complex)."""
    if _is_protein_complex_place(name):
        return False
    if name.endswith(".protein"):
        return True
    try:
        return classify_place(name) == SpnPlaceKind.PROTEIN
    except ValueError:
        return False


def _is_protein_place(name: str) -> bool:
    """Any protein-like place (monomer or complex)."""
    if _is_protein_complex_place(name) or name.endswith(".protein"):
        return True
    try:
        return classify_place(name) in (SpnPlaceKind.PROTEIN, SpnPlaceKind.PROTEIN_COMPLEX)
    except ValueError:
        return False


def resolve_fire_mode(
    source: str,
    target: str,
    kind: str,
    *,
    transition_id: str | None = None,
    layer_rules: str = LAYER_RULES_OFF,
) -> FireMode:
    """Map an SPN interaction to a fire mode."""
    if layer_rules != LAYER_RULES_V2:
        return FireMode.LEGACY_CONSUME

    tid = (transition_id or "").lower()
    if kind == "inhibition" and source == target:
        return FireMode.DEGRADE
    if "degradation" in tid:
        return FireMode.DEGRADE
    if "transcription" in tid:
        return FireMode.PRODUCE_TRANSCRIPTION
    if "translation" in tid:
        return FireMode.PRODUCE_TRANSLATION

    # Histone marks or TF/protein → promoter / chromatin: regulatory
    if _is_histone_mark(source):
        return FireMode.REGULATORY
    if _is_protein_place(source) and (
        _is_promoter_active(target) or ".promoter" in target.lower()
    ):
        return FireMode.REGULATORY

    # Complex assembly / dissociation (A_B): consume/produce tokens.
    if _is_protein_complex_place(source) or _is_protein_complex_place(target):
        return FireMode.TRANSFER_CONSERVING

    # Enzyme–metabolite: protein catalyst is not consumed (e.g. PTEN: PIP3→PIP2).
    if _is_protein_monomer(source) and _is_metabolite(target):
        return FireMode.CATALYST
    # Kinase / enzyme acting on protein substrate: catalyst not consumed.
    if _is_protein_monomer(source) and _is_protein_monomer(target) and source != target:
        return FireMode.CATALYST

    # Metabolite↔metabolite (and other met-involved non-catalyst): conserving transfer
    if _is_metabolite(source) and _is_metabolite(target):
        return FireMode.TRANSFER_CONSERVING
    if _is_metabolite(source) or _is_metabolite(target):
        return FireMode.TRANSFER_CONSERVING

    # Default internal edges: conserving transfer (safer than minting)
    return FireMode.TRANSFER_CONSERVING


SINK_RNA = "__SPN_sink__.RNA"
SINK_PROTEIN = "__SPN_sink__.protein"
SINK_METABOLITE = "__SPN_sink__.metabolite"
SINK_HISTONE = "__SPN_sink__.histone"
SINK_CHROMATIN = "__SPN_sink__.chromatin"


def logistic_sigmoid(z: float) -> float:
    """Numerically stable σ(z) = 1/(1+e^{-z})."""
    if z >= 0.0:
        ez = math.exp(-z)
        return 1.0 / (1.0 + ez)
    ez = math.exp(z)
    return ez / (1.0 + ez)


def logit01(p: float, *, eps: float = 1e-9) -> float:
    p = min(1.0 - eps, max(eps, float(p)))
    return math.log(p / (1.0 - p))


def regulatory_sigmoid_update(
    level: float,
    delta: float,
    *,
    activate: bool,
    m_max: float = DEFAULT_ACTIVITY_SIGMOID_M_MAX,
    steepness: float = DEFAULT_ACTIVITY_SIGMOID_K,
    min_tokens: float = 0.0,
) -> float:
    """Update activity in ``[min_tokens, m_max]`` via logit-space ± drive.

    Biology: promoter/chromatin activity saturates; successive TF/histone hits
    yield diminishing returns (logistic). ``delta`` is scaled by ``m_max`` so
    typical token drives (~1..high_tokens) stay O(1) on the logit axis.
    """
    m_max = float(m_max)
    if m_max <= 0.0:
        raw = float(level) + (float(delta) if activate else -float(delta))
        return max(min_tokens, raw)
    # Near-empty: logit(ε) ≈ −∞ so O(1) drives cannot leave the floor. Warm-start
    # with a clipped linear step, then subsequent hits use logistic saturation.
    if float(level) <= max(float(min_tokens), 1e-12 * m_max):
        step = float(delta) if activate else -float(delta)
        return max(min_tokens, min(m_max, float(level) + step))
    frac = min(1.0 - 1e-9, max(1e-9, float(level) / m_max))
    z = logit01(frac)
    drive = float(steepness) * (float(delta) / m_max)
    z_new = z + drive if activate else z - drive
    return max(min_tokens, min(m_max, m_max * logistic_sigmoid(z_new)))


def regulatory_relu_update(
    level: float,
    delta: float,
    *,
    activate: bool,
    m_max: float = 0.0,
    min_tokens: float = 0.0,
) -> float:
    """Linear ±δ with ReLU floor: ``m ← max(min_tokens, m ± δ)``.

    If ``m_max > 0``, also clamp to ``m_max`` (optional ceiling). ``m_max <= 0``
    means no upper bound (pure ReLU / unbounded linear).
    """
    step = float(delta) if activate else -float(delta)
    new = float(level) + step
    if float(m_max) > 0.0:
        new = min(float(m_max), new)
    return max(float(min_tokens), new)


ACTIVITY_UPDATE_SIGMOID = "sigmoid"
ACTIVITY_UPDATE_RELU = "relu"
ACTIVITY_UPDATE_RELU_ACTIVE = "relu_active"
"""ReLU only on ``promoter.active`` / ``promoter;state``; histone marks stay sigmoid."""


def is_promoter_activity_signal(name: str) -> bool:
    """Public alias: ``promoter.active`` / ``promoter;state`` (not histone marks)."""
    return _is_promoter_activity_signal(name)


def is_tf_dna_regulatory_edge(
    source: str,
    target: str,
    tf_genes: set[str] | frozenset[str],
) -> bool:
    """TF protein → promoter hub (DNA-facing); exempt from REGULATORY out-route fan-out."""
    if not _is_protein_place(source):
        return False
    gene = source[: -len(".protein")] if source.endswith(".protein") else source
    if gene not in tf_genes:
        return False
    return _is_promoter_activity_signal(target)


def _is_promoter_activity_signal(name: str) -> bool:
    """``promoter.active`` / ``promoter;state`` (not histone mark places)."""
    low = str(name).lower()
    return (
        ".promoter.active" in low
        or low.endswith("promoter.active")
        or ".promoter;state" in low
        or low.endswith("promoter;state")
    )


def resolve_place_activity_update(
    place: str,
    mode: str,
    *,
    m_max: float,
) -> tuple[str, float]:
    """Map global ``activity_update`` + place → ``(effective_mode, effective_m_max)``.

    - ``sigmoid``: logit update; ``m_max<=0`` forced to default capacity by caller if needed
    - ``relu``: ReLU on all activity places (hist + active)
    - ``relu_active``: ReLU on promoter.active/;state with **no ceiling** (``m_max=0``);
      histone marks keep sigmoid with the configured ``m_max``
    """
    m = (mode or ACTIVITY_UPDATE_SIGMOID).strip().lower()
    if m == ACTIVITY_UPDATE_RELU:
        return ACTIVITY_UPDATE_RELU, float(m_max)
    if m == ACTIVITY_UPDATE_RELU_ACTIVE:
        if _is_promoter_activity_signal(place):
            return ACTIVITY_UPDATE_RELU, 0.0
        return ACTIVITY_UPDATE_SIGMOID, float(m_max)
    return ACTIVITY_UPDATE_SIGMOID, float(m_max)


def scale_promoter_active_relu_drive(
    delta: float,
    *,
    place: str,
    mode: str,
    scale: float = 1.0,
    scale_mode: str = "drive",
) -> float:
    """Apply empirical α to REGULATORY δ when ``scale_mode=='drive'``.

    Histone marks are never scaled. ``scale=1`` or ``scale_mode=='net'`` is a no-op
    here (net scaling is applied end-of-block in the simulator).
    """
    if str(scale_mode or "drive").strip().lower() != "drive":
        return float(delta)
    s = float(scale)
    if s == 1.0 or not _is_promoter_activity_signal(place):
        return float(delta)
    eff_mode, _ = resolve_place_activity_update(place, mode, m_max=1.0)
    if eff_mode != ACTIVITY_UPDATE_RELU:
        return float(delta)
    return float(delta) * s


def apply_promoter_active_relu_net_scale(
    marking: Any,
    idx: Mapping[str, int],
    places: Sequence[str],
    marking_start: Any,
    *,
    mode: str,
    scale: float,
    min_tokens: float = 0.0,
) -> int:
    """End-of-block: ``m ← m0 + α·(m−m0)`` for promoter.active under ReLU write mode.

    Returns number of places adjusted. No-op if ``scale==1``.
    """
    s = float(scale)
    if s == 1.0:
        return 0
    n = 0
    for p in places:
        if p not in idx:
            continue
        if not _is_promoter_activity_signal(p):
            continue
        eff_mode, _ = resolve_place_activity_update(p, mode, m_max=1.0)
        if eff_mode != ACTIVITY_UPDATE_RELU:
            continue
        i = idx[p]
        m0 = float(marking_start[i])
        m1 = float(marking[i])
        marking[i] = max(float(min_tokens), m0 + s * (m1 - m0))
        n += 1
    return n


def regulatory_activity_update(
    level: float,
    delta: float,
    *,
    activate: bool,
    mode: str = ACTIVITY_UPDATE_SIGMOID,
    m_max: float = DEFAULT_ACTIVITY_SIGMOID_M_MAX,
    steepness: float = DEFAULT_ACTIVITY_SIGMOID_K,
    min_tokens: float = 0.0,
    place: str | None = None,
) -> float:
    """Dispatch REGULATORY writes: ``sigmoid``, ``relu``, or ``relu_active`` (+ place)."""
    eff_mode, eff_m = (mode or ACTIVITY_UPDATE_SIGMOID), float(m_max)
    if place is not None:
        eff_mode, eff_m = resolve_place_activity_update(place, mode, m_max=m_max)
    if eff_mode == ACTIVITY_UPDATE_RELU:
        return regulatory_relu_update(
            level, delta, activate=activate, m_max=eff_m, min_tokens=min_tokens
        )
    return regulatory_sigmoid_update(
        level,
        delta,
        activate=activate,
        m_max=eff_m,
        steepness=steepness,
        min_tokens=min_tokens,
    )


def activity_readout(
    level: float,
    *,
    m_max: float = DEFAULT_ACTIVITY_SIGMOID_M_MAX,
    steepness: float = DEFAULT_ACTIVITY_SIGMOID_K,
    mode: str = ACTIVITY_UPDATE_SIGMOID,
    place: str | None = None,
) -> float:
    """Readout for transcription flux.

    Sigmoid: ``m_max · σ(k·(2m/M−1))``. ReLU / relu_active on active: identity ``max(0, m)``.
    """
    eff_mode, _eff_m = (mode or ACTIVITY_UPDATE_SIGMOID), float(m_max)
    if place is not None:
        eff_mode, _eff_m = resolve_place_activity_update(place, mode, m_max=m_max)
    if eff_mode == ACTIVITY_UPDATE_RELU:
        return max(0.0, float(level))
    m_max = float(m_max)
    if m_max <= 0.0:
        return max(0.0, float(level))
    z = float(steepness) * (2.0 * float(level) / m_max - 1.0)
    return m_max * logistic_sigmoid(z)


def production_amount(
    level: float,
    *,
    alpha: float,
    m_max: float,
    rng: Any,
    distribution: str = "uniform",
) -> float:
    """Δ = α · sample(1..hi); hi = floor(level) or floor(min(level, Mmax)) when Mmax>0."""
    import numpy as np

    capped = min(float(level), float(m_max)) if m_max > 0 else float(level)
    if capped <= 0 or alpha <= 0 or not np.isfinite(capped):
        return 0.0
    # Keep hi in exact float64 integer range (not int64.max via float).
    max_hi = (1 << 53) - 1
    hi = int(np.floor(capped))
    if hi > max_hi:
        hi = max_hi
    if hi < 1:
        return 0.0
    if distribution == "uniform":
        base = float(rng.integers(1, hi + 1))
    else:
        base = float(min(hi, max(1, int(rng.exponential(hi / 3.0)))))
    return float(alpha) * base


def effective_protein_production_output_scale(
    *,
    protein_production_output_scale: float = 1.0,
    translation_output_scale: float = 1.0,
) -> float:
    """Global multiplier on activating ``PRODUCE_TRANSLATION`` mint (layer_rules v2).

    Applies to every RNA→protein translation edge; miRNA/inhibitory repress unchanged.
    ``translation_output_scale`` is retained as a legacy alias (product of both).
    """
    return float(protein_production_output_scale) * float(translation_output_scale)


def sink_name_for_place(place: str) -> str:
    """Choose degradation sink place for a real place name."""
    if _is_rna(place):
        return SINK_RNA
    if _is_histone_mark(place):
        return SINK_HISTONE
    if _is_promoter_active(place) or ".promoter" in place.lower():
        return SINK_CHROMATIN
    if _is_metabolite(place):
        return SINK_METABOLITE
    return SINK_PROTEIN


def empty_ledger() -> dict[str, float]:
    return {
        "seed_inject": 0.0,
        "prod_transcription": 0.0,
        "prod_translation": 0.0,
        # sink_remove kept as alias of tokens_destroyed (legacy ledger key)
        "sink_remove": 0.0,
        "tokens_destroyed": 0.0,
        "transfer_moved": 0.0,
        "regulatory_delta": 0.0,
        "catalyst_gate": 0.0,
        "regulatory_skipped": 0.0,
        "out_route_events": 0.0,
        "out_route_flux": 0.0,
        "illegal_create": 0.0,
        "tf_proteolysis": 0.0,
        "translation_repress": 0.0,
        "hme_regulatory_scaled": 0.0,
        "chromhmm_emw_applied": 0.0,
    }

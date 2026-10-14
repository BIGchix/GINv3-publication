"""GIN-aligned SPN simulator: token-consumption and layer-rules v2 firing."""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import AbstractSet, Any, Literal, Mapping, Sequence

import numpy as np

from gin_tool.subgraph_pipeline.spn_layer_rules import (
    ACTIVITY_UPDATE_RELU,
    ACTIVITY_UPDATE_SIGMOID,
    DEFAULT_ACTIVITY_SIGMOID_K,
    DEFAULT_ACTIVITY_SIGMOID_M_MAX,
    LAYER_RULES_OFF,
    LAYER_RULES_V2,
    FireMode,
    activity_readout,
    apply_promoter_active_relu_net_scale,
    empty_ledger,
    is_promoter_activity_signal,
    is_sigmoid_activity_place,
    is_tf_dna_regulatory_edge,
    production_amount,
    regulatory_activity_update,
    resolve_fire_mode,
    resolve_place_activity_update,
    scale_promoter_active_relu_drive,
    effective_protein_production_output_scale,
)
from gin_tool.subgraph_pipeline.interaction_weight_schedule import (
    InteractionWeightPulse,
    scale_interaction_tuples_for_block,
)
from gin_tool.subgraph_pipeline.spn_degradation import first_order_turnover_weight
from gin_tool.spn_export import is_mirna_place

Distribution = Literal["uniform", "normal"]


def default_hme_protein_places(*, node_list_dir: Path | None = None) -> frozenset[str]:
    """``GENE.protein`` places for HME enzymes from ``nodeList/nodeList_HME.txt``."""
    from gin_tool.hme_chain_repair import load_hme_symbols

    root = Path(__file__).resolve().parents[2]
    nld = Path(node_list_dir) if node_list_dir is not None else root / "nodeList"
    try:
        syms = load_hme_symbols(nld)
    except OSError:
        return frozenset()
    return frozenset(f"{s}.protein" for s in syms if s)
GIN_FIRING_POLICY = "gin_consume_v1"
FIRE_PROPENSITY_GEOM_MEAN = "geom_mean"
# float64 cannot uniquely represent int64 max; keep hi in exact float integer range.
_MAX_SAMPLE_HI = (1 << 53) - 1  # 9007199254740991
GIN_FIRING_POLICY_V2 = "layer_rules_v2"


@dataclass
class SimulationResult:
    place_names: list[str]
    mean_markings: np.ndarray
    std_markings: np.ndarray
    n_runs: int
    n_blocks: int
    ledger_totals: dict[str, float] = field(default_factory=dict)
    firing_policy: str = GIN_FIRING_POLICY
    nf_passrate: dict | None = None
    """Optional per-run NF phase summary when nf_qc_pair was set."""
    nf_per_run: list[dict] | None = None
    """Per-run phase dicts (same order as MC runs)."""
    pair_block_samples: np.ndarray | None = None
    """Optional (n_runs, n_blocks+1, 2) for nf_qc_pair places only."""
    pair_place_names: tuple[str, str] | None = None

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(
                {
                    "places": self.place_names,
                    "n_runs": self.n_runs,
                    "n_blocks": self.n_blocks,
                    "firing_policy": self.firing_policy,
                    "ledger_totals": self.ledger_totals,
                    "mean_markings": self.mean_markings.tolist(),
                    "std_markings": self.std_markings.tolist(),
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    def save_csv(self, path: str | Path) -> None:
        out = Path(path)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["block"] + self.place_names)
            for b in range(self.mean_markings.shape[0]):
                writer.writerow([b] + [f"{v:.4f}" for v in self.mean_markings[b]])


def _sample_amount(
    rng: np.random.Generator,
    upper: float,
    *,
    distribution: Distribution,
) -> float:
    if upper <= 0 or not np.isfinite(upper):
        return 0.0
    # Guard against token blow-ups overflowing rng.integers int64 high.
    # Do not min() via float(_MAX_SAMPLE_HI): float(int64.max) rounds above max.
    hi = int(np.floor(upper))
    if hi > _MAX_SAMPLE_HI:
        hi = _MAX_SAMPLE_HI
    if hi < 1:
        return 0.0
    if distribution == "normal":
        val = int(round(rng.normal(loc=(1 + hi) / 2, scale=max(1.0, hi / 3))))
        val = int(np.clip(val, 1, hi))
        return float(val)
    return float(rng.integers(1, hi + 1))


def _credit_sink(
    marking: np.ndarray,
    idx: dict[str, int],
    place: str,
    amount: float,
    ledger: dict[str, float] | None,
) -> None:
    """Record destroyed mass in the ledger only (no sink places).

    ``marking`` / ``idx`` / ``place`` are kept for call-site compatibility;
    tokens are destroyed at the source place before this helper runs.
    Legacy ``__SPN_sink__.*`` places are ignored even if present in a model.
    """
    if amount <= 0:
        return
    if ledger is not None:
        ledger["sink_remove"] = float(ledger.get("sink_remove", 0.0)) + amount
        ledger["tokens_destroyed"] = float(ledger.get("tokens_destroyed", 0.0)) + amount


def _record_target_flux(
    flux_events: list[dict[str, Any]] | None,
    flux_watch_targets: frozenset[str] | set[str] | None,
    *,
    block: int,
    source: str,
    target: str,
    kind: str,
    mode: FireMode,
    before: float,
    after: float,
) -> None:
    if flux_events is None or not flux_watch_targets or target not in flux_watch_targets:
        return
    delta = float(after) - float(before)
    if abs(delta) < 1e-15:
        return
    flux_events.append(
        {
            "block": int(block),
            "source": str(source),
            "target": str(target),
            "kind": str(kind),
            "mode": mode.value if isinstance(mode, FireMode) else str(mode),
            "delta_target": delta,
        }
    )


def _fire_one(
    marking: np.ndarray,
    idx: dict[str, int],
    source: str,
    target: str,
    kind: str,
    weight: float,
    rng: np.random.Generator,
    *,
    distribution: Distribution,
    min_tokens: float,
    fire_mode: FireMode = FireMode.LEGACY_CONSUME,
    ledger: dict[str, float] | None = None,
    prod_alpha_transcription: float = 1.0,
    prod_alpha_translation: float = 1.0,
    prod_m_max: float = 50.0,
    translation_output_scale: float = 1.0,
    protein_production_output_scale: float = 1.0,
    tf_proteolysis_eps: float = 0.0,
    activity_sigmoid_m_max: float = DEFAULT_ACTIVITY_SIGMOID_M_MAX,
    activity_sigmoid_k: float = DEFAULT_ACTIVITY_SIGMOID_K,
    activity_update_mode: str = ACTIVITY_UPDATE_SIGMOID,
    activity_relu_drive_scale: float = 1.0,
    activity_relu_scale_mode: str = "drive",
    force_delta: float | None = None,
    hme_regulatory_budget_scale: float = 1.0,
    hme_protein_places: frozenset[str] | set[str] | None = None,
    chromhmm_mark_weights: Mapping[str, float] | None = None,
    flux_events: list[dict[str, Any]] | None = None,
    flux_watch_targets: frozenset[str] | set[str] | None = None,
    current_block: int = -1,
) -> None:
    if source not in idx or target not in idx:
        return
    u = idx[source]
    v = idx[target]
    w = max(float(weight), 1e-8)
    mode = fire_mode if isinstance(fire_mode, FireMode) else FireMode(str(fire_mode))

    if mode == FireMode.LEGACY_CONSUME:
        src_tokens = marking[u] * w
        delta = _sample_amount(rng, src_tokens, distribution=distribution)
        if delta <= 0:
            return
        if kind == "activation":
            marking[v] += delta
            marking[u] = max(min_tokens, marking[u] - delta)
        else:
            marking[v] = max(min_tokens, marking[v] - delta)
            if u != v:
                marking[u] = max(min_tokens, marking[u] - delta)
        return

    if mode == FireMode.DEGRADE:
        src_tokens = marking[u] * w
        delta = _sample_amount(rng, src_tokens, distribution=distribution)
        if delta <= 0:
            return
        before = float(marking[u])
        marking[u] = max(min_tokens, marking[u] - delta)
        if is_sigmoid_activity_place(source):
            _mode, eff_m = resolve_place_activity_update(
                source, activity_update_mode, m_max=float(activity_sigmoid_m_max)
            )
            if eff_m > 0:
                marking[u] = min(float(eff_m), float(marking[u]))
        destroyed = before - float(marking[u])
        _credit_sink(marking, idx, source, destroyed, ledger)
        return

    if mode == FireMode.PRODUCE_TRANSCRIPTION:
        raw = float(marking[u])
        if is_sigmoid_activity_place(source):
            eff = activity_readout(
                raw,
                m_max=activity_sigmoid_m_max,
                steepness=activity_sigmoid_k,
                mode=activity_update_mode,
                place=source,
            )
        else:
            eff = raw
        amount = production_amount(
            eff,
            alpha=prod_alpha_transcription,
            m_max=prod_m_max,
            rng=rng,
            distribution=distribution,
        )
        if amount <= 0:
            return
        marking[v] += amount
        if ledger is not None:
            ledger["prod_transcription"] = float(ledger.get("prod_transcription", 0.0)) + amount
        return

    if mode == FireMode.PRODUCE_TRANSLATION:
        # Inhibitory / miRNA inputs repress translation (never mint from repressor level).
        if kind == "inhibition" or is_mirna_place(source):
            amount = production_amount(
                float(marking[u]),
                alpha=prod_alpha_translation,
                m_max=prod_m_max,
                rng=rng,
                distribution=distribution,
            )
            if amount <= 0:
                return
            available = max(0.0, float(marking[v]) - min_tokens)
            removed = min(available, amount)
            if removed <= 0:
                return
            marking[v] = max(min_tokens, marking[v] - removed)
            _credit_sink(marking, idx, target, removed, ledger)
            if ledger is not None:
                ledger["translation_repress"] = float(
                    ledger.get("translation_repress", 0.0)
                ) + removed
            return
        amount = production_amount(
            float(marking[u]),
            alpha=prod_alpha_translation,
            m_max=prod_m_max,
            rng=rng,
            distribution=distribution,
        )
        if amount <= 0:
            return
        tl_scale = effective_protein_production_output_scale(
            protein_production_output_scale=protein_production_output_scale,
            translation_output_scale=translation_output_scale,
        )
        if tl_scale != 1.0:
            amount = float(amount) * tl_scale
        if amount <= 0:
            return
        marking[v] += amount
        if ledger is not None:
            ledger["prod_translation"] = float(ledger.get("prod_translation", 0.0)) + amount
        return

    if mode == FireMode.CATALYST:
        # Enzyme / kinase presence: do not consume the catalyst place.
        # Inhibition must still destroy substrate tokens (MDM2⊣TP53 etc.).
        # Activation remains a pure gate (no mint) — chemistry stays on TRANSFER.
        if kind == "inhibition" or is_mirna_place(source):
            src_tokens = marking[u] * w
            delta = _sample_amount(rng, src_tokens, distribution=distribution)
            if delta <= 0:
                return
            available = max(0.0, float(marking[v]) - min_tokens)
            removed = min(available, float(delta))
            if removed <= 0:
                return
            marking[v] = max(min_tokens, marking[v] - removed)
            _credit_sink(marking, idx, target, removed, ledger)
            if ledger is not None:
                ledger["catalyst_inhibition"] = float(
                    ledger.get("catalyst_inhibition", 0.0)
                ) + removed
            return
        if ledger is not None and float(marking[u]) > 0.0:
            ledger["catalyst_gate"] = float(ledger.get("catalyst_gate", 0.0)) + 1.0
        return

    if mode == FireMode.REGULATORY:
        # Chromatin/activity occupancy only — never unbounded token minting.
        if not is_sigmoid_activity_place(target):
            # Non-activity targets: still honour inhibition by draining tokens.
            if kind == "inhibition" or is_mirna_place(source):
                src_tokens = marking[u] * w
                delta = _sample_amount(rng, src_tokens, distribution=distribution)
                if delta <= 0:
                    return
                available = max(0.0, float(marking[v]) - min_tokens)
                removed = min(available, float(delta))
                if removed <= 0:
                    return
                marking[v] = max(min_tokens, marking[v] - removed)
                _credit_sink(marking, idx, target, removed, ledger)
                if ledger is not None:
                    ledger["regulatory_inhibition"] = float(
                        ledger.get("regulatory_inhibition", 0.0)
                    ) + removed
                return
            if ledger is not None:
                ledger["regulatory_skipped"] = float(ledger.get("regulatory_skipped", 0.0)) + 1.0
            return
        upd, m_max = resolve_place_activity_update(
            target, activity_update_mode, m_max=float(activity_sigmoid_m_max)
        )
        # Sigmoid needs a positive capacity; ReLU allows m_max<=0 = no ceiling.
        if upd != ACTIVITY_UPDATE_RELU and m_max <= 0.0:
            m_max = float(DEFAULT_ACTIVITY_SIGMOID_M_MAX)
        if force_delta is not None:
            # Out-route supplies a raw drive already sampled as f(source level).
            delta = float(force_delta)
        else:
            src_tokens = marking[u] * w
            delta = _sample_amount(rng, src_tokens, distribution=distribution)
        if delta <= 0:
            return
        # HME writer/eraser: scale extracted enzyme drive; source tokens unchanged.
        hme_set = hme_protein_places or ()
        scale = float(hme_regulatory_budget_scale)
        if scale != 1.0 and source in hme_set:
            delta = float(delta) * scale
            if ledger is not None:
                ledger["hme_regulatory_scaled"] = float(
                    ledger.get("hme_regulatory_scaled", 0.0)
                ) + 1.0
        if delta <= 0:
            return
        # SPN v8: chromHMM signed mark weight on hist → promoter.active only.
        from gin_tool.subgraph_pipeline.chromhmm_emission_weights import (
            apply_signed_chromhmm_drive,
        )

        delta, kind, emw_applied = apply_signed_chromhmm_drive(
            delta,
            source=source,
            target=target,
            kind=kind,
            weights=chromhmm_mark_weights,
        )
        if emw_applied and ledger is not None:
            ledger["chromhmm_emw_applied"] = float(
                ledger.get("chromhmm_emw_applied", 0.0)
            ) + 1.0
        if delta <= 0:
            return
        before = float(marking[v])
        delta = scale_promoter_active_relu_drive(
            delta,
            place=target,
            mode=activity_update_mode,
            scale=activity_relu_drive_scale,
            scale_mode=activity_relu_scale_mode,
        )
        if delta <= 0:
            return
        marking[v] = regulatory_activity_update(
            before,
            delta,
            activate=(kind == "activation"),
            mode=upd,
            m_max=m_max,
            steepness=activity_sigmoid_k,
            min_tokens=min_tokens,
            place=target,
        )
        if ledger is not None:
            ledger["regulatory_delta"] = float(ledger.get("regulatory_delta", 0.0)) + abs(
                float(marking[v]) - before
            )
        _record_target_flux(
            flux_events,
            flux_watch_targets,
            block=current_block,
            source=source,
            target=target,
            kind=kind,
            mode=mode,
            before=before,
            after=float(marking[v]),
        )
        # Optional TF→promoter UPS timer (opt-in; default recycle).
        if tf_proteolysis_eps > 0 and u != v:
            removed = min(float(marking[u]) - min_tokens, float(tf_proteolysis_eps) * delta)
            if removed > 0:
                marking[u] = max(min_tokens, marking[u] - removed)
                _credit_sink(marking, idx, source, removed, ledger)
                if ledger is not None:
                    ledger["tf_proteolysis"] = float(ledger.get("tf_proteolysis", 0.0)) + removed
        return

    # TRANSFER_CONSERVING (default metabolic / internal under v2)
    available = max(0.0, float(marking[u]) - min_tokens)
    if available <= 0:
        return
    if force_delta is not None:
        delta = min(float(force_delta), available)
    else:
        # Clamp sample upper by available tokens so w>1 cannot mint.
        if w >= 1.0:
            raw_upper = available
            if w > 1.0 + 1e-12 and ledger is not None:
                # Would-have-minted mass under legacy w>1 sampling (diagnostic only).
                ledger["illegal_create"] = float(ledger.get("illegal_create", 0.0)) + 0.0
        else:
            raw_upper = available * w
        delta = _sample_amount(rng, raw_upper, distribution=distribution)
        if delta <= 0:
            return
        delta = min(delta, available)
    if delta <= 0:
        return
    if kind == "activation":
        if is_sigmoid_activity_place(target):
            before = float(marking[v])
            marking[u] = max(min_tokens, marking[u] - delta)
            d_use = scale_promoter_active_relu_drive(
                delta,
                place=target,
                mode=activity_update_mode,
                scale=activity_relu_drive_scale,
                scale_mode=activity_relu_scale_mode,
            )
            marking[v] = regulatory_activity_update(
                before,
                d_use,
                activate=True,
                mode=activity_update_mode,
                m_max=activity_sigmoid_m_max,
                steepness=activity_sigmoid_k,
                min_tokens=min_tokens,
                place=target,
            )
            _record_target_flux(
                flux_events,
                flux_watch_targets,
                block=current_block,
                source=source,
                target=target,
                kind=kind,
                mode=mode,
                before=before,
                after=float(marking[v]),
            )
        else:
            before = float(marking[v])
            marking[v] += delta
            marking[u] = max(min_tokens, marking[u] - delta)
            _record_target_flux(
                flux_events,
                flux_watch_targets,
                block=current_block,
                source=source,
                target=target,
                kind=kind,
                mode=mode,
                before=before,
                after=float(marking[v]),
            )
        if ledger is not None:
            ledger["transfer_moved"] = float(ledger.get("transfer_moved", 0.0)) + delta
    else:
        # Inhibition: destroy on target (and source if distinct); credit sinks.
        before_v = float(marking[v])
        if is_sigmoid_activity_place(target):
            d_use = scale_promoter_active_relu_drive(
                delta,
                place=target,
                mode=activity_update_mode,
                scale=activity_relu_drive_scale,
                scale_mode=activity_relu_scale_mode,
            )
            marking[v] = regulatory_activity_update(
                before_v,
                d_use,
                activate=False,
                mode=activity_update_mode,
                m_max=activity_sigmoid_m_max,
                steepness=activity_sigmoid_k,
                min_tokens=min_tokens,
                place=target,
            )
            destroyed_v = 0.0  # signal down, not mass to sink
        else:
            marking[v] = max(min_tokens, marking[v] - delta)
            destroyed_v = before_v - float(marking[v])
            _credit_sink(marking, idx, target, destroyed_v, ledger)
        if u != v:
            before_u = float(marking[u])
            marking[u] = max(min_tokens, marking[u] - delta)
            destroyed_u = before_u - float(marking[u])
            _credit_sink(marking, idx, source, destroyed_u, ledger)


InteractionTuple = tuple[str, str, str, float, str, str]  # src,tgt,kind,w,mode,tid

# Modes whose out-edges compete: weights → categorical pick; flux f(m) shared / source / block.
_ROUTED_MODES = frozenset({FireMode.TRANSFER_CONSERVING.value, FireMode.REGULATORY.value})
DEFAULT_OUT_ROUTE_HALF_LIFE_BLOCKS = 20.0

REG_NORM_OFF = "off"
REG_NORM_TARGET_SOFTMAX = "target_softmax"
REG_NORM_HILL_TARGET_SOFTMAX = "hill_target_softmax"
_REG_NORM_HUB_MODES = frozenset({REG_NORM_TARGET_SOFTMAX, REG_NORM_HILL_TARGET_SOFTMAX})


def _hill_occupancy(level: float, k: float) -> float:
    m = max(0.0, float(level))
    kk = max(1e-9, float(k))
    return m / (kk + m)


def _softmax_probs(scores: Sequence[float]) -> np.ndarray:
    x = np.asarray(scores, dtype=np.float64)
    if x.size == 0:
        return x
    x = np.where(np.isfinite(x), x, -1e9)
    x = x - float(np.max(x))
    e = np.exp(np.clip(x, -40.0, 40.0))
    s = float(e.sum())
    if s <= 0.0 or not np.isfinite(s):
        return np.ones(x.size, dtype=np.float64) / float(x.size)
    return e / s


def _split_hub_regulatory(
    interactions: Sequence[InteractionTuple],
) -> tuple[list[InteractionTuple], list[InteractionTuple]]:
    """Peel REGULATORY→promoter.active/;state for target-side competition."""
    hub: list[InteractionTuple] = []
    rest: list[InteractionTuple] = []
    for ix in interactions:
        mode_s = ix[4]
        tgt = ix[1]
        if mode_s == FireMode.REGULATORY.value and is_promoter_activity_signal(tgt):
            hub.append(ix)
        else:
            rest.append(ix)
    return hub, rest


def _fire_hub_regulatory_target_softmax(
    marking: np.ndarray,
    idx: dict[str, int],
    hub_ixs: Sequence[InteractionTuple],
    rng: np.random.Generator,
    *,
    regulatory_norm: str,
    regulatory_hill_k: float,
    regulatory_target_budget: float,
    distribution: Distribution,
    min_tokens: float,
    ledger: dict[str, float] | None,
    out_route_half_life_blocks: float,
    fire_kwargs: dict[str, Any],
) -> None:
    """Compete REGULATORY in-edges on each promoter hub: shared budget × softmax share."""
    if not hub_ixs:
        return
    groups: dict[tuple[str, str], list[InteractionTuple]] = {}
    for ix in hub_ixs:
        groups.setdefault((ix[1], ix[2]), []).append(ix)

    keys = list(groups.keys())
    rng.shuffle(keys)
    use_hill = regulatory_norm == REG_NORM_HILL_TARGET_SOFTMAX
    hl = float(out_route_half_life_blocks)
    budget_cap = max(0.0, float(regulatory_target_budget))
    eps = 1e-9

    for tgt, kind in keys:
        ixs = groups[(tgt, kind)]
        scores: list[float] = []
        avail: list[float] = []
        for src, _t, _k, weight, _m, _tid in ixs:
            if src not in idx:
                scores.append(-1e9)
                avail.append(0.0)
                continue
            a = max(0.0, float(marking[idx[src]]) - min_tokens)
            avail.append(a)
            w = abs(float(weight))
            if a <= 0.0:
                scores.append(-1e9)
                continue
            if use_hill:
                # 1+2: saturate occupancy then compete (∝ hill(m), not raw m)
                occ = _hill_occupancy(a, regulatory_hill_k)
                scores.append(math.log(occ + eps) + math.log(w + eps))
            else:
                # mode2: equal footing among active sources (budget share ≉ abundance).
                # Edge |w| still modulates; presence is binary so high-token hubs cannot
                # monopolize via raw m (softmax(log m) would collapse to ∝ m).
                scores.append(math.log(w + eps))
        probs = _softmax_probs(scores)
        if budget_cap <= 0.0:
            continue
        if hl > 0.0:
            B = _sample_first_order_route_amount(
                rng, budget_cap, half_life_blocks=hl, distribution=distribution
            )
        else:
            B = _sample_amount(rng, budget_cap, distribution=distribution)
        if B <= 0.0:
            continue

        order = rng.permutation(len(ixs))
        for j in order:
            j = int(j)
            if float(probs[j]) <= 0.0 or avail[j] <= 0.0:
                continue
            delta = float(B) * float(probs[j])
            if delta <= 0.0:
                continue
            src, tgt_i, kind_i, weight, mode_s, _tid = ixs[j]
            _fire_one(
                marking,
                idx,
                src,
                tgt_i,
                kind_i,
                weight,
                rng,
                distribution=distribution,
                min_tokens=min_tokens,
                fire_mode=FireMode(mode_s),
                ledger=ledger,
                force_delta=delta,
                **fire_kwargs,
            )
        if ledger is not None:
            ledger["reg_norm_hub_events"] = float(
                ledger.get("reg_norm_hub_events", 0.0)
            ) + 1.0
            ledger["reg_norm_hub_flux"] = float(ledger.get("reg_norm_hub_flux", 0.0)) + float(
                B
            )


def _partition_routed_interactions(
    interactions: Sequence[InteractionTuple],
) -> tuple[dict[tuple[str, str, str], list[InteractionTuple]], list[InteractionTuple]]:
    """Split into out-route groups ``(src, mode, kind)`` vs independently fired edges."""
    routed: dict[tuple[str, str, str], list[InteractionTuple]] = {}
    independent: list[InteractionTuple] = []
    for ix in interactions:
        mode_s = ix[4]
        if mode_s in _ROUTED_MODES:
            key = (ix[0], mode_s, ix[2])
            routed.setdefault(key, []).append(ix)
        else:
            independent.append(ix)
    return routed, independent


def _split_tf_dna_regulatory(
    interactions: Sequence[InteractionTuple],
    tf_genes: set[str] | frozenset[str],
) -> tuple[list[InteractionTuple], list[InteractionTuple]]:
    """Peel TF.protein→promoter hub REGULATORY edges for fixed-drive firing."""
    tf_dna: list[InteractionTuple] = []
    rest: list[InteractionTuple] = []
    for ix in interactions:
        src, tgt, _kind, _w, mode_s, _tid = ix
        if (
            mode_s == FireMode.REGULATORY.value
            and is_tf_dna_regulatory_edge(src, tgt, tf_genes)
        ):
            tf_dna.append(ix)
        else:
            rest.append(ix)
    return tf_dna, rest


def _fire_tf_dna_regulatory_independent(
    marking: np.ndarray,
    idx: dict[str, int],
    ixs: Sequence[InteractionTuple],
    rng: np.random.Generator,
    *,
    tf_dna_edge_weight: float,
    distribution: Distribution,
    min_tokens: float,
    ledger: dict[str, float] | None,
    fire_kwargs: Mapping[str, Any],
    use_geom: bool,
    fire_propensity_scale: float,
    fire_propensity_ref: float,
) -> None:
    """Each TF→DNA edge fires independently with δ ~ sample(m_TF × w), w = tf_dna_edge_weight."""
    w_fixed = float(tf_dna_edge_weight)
    if w_fixed <= 0.0 or not ixs:
        return
    order = rng.permutation(len(ixs))
    for j in order:
        src, tgt, kind, _weight, mode_s, _tid = ixs[int(j)]
        if use_geom:
            p = _fire_probability_geom_mean(
                marking,
                idx,
                src,
                min_tokens=min_tokens,
                scale=fire_propensity_scale,
                ref=fire_propensity_ref,
            )
            if float(rng.random()) >= p:
                continue
        before_tgt = float(marking[idx[tgt]]) if tgt in idx else 0.0
        _fire_one(
            marking,
            idx,
            src,
            tgt,
            kind,
            w_fixed,
            rng,
            distribution=distribution,
            min_tokens=min_tokens,
            fire_mode=FireMode(mode_s),
            ledger=ledger,
            **dict(fire_kwargs),
        )
        if ledger is not None:
            ledger["tf_dna_indep_events"] = float(
                ledger.get("tf_dna_indep_events", 0.0)
            ) + 1.0
            if tgt in idx:
                delta_flux = abs(float(marking[idx[tgt]]) - before_tgt)
                ledger["tf_dna_indep_flux"] = float(
                    ledger.get("tf_dna_indep_flux", 0.0)
                ) + delta_flux


def _categorical_pick(rng: np.random.Generator, weights: Sequence[float]) -> int:
    w = np.asarray(weights, dtype=np.float64)
    w = np.maximum(w, 0.0)
    total = float(w.sum())
    if total <= 0.0 or not np.isfinite(total):
        return int(rng.integers(0, len(w)))
    return int(rng.choice(len(w), p=w / total))


def _sample_first_order_route_amount(
    rng: np.random.Generator,
    available: float,
    *,
    half_life_blocks: float,
    distribution: Distribution,
) -> float:
    """First-order reaction budget: δ ~ f(m), same weight map as discrete half-life (not a sink)."""
    if available <= 0 or half_life_blocks <= 0:
        return 0.0
    w = first_order_turnover_weight(half_life_blocks)
    upper = float(available) * float(w)
    if upper <= 0 or not np.isfinite(upper):
        return 0.0
    hi = int(np.floor(upper))
    if hi < 1:
        # Fractional propensity: Bernoulli(upper) for one token.
        return 1.0 if float(rng.random()) < upper else 0.0
    return min(float(available), _sample_amount(rng, upper, distribution=distribution))


def _geom_mean_available(
    marking: np.ndarray,
    idx: dict[str, int],
    substrate_places: Sequence[str],
    *,
    min_tokens: float,
) -> float:
    """Geometric mean of available substrate marking; 0 if any substrate is depleted."""
    if not substrate_places:
        return 0.0
    log_sum = 0.0
    n = 0
    for place in substrate_places:
        if place not in idx:
            return 0.0
        available = max(0.0, float(marking[idx[place]]) - min_tokens)
        if available <= 0.0:
            return 0.0
        log_sum += math.log(available)
        n += 1
    return math.exp(log_sum / n)


def _geom_propensity_applies(fire_mode: FireMode) -> bool:
    """PPI (catalyst) and TF/histone regulatory edges only (phase 1)."""
    return fire_mode in (FireMode.CATALYST, FireMode.REGULATORY)


def _fire_probability_geom_mean(
    marking: np.ndarray,
    idx: dict[str, int],
    source: str,
    *,
    min_tokens: float,
    scale: float,
    ref: float,
) -> float:
    """p_fire = min(1, scale × geom_mean(substrates) / ref); unary → source token."""
    if scale <= 0.0 or ref <= 0.0:
        return 1.0
    gm = _geom_mean_available(marking, idx, [source], min_tokens=min_tokens)
    if gm <= 0.0:
        return 0.0
    return min(1.0, float(scale) * gm / float(ref))


def _fire_block_interactions(
    marking: np.ndarray,
    idx: dict[str, int],
    interactions: Sequence[InteractionTuple],
    rng: np.random.Generator,
    *,
    layer_rules: str,
    distribution: Distribution,
    min_tokens: float,
    ledger: dict[str, float] | None,
    prod_alpha_transcription: float,
    prod_alpha_translation: float,
    prod_m_max: float,
    translation_output_scale: float = 1.0,
    protein_production_output_scale: float = 1.0,
    tf_proteolysis_eps: float,
    activity_sigmoid_m_max: float,
    activity_sigmoid_k: float,
    activity_update_mode: str,
    activity_relu_drive_scale: float = 1.0,
    activity_relu_scale_mode: str = "drive",
    out_route_half_life_blocks: float,
    hme_regulatory_budget_scale: float = 1.0,
    hme_protein_places: frozenset[str] | set[str] | None = None,
    chromhmm_mark_weights: Mapping[str, float] | None = None,
    fire_propensity_mode: str = "off",
    fire_propensity_scale: float = 1.0,
    fire_propensity_ref: float = 10.0,
    regulatory_norm: str = REG_NORM_OFF,
    regulatory_hill_k: float = 5.0,
    regulatory_target_budget: float = 10.0,
    tf_dna_route_independent: bool = False,
    tf_dna_edge_weight: float = 0.5,
    tf_gene_symbols: frozenset[str] | set[str] | None = None,
    flux_events: list[dict[str, Any]] | None = None,
    flux_watch_targets: frozenset[str] | set[str] | None = None,
    current_block: int = -1,
) -> None:
    """One MC block: independent edges all fire; routed outs share f(m) flux on one random edge."""
    use_route = layer_rules == LAYER_RULES_V2 and float(out_route_half_life_blocks) > 0.0
    use_geom = fire_propensity_mode == FIRE_PROPENSITY_GEOM_MEAN
    marking_start = marking.copy()
    hme_kw = dict(
        hme_regulatory_budget_scale=hme_regulatory_budget_scale,
        hme_protein_places=hme_protein_places,
        chromhmm_mark_weights=chromhmm_mark_weights,
        activity_update_mode=activity_update_mode,
        activity_relu_drive_scale=activity_relu_drive_scale,
        activity_relu_scale_mode=activity_relu_scale_mode,
    )
    fire_common = dict(
        prod_alpha_transcription=prod_alpha_transcription,
        prod_alpha_translation=prod_alpha_translation,
        prod_m_max=prod_m_max,
        translation_output_scale=translation_output_scale,
        protein_production_output_scale=protein_production_output_scale,
        tf_proteolysis_eps=tf_proteolysis_eps,
        activity_sigmoid_m_max=activity_sigmoid_m_max,
        activity_sigmoid_k=activity_sigmoid_k,
        flux_events=flux_events,
        flux_watch_targets=flux_watch_targets,
        current_block=current_block,
        **hme_kw,
    )
    norm = str(regulatory_norm or REG_NORM_OFF).strip().lower()
    hub_ixs: list[InteractionTuple] = []
    work_ixs: Sequence[InteractionTuple] = interactions
    if norm in _REG_NORM_HUB_MODES:
        hub_ixs, rest = _split_hub_regulatory(interactions)
        work_ixs = rest

    tf_dna_ixs: list[InteractionTuple] = []
    tf_genes = frozenset(tf_gene_symbols or ())
    if tf_dna_route_independent and tf_genes and layer_rules == LAYER_RULES_V2:
        tf_dna_ixs, work_ixs = _split_tf_dna_regulatory(list(work_ixs), tf_genes)

    if tf_dna_ixs:
        _fire_tf_dna_regulatory_independent(
            marking,
            idx,
            tf_dna_ixs,
            rng,
            tf_dna_edge_weight=float(tf_dna_edge_weight),
            distribution=distribution,
            min_tokens=min_tokens,
            ledger=ledger,
            fire_kwargs=fire_common,
            use_geom=use_geom,
            fire_propensity_scale=fire_propensity_scale,
            fire_propensity_ref=fire_propensity_ref,
        )

    if not use_route:
        order = rng.permutation(len(work_ixs))
        for t_idx in order:
            src, tgt, kind, weight, mode_s, _tid = work_ixs[int(t_idx)]
            mode = FireMode(mode_s)
            if use_geom and _geom_propensity_applies(mode):
                p = _fire_probability_geom_mean(
                    marking,
                    idx,
                    src,
                    min_tokens=min_tokens,
                    scale=fire_propensity_scale,
                    ref=fire_propensity_ref,
                )
                if float(rng.random()) >= p:
                    continue
            _fire_one(
                marking,
                idx,
                src,
                tgt,
                kind,
                weight,
                rng,
                distribution=distribution,
                min_tokens=min_tokens,
                fire_mode=mode,
                ledger=ledger,
                **fire_common,
            )
    else:
        routed, independent = _partition_routed_interactions(work_ixs)
        if independent:
            order = rng.permutation(len(independent))
            for j in order:
                src, tgt, kind, weight, mode_s, _tid = independent[int(j)]
                mode = FireMode(mode_s)
                if use_geom and _geom_propensity_applies(mode):
                    p = _fire_probability_geom_mean(
                        marking,
                        idx,
                        src,
                        min_tokens=min_tokens,
                        scale=fire_propensity_scale,
                        ref=fire_propensity_ref,
                    )
                    if float(rng.random()) >= p:
                        continue
                _fire_one(
                    marking,
                    idx,
                    src,
                    tgt,
                    kind,
                    weight,
                    rng,
                    distribution=distribution,
                    min_tokens=min_tokens,
                    fire_mode=mode,
                    ledger=ledger,
                    **fire_common,
                )

        keys = list(routed.keys())
        rng.shuffle(keys)
        hl = float(out_route_half_life_blocks)
        for key in keys:
            ixs = routed[key]
            src = key[0]
            mode_s = key[1]
            if src not in idx:
                continue
            available = max(0.0, float(marking[idx[src]]) - min_tokens)
            if available <= 0:
                continue
            budget = _sample_first_order_route_amount(
                rng, available, half_life_blocks=hl, distribution=distribution
            )
            if budget <= 0:
                continue
            pick = _categorical_pick(rng, [float(ix[3]) for ix in ixs])
            _src, tgt, kind, weight, _mode, _tid = ixs[pick]
            _fire_one(
                marking,
                idx,
                src,
                tgt,
                kind,
                weight,
                rng,
                distribution=distribution,
                min_tokens=min_tokens,
                fire_mode=FireMode(mode_s),
                ledger=ledger,
                force_delta=budget,
                **fire_common,
            )
            if ledger is not None:
                ledger["out_route_events"] = float(ledger.get("out_route_events", 0.0)) + 1.0
                ledger["out_route_flux"] = float(ledger.get("out_route_flux", 0.0)) + float(
                    budget
                )

    if hub_ixs:
        _fire_hub_regulatory_target_softmax(
            marking,
            idx,
            hub_ixs,
            rng,
            regulatory_norm=norm,
            regulatory_hill_k=float(regulatory_hill_k),
            regulatory_target_budget=float(regulatory_target_budget),
            distribution=distribution,
            min_tokens=min_tokens,
            ledger=ledger,
            out_route_half_life_blocks=float(out_route_half_life_blocks),
            fire_kwargs=fire_common,
        )

    if str(activity_relu_scale_mode or "drive").strip().lower() == "net":
        n_adj = apply_promoter_active_relu_net_scale(
            marking,
            idx,
            list(idx.keys()),
            marking_start,
            mode=activity_update_mode,
            scale=activity_relu_drive_scale,
            min_tokens=min_tokens,
        )
        if ledger is not None and n_adj > 0:
            ledger["relu_net_scale_places"] = float(
                ledger.get("relu_net_scale_places", 0.0)
            ) + float(n_adj)


def _run_one_trajectory(payload: dict[str, Any]) -> tuple[np.ndarray, dict[str, float], list[dict[str, Any]]]:
    """Module-level worker for ProcessPool (must be picklable)."""
    places: list[str] = payload["places"]
    interactions: list[InteractionTuple] = payload["interactions"]
    fixed_high: set[str] = set(payload["fixed_high"])
    fixed_low: set[str] = set(payload["fixed_low"])
    n_blocks: int = payload["n_blocks"]
    high_tokens: int = payload["high_tokens"]
    low_tokens: int = payload["low_tokens"]
    distribution: Distribution = payload["distribution"]
    min_tokens: float = payload["min_tokens"]
    maintain_seeds: bool = payload["maintain_seeds"]
    maintain_every_n_blocks: int = max(1, int(payload.get("maintain_every_n_blocks") or 1))
    maintain_places: set[str] = set(payload.get("maintain_places") or ())
    maintain_levels: dict[str, float] = dict(payload.get("maintain_levels") or {})
    initial_marking: dict[str, float] | None = payload.get("initial_marking")
    run_seed: int = payload["run_seed"]
    layer_rules: str = payload.get("layer_rules") or LAYER_RULES_OFF
    prod_alpha_transcription: float = float(payload.get("prod_alpha_transcription", 1.0))
    prod_alpha_translation: float = float(payload.get("prod_alpha_translation", 1.0))
    prod_m_max: float = float(payload.get("prod_m_max", 50.0))
    translation_output_scale: float = float(payload.get("translation_output_scale", 1.0))
    protein_production_output_scale: float = float(
        payload.get("protein_production_output_scale", 1.0)
    )
    tf_proteolysis_eps: float = float(payload.get("tf_proteolysis_eps", 0.0))
    activity_sigmoid_m_max: float = float(
        payload.get("activity_sigmoid_m_max", DEFAULT_ACTIVITY_SIGMOID_M_MAX)
    )
    activity_sigmoid_k: float = float(
        payload.get("activity_sigmoid_k", DEFAULT_ACTIVITY_SIGMOID_K)
    )
    activity_update_mode: str = str(
        payload.get("activity_update_mode", ACTIVITY_UPDATE_SIGMOID) or ACTIVITY_UPDATE_SIGMOID
    )
    activity_relu_drive_scale: float = float(
        payload.get("activity_relu_drive_scale", 1.0)
    )
    activity_relu_scale_mode: str = str(
        payload.get("activity_relu_scale_mode", "drive") or "drive"
    )
    fire_propensity_mode = str(payload.get("fire_propensity_mode") or "off")
    fire_propensity_scale = float(payload.get("fire_propensity_scale", 1.0))
    fire_propensity_ref = float(payload.get("fire_propensity_ref", 10.0))
    regulatory_norm = str(payload.get("regulatory_norm") or REG_NORM_OFF)
    regulatory_hill_k = float(payload.get("regulatory_hill_k", 5.0))
    regulatory_target_budget = float(payload.get("regulatory_target_budget", 10.0))
    tf_dna_route_independent = bool(payload.get("tf_dna_route_independent", False))
    tf_dna_edge_weight = float(
        payload.get("tf_dna_edge_weight", payload.get("tf_dna_fixed_drive", 0.5))
    )
    tf_gene_symbols = frozenset(payload.get("tf_gene_symbols") or ())
    out_route_half_life_blocks: float = float(
        payload.get("out_route_half_life_blocks", DEFAULT_OUT_ROUTE_HALF_LIFE_BLOCKS)
    )
    hme_regulatory_budget_scale: float = float(
        payload.get("hme_regulatory_budget_scale", 1.0)
    )
    hme_protein_places = frozenset(payload.get("hme_protein_places") or ())
    chromhmm_mark_weights = dict(payload.get("chromhmm_mark_weights") or {}) or None
    weight_pulse = InteractionWeightPulse.from_dict(payload.get("interaction_weight_pulse"))
    flux_watch_raw = payload.get("flux_watch_targets")
    flux_events: list[dict[str, Any]] | None = (
        [] if flux_watch_raw else None
    )
    flux_watch_targets = frozenset(flux_watch_raw or ())

    idx = {p: i for i, p in enumerate(places)}
    n_places = len(places)
    rng = np.random.default_rng(run_seed)
    ledger = empty_ledger()

    def _level_for(place: str) -> float:
        if place in maintain_levels:
            level = float(maintain_levels[place])
        else:
            level = float(high_tokens)
        if layer_rules == LAYER_RULES_V2 and is_sigmoid_activity_place(place):
            _mode, eff_m = resolve_place_activity_update(
                place, activity_update_mode, m_max=float(activity_sigmoid_m_max)
            )
            if eff_m > 0:
                level = min(level, float(eff_m))
        return level

    def _inject_seeds(m: np.ndarray) -> None:
        for p in fixed_high:
            if p not in idx:
                continue
            i = idx[p]
            level = _level_for(p)
            before = float(m[i])
            m[i] = level
            if level > before:
                ledger["seed_inject"] += level - before
        for p in fixed_low:
            if p in idx:
                m[idx[p]] = low_tokens

    def _inject_places(m: np.ndarray, place_names: AbstractSet[str]) -> None:
        for p in place_names:
            if p not in idx:
                continue
            i = idx[p]
            level = _level_for(p)
            before = float(m[i])
            m[i] = level
            if level > before:
                ledger["seed_inject"] += level - before

    m = np.zeros(n_places, dtype=np.float64)
    _inject_seeds(m)
    if initial_marking:
        for p, val in initial_marking.items():
            if p in idx:
                m[idx[p]] = val

    out = np.zeros((n_blocks + 1, n_places), dtype=np.float64)
    out[0] = m
    for block in range(n_blocks):
        if block % maintain_every_n_blocks == 0:
            if maintain_seeds:
                _inject_seeds(m)
            elif maintain_places:
                _inject_places(m, maintain_places)
        block_ix = scale_interaction_tuples_for_block(interactions, weight_pulse, block)
        _fire_block_interactions(
            m,
            idx,
            block_ix,
            rng,
            layer_rules=layer_rules,
            distribution=distribution,
            min_tokens=min_tokens,
            ledger=ledger,
            prod_alpha_transcription=prod_alpha_transcription,
            prod_alpha_translation=prod_alpha_translation,
            prod_m_max=prod_m_max,
            translation_output_scale=translation_output_scale,
            protein_production_output_scale=protein_production_output_scale,
            tf_proteolysis_eps=tf_proteolysis_eps,
            activity_sigmoid_m_max=activity_sigmoid_m_max,
            activity_sigmoid_k=activity_sigmoid_k,
            activity_update_mode=activity_update_mode,
            activity_relu_drive_scale=activity_relu_drive_scale,
            activity_relu_scale_mode=activity_relu_scale_mode,
            out_route_half_life_blocks=out_route_half_life_blocks,
            hme_regulatory_budget_scale=hme_regulatory_budget_scale,
            hme_protein_places=hme_protein_places,
            chromhmm_mark_weights=chromhmm_mark_weights,
            fire_propensity_mode=fire_propensity_mode,
            fire_propensity_scale=fire_propensity_scale,
            fire_propensity_ref=fire_propensity_ref,
            regulatory_norm=regulatory_norm,
            regulatory_hill_k=regulatory_hill_k,
            regulatory_target_budget=regulatory_target_budget,
            tf_dna_route_independent=tf_dna_route_independent,
            tf_dna_edge_weight=tf_dna_edge_weight,
            tf_gene_symbols=tf_gene_symbols,
            flux_events=flux_events,
            flux_watch_targets=flux_watch_targets,
            current_block=block,
        )
        out[block + 1] = m
    return out, ledger, flux_events or []


def _interactions_as_tuples(
    interactions: Sequence[Any],
    *,
    layer_rules: str = LAYER_RULES_OFF,
) -> list[InteractionTuple]:
    out: list[InteractionTuple] = []
    for edge in interactions:
        src = str(edge.source)
        tgt = str(edge.target)
        kind = str(edge.kind)
        weight = float(getattr(edge, "weight", 1.0) or 1.0)
        tid = str(getattr(edge, "transition_id", "") or "")
        mode = resolve_fire_mode(
            src, tgt, kind, transition_id=tid, layer_rules=layer_rules
        )
        out.append((src, tgt, kind, weight, mode.value, tid))
    return out


@dataclass
class GINSPNSimulator:
    """
    Ruths-style block MC with GIN firing.

    Default (layer_rules=off): legacy gin_consume_v1.
    With layer_rules=v2: dogma production, regulatory TF/histone, conserving transfer.
    """

    n_blocks: int = 200
    n_runs: int = 200
    high_tokens: int = 10
    low_tokens: int = 0
    distribution: Distribution = "uniform"
    seed: int | None = 42
    min_tokens: float = 0.0
    maintain_seeds: bool = True
    """If True, re-inject fixed_high places on a schedule (see maintain_every_n_blocks)."""
    maintain_every_n_blocks: int = 1
    """Clamp seeds every N blocks (1=every block ≈ full rate; 2=half injection frequency)."""
    maintain_places: frozenset[str] | None = None
    """Optional subset to clamp each block when maintain_seeds=False (e.g. met sources)."""
    maintain_levels: dict[str, float] | None = None
    """Optional per-place token levels (e.g. expression-scaled); overrides high_tokens."""
    n_jobs: int = 20
    """Parallel workers for independent MC runs (multiprocess; 1 = serial)."""
    layer_rules: str = LAYER_RULES_OFF
    prod_alpha_transcription: float = 1.0
    prod_alpha_translation: float = 1.0
    prod_m_max: float = 50.0
    translation_output_scale: float = 1.0
    protein_production_output_scale: float = 1.0
    tf_proteolysis_eps: float = 0.0
    activity_sigmoid_m_max: float = DEFAULT_ACTIVITY_SIGMOID_M_MAX
    activity_sigmoid_k: float = DEFAULT_ACTIVITY_SIGMOID_K
    activity_update_mode: str = ACTIVITY_UPDATE_SIGMOID
    activity_relu_drive_scale: float = 1.0
    activity_relu_scale_mode: str = "drive"
    out_route_half_life_blocks: float = DEFAULT_OUT_ROUTE_HALF_LIFE_BLOCKS
    """Effective first-order half-life (blocks) for routed TRANSFER/REGULATORY out-flux f(m)."""
    hme_regulatory_budget_scale: float = 1.0
    """Multiply REGULATORY δ when source ∈ HME.protein (enzyme not consumed)."""
    hme_protein_places: frozenset[str] | None = None
    """If None and scale≠1, load defaults from nodeList_HME at run()."""
    chromhmm_mark_weights: dict[str, float] | None = None
    """Optional signed mark weights (normalized keys) for hist→promoter.active (SPN v8)."""
    fire_propensity_mode: str = "off"
    """off | geom_mean — Bernoulli skip for CATALYST/REGULATORY before _fire_one."""
    fire_propensity_scale: float = 1.0
    fire_propensity_ref: float = 10.0
    """Reference token level for p = min(1, scale×geom_mean/ref). Default = high_tokens."""
    interaction_weight_pulse: InteractionWeightPulse | dict | None = None
    """Optional MDM2⊣TP53-style weight pulse (dict or InteractionWeightPulse)."""
    nf_qc_pair: tuple[str, str] | None = None
    """If set, compute per-run lagged phase for this place pair before discarding samples."""
    nf_warm: int = 20
    nf_max_lag: int = 40
    regulatory_norm: str = REG_NORM_OFF
    """off | target_softmax | hill_target_softmax — promoter hub in-edge competition."""
    regulatory_hill_k: float = 5.0
    """Michaelis K for hill_target_softmax occupancy m/(K+m)."""
    regulatory_target_budget: float = 10.0
    """Shared REGULATORY drive budget tokens per (promoter hub, kind) per block."""
    tf_dna_route_independent: bool = True
    """TF.protein→promoter hub: independent REGULATORY fire, exempt from out-route fan-out."""
    tf_dna_edge_weight: float = 0.5
    """Fixed multiplicative edge weight w for TF→DNA REGULATORY (δ ~ sample(m_TF × w))."""
    tf_gene_symbols: frozenset[str] | None = None
    """TF gene symbols; None → load nodeList_TF at run()."""

    def _sample_amount(self, rng: np.random.Generator, upper: float) -> float:
        return _sample_amount(rng, upper, distribution=self.distribution)

    def _level_for(self, place: str) -> float:
        levels = self.maintain_levels or {}
        if place in levels:
            level = float(levels[place])
        else:
            level = float(self.high_tokens)
        if (
            self.layer_rules == LAYER_RULES_V2
            and is_sigmoid_activity_place(place)
        ):
            _mode, eff_m = resolve_place_activity_update(
                place, self.activity_update_mode, m_max=float(self.activity_sigmoid_m_max)
            )
            if eff_m > 0:
                level = min(level, float(eff_m))
        return level

    def _apply_place_marking(
        self,
        marking: np.ndarray,
        places: AbstractSet[str],
        idx: dict[str, int],
        ledger: dict[str, float] | None = None,
    ) -> None:
        for p in places:
            if p not in idx:
                continue
            i = idx[p]
            level = self._level_for(p)
            before = float(marking[i])
            marking[i] = level
            if ledger is not None and level > before:
                ledger["seed_inject"] += level - before

    def _apply_seed_marking(
        self,
        marking: np.ndarray,
        model,
        idx: dict[str, int],
        ledger: dict[str, float] | None = None,
    ) -> None:
        for p in model.fixed_high:
            if p not in idx:
                continue
            i = idx[p]
            level = self._level_for(p)
            before = float(marking[i])
            marking[i] = level
            if ledger is not None and level > before:
                ledger["seed_inject"] += level - before
        for p in model.fixed_low:
            if p in idx:
                marking[idx[p]] = self.low_tokens

    def _initial_marking(
        self,
        model,
        places: list[str],
        idx: dict[str, int],
        ledger: dict[str, float] | None = None,
    ) -> np.ndarray:
        m = np.zeros(len(places), dtype=np.float64)
        self._apply_seed_marking(m, model, idx, ledger=ledger)
        return m

    def _fire_transition(
        self,
        marking: np.ndarray,
        idx: dict[str, int],
        edge,
        rng: np.random.Generator,
        ledger: dict[str, float] | None = None,
    ) -> None:
        tid = str(getattr(edge, "transition_id", "") or "")
        mode = resolve_fire_mode(
            edge.source,
            edge.target,
            edge.kind,
            transition_id=tid,
            layer_rules=self.layer_rules,
        )
        _fire_one(
            marking,
            idx,
            edge.source,
            edge.target,
            edge.kind,
            float(getattr(edge, "weight", 1.0) or 1.0),
            rng,
            distribution=self.distribution,
            min_tokens=self.min_tokens,
            fire_mode=mode,
            ledger=ledger,
            prod_alpha_transcription=self.prod_alpha_transcription,
            prod_alpha_translation=self.prod_alpha_translation,
            prod_m_max=self.prod_m_max,
            translation_output_scale=self.translation_output_scale,
            protein_production_output_scale=self.protein_production_output_scale,
            tf_proteolysis_eps=self.tf_proteolysis_eps,
            activity_sigmoid_m_max=self.activity_sigmoid_m_max,
            activity_sigmoid_k=self.activity_sigmoid_k,
            activity_update_mode=self.activity_update_mode,
            activity_relu_drive_scale=self.activity_relu_drive_scale,
            activity_relu_scale_mode=self.activity_relu_scale_mode,
        )

    def _run_seed_for(self, run: int) -> int:
        base = 0 if self.seed is None else int(self.seed)
        return base + 1_000_003 * run + 17

    @property
    def firing_policy_id(self) -> str:
        base = GIN_FIRING_POLICY_V2 if self.layer_rules == LAYER_RULES_V2 else GIN_FIRING_POLICY
        if self.fire_propensity_mode == FIRE_PROPENSITY_GEOM_MEAN:
            return f"{base}+{FIRE_PROPENSITY_GEOM_MEAN}"
        return base

    def run(self, model, *, initial_marking: dict[str, float] | None = None, stimulus=None) -> SimulationResult:
        places = list(model.places)
        interactions = list(model.interactions)
        if not interactions:
            raise ValueError("SPN model has no interactions")

        n_places = len(places)
        n_jobs = max(1, int(self.n_jobs))
        inter_tuples = _interactions_as_tuples(interactions, layer_rules=self.layer_rules)
        fixed_high = list(getattr(model, "fixed_high", set()) or [])
        fixed_low = list(getattr(model, "fixed_low", set()) or [])
        hme_scale = float(self.hme_regulatory_budget_scale)
        if self.hme_protein_places is not None:
            hme_places = frozenset(self.hme_protein_places)
        elif hme_scale != 1.0:
            hme_places = default_hme_protein_places()
        else:
            hme_places = frozenset()
        emw = dict(self.chromhmm_mark_weights or {})
        tf_genes = self.tf_gene_symbols
        if tf_genes is None and self.tf_dna_route_independent:
            from gin_tool.active_tf import load_tf_symbol_pool

            root = Path(__file__).resolve().parents[2]
            tf_genes = frozenset(load_tf_symbol_pool(root / "nodeList"))
        elif tf_genes is not None:
            tf_genes = frozenset(tf_genes)
        else:
            tf_genes = frozenset()
        pulse = self.interaction_weight_pulse
        if isinstance(pulse, InteractionWeightPulse):
            pulse_payload = pulse.to_dict()
        elif isinstance(pulse, dict):
            pulse_payload = dict(pulse)
        else:
            pulse_payload = None

        levels = dict(self.maintain_levels or {})
        maintain_places = sorted(self.maintain_places or ())
        payloads = [
            {
                "places": places,
                "interactions": inter_tuples,
                "fixed_high": fixed_high,
                "fixed_low": fixed_low,
                "n_blocks": self.n_blocks,
                "high_tokens": self.high_tokens,
                "low_tokens": self.low_tokens,
                "distribution": self.distribution,
                "min_tokens": self.min_tokens,
                "maintain_seeds": self.maintain_seeds,
                "maintain_every_n_blocks": max(1, int(self.maintain_every_n_blocks)),
                "maintain_places": maintain_places,
                "maintain_levels": levels,
                "initial_marking": initial_marking,
                "run_seed": self._run_seed_for(run),
                "layer_rules": self.layer_rules,
                "prod_alpha_transcription": self.prod_alpha_transcription,
                "prod_alpha_translation": self.prod_alpha_translation,
                "prod_m_max": self.prod_m_max,
                "translation_output_scale": self.translation_output_scale,
                "protein_production_output_scale": self.protein_production_output_scale,
                "tf_proteolysis_eps": self.tf_proteolysis_eps,
                "activity_sigmoid_m_max": self.activity_sigmoid_m_max,
                "activity_sigmoid_k": self.activity_sigmoid_k,
                "activity_update_mode": self.activity_update_mode,
                "activity_relu_drive_scale": self.activity_relu_drive_scale,
                "activity_relu_scale_mode": self.activity_relu_scale_mode,
                "fire_propensity_mode": self.fire_propensity_mode,
                "fire_propensity_scale": float(self.fire_propensity_scale),
                "fire_propensity_ref": float(
                    self.fire_propensity_ref
                    if self.fire_propensity_ref > 0
                    else float(self.high_tokens)
                ),
                "out_route_half_life_blocks": self.out_route_half_life_blocks,
                "hme_regulatory_budget_scale": hme_scale,
                "hme_protein_places": list(hme_places),
                "chromhmm_mark_weights": emw,
                "interaction_weight_pulse": pulse_payload,
                "regulatory_norm": self.regulatory_norm,
                "regulatory_hill_k": float(self.regulatory_hill_k),
                "regulatory_target_budget": float(self.regulatory_target_budget),
                "tf_dna_route_independent": bool(self.tf_dna_route_independent),
                "tf_dna_edge_weight": float(self.tf_dna_edge_weight),
                "tf_gene_symbols": list(tf_genes),
            }
            for run in range(self.n_runs)
        ]

        block_samples = np.zeros((self.n_runs, self.n_blocks + 1, n_places), dtype=np.float64)
        ledger_sum = empty_ledger()
        if n_jobs == 1:
            rng = np.random.default_rng(self.seed)
            idx = {p: i for i, p in enumerate(places)}
            weight_pulse = InteractionWeightPulse.from_dict(pulse_payload)
            for run in range(self.n_runs):
                ledger = empty_ledger()
                m = self._initial_marking(model, places, idx, ledger=ledger)
                if initial_marking:
                    for p, val in initial_marking.items():
                        if p in idx:
                            m[idx[p]] = val
                block_samples[run, 0] = m.copy()
                for block in range(self.n_blocks):
                    every = max(1, int(self.maintain_every_n_blocks))
                    if block % every == 0:
                        if self.maintain_seeds:
                            self._apply_seed_marking(m, model, idx, ledger=ledger)
                        elif self.maintain_places:
                            self._apply_place_marking(
                                m, self.maintain_places, idx, ledger=ledger
                            )
                    block_ix = scale_interaction_tuples_for_block(
                        inter_tuples, weight_pulse, block
                    )
                    _fire_block_interactions(
                        m,
                        idx,
                        block_ix,
                        rng,
                        layer_rules=self.layer_rules,
                        distribution=self.distribution,
                        min_tokens=self.min_tokens,
                        ledger=ledger,
                        prod_alpha_transcription=self.prod_alpha_transcription,
                        prod_alpha_translation=self.prod_alpha_translation,
                        prod_m_max=self.prod_m_max,
                        translation_output_scale=self.translation_output_scale,
                        protein_production_output_scale=self.protein_production_output_scale,
                        tf_proteolysis_eps=self.tf_proteolysis_eps,
                        activity_sigmoid_m_max=self.activity_sigmoid_m_max,
                        activity_sigmoid_k=self.activity_sigmoid_k,
                        activity_update_mode=self.activity_update_mode,
                        activity_relu_drive_scale=self.activity_relu_drive_scale,
                        activity_relu_scale_mode=self.activity_relu_scale_mode,
                        out_route_half_life_blocks=self.out_route_half_life_blocks,
                        hme_regulatory_budget_scale=hme_scale,
                        hme_protein_places=hme_places,
                        chromhmm_mark_weights=emw or None,
                        fire_propensity_mode=self.fire_propensity_mode,
                        fire_propensity_scale=float(self.fire_propensity_scale),
                        fire_propensity_ref=float(
                            self.fire_propensity_ref
                            if self.fire_propensity_ref > 0
                            else float(self.high_tokens)
                        ),
                        regulatory_norm=self.regulatory_norm,
                        regulatory_hill_k=float(self.regulatory_hill_k),
                        regulatory_target_budget=float(self.regulatory_target_budget),
                        tf_dna_route_independent=bool(self.tf_dna_route_independent),
                        tf_dna_edge_weight=float(self.tf_dna_edge_weight),
                        tf_gene_symbols=tf_genes,
                    )
                    block_samples[run, block + 1] = m.copy()
                for k, v in ledger.items():
                    ledger_sum[k] = float(ledger_sum.get(k, 0.0)) + float(v)
        else:
            workers = min(n_jobs, self.n_runs)
            print(
                f"[spn] ProcessPoolExecutor workers={workers} runs={self.n_runs} "
                f"blocks={self.n_blocks} places={n_places} layer_rules={self.layer_rules}",
                flush=True,
            )
            with ProcessPoolExecutor(max_workers=workers) as pool:
                for run, (traj, ledger, _flux) in enumerate(
                    pool.map(_run_one_trajectory, payloads, chunksize=1)
                ):
                    block_samples[run] = traj
                    for k, v in ledger.items():
                        ledger_sum[k] = float(ledger_sum.get(k, 0.0)) + float(v)

        mean_m = block_samples.mean(axis=0)
        std_m = block_samples.std(axis=0, ddof=1) if self.n_runs > 1 else np.zeros_like(mean_m)

        nf_passrate = None
        nf_per_run = None
        pair_block_samples = None
        pair_place_names = None
        if self.nf_qc_pair is not None:
            from gin_tool.subgraph_pipeline.spn_nf_phase import (
                per_run_phase_from_block_samples,
                summarize_nf_passrate,
            )

            a_name, b_name = self.nf_qc_pair
            ia = places.index(a_name)
            ib = places.index(b_name)
            pair_place_names = (a_name, b_name)
            pair_block_samples = np.ascontiguousarray(
                block_samples[:, :, [ia, ib]], dtype=np.float64
            )
            nf_per_run = per_run_phase_from_block_samples(
                block_samples,
                places,
                self.nf_qc_pair,
                warm=int(self.nf_warm),
                max_lag=int(self.nf_max_lag),
            )
            nf_passrate = summarize_nf_passrate(nf_per_run)

        return SimulationResult(
            place_names=places,
            mean_markings=mean_m,
            std_markings=std_m,
            n_runs=self.n_runs,
            n_blocks=self.n_blocks,
            ledger_totals=ledger_sum,
            firing_policy=self.firing_policy_id,
            nf_passrate=nf_passrate,
            nf_per_run=nf_per_run,
            pair_block_samples=pair_block_samples,
            pair_place_names=pair_place_names,
        )

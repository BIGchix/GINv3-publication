"""H1 B0 standard SPN policy (fast×e + e3ewe + uncapped dogma + met-source maintain).

Frozen from Job 382 ``v15cfea10e3ewe_o150b350_metsrc`` (user-validated 2026-08-23).
"""

from __future__ import annotations

import math
from typing import Any

from gin_tool.subgraph_pipeline.spn_layer_rules import LAYER_RULES_V2
from gin_tool.subgraph_pipeline.spn_policy import (
    SPN_PIPELINE_VERSION_STD,
    SpnSimPolicy,
    WEAK_E3_EDGE_SCALE_BASELINE,
)

MATH_E = math.e
STD_FAST_E_TIMESCALE = MATH_E

# v15c reference half-lives (blocks) before fast×e.
STD_V15C_PROTEIN_HL = 10.0
STD_V15C_RNA_HL = 5.0
STD_V15C_PROMOTER_ACTIVE_HL = 5.0
STD_V15C_OUT_ROUTE_HL = 20.0

STD_PROMOTER_INCOMING_SCALE = 0.1
STD_PROD_M_MAX = 0.0
STD_ALPHA_TX = 1.0
STD_STRONG_E3 = MATH_E
STD_WEAK_E3 = WEAK_E3_EDGE_SCALE_BASELINE  # 1/e
STD_TF_DNA_EDGE_W_BASE = 0.5
STD_TRANSLATION_SCALE = 1.0
STD_PROTEIN_PROD_SCALE = 1.0

STD_ACTIVITY_UPDATE_MODE = "relu_active"
STD_ACTIVITY_RELU_DRIVE_SCALE = 0.2
STD_ACTIVITY_RELU_SCALE_MODE = "net"
STD_ACTIVITY_SIGMOID_M_MAX = 40.0
STD_HISTONE_TO_PROMOTER_STATE_SCALE = 0.1

# B0 MDM2⊣TP53 pulse (plateau onset).
STD_B0_BLOCKS = 350
STD_B0_RUNS = 350
STD_B0_ONSET = 150
STD_B0_HOLD = 40
STD_B0_PULSE_LOW_FACTOR = 0.1
STD_B0_PULSE_RAMP_BLOCKS = 20


def fast_e_half_life(v15c_blocks: float) -> float:
    """Scale a v15c half-life by fast×e (÷e)."""
    return float(v15c_blocks) / STD_FAST_E_TIMESCALE


def standard_h1_b0_half_lives() -> dict[str, float]:
    return {
        "protein_half_life_blocks": fast_e_half_life(STD_V15C_PROTEIN_HL),
        "rna_half_life_blocks": fast_e_half_life(STD_V15C_RNA_HL),
        "promoter_active_half_life_blocks": fast_e_half_life(STD_V15C_PROMOTER_ACTIVE_HL),
        "out_route_half_life_blocks": fast_e_half_life(STD_V15C_OUT_ROUTE_HL),
    }


def standard_h1_b0_spn_policy(
    *,
    chromhmm_emission_tsv: str | None = None,
    interaction_weight_pulse: dict | None = None,
    blocks: int = STD_B0_BLOCKS,
    runs: int = STD_B0_RUNS,
    **overrides: Any,
) -> SpnSimPolicy:
    """Return the validated H1 B0 standard ``SpnSimPolicy`` (override knobs sparingly)."""
    hl = standard_h1_b0_half_lives()
    base = dict(
        version=SPN_PIPELINE_VERSION_STD,
        layer_rules=LAYER_RULES_V2,
        blocks=int(blocks),
        runs=int(runs),
        maintain_seeds=False,
        maintain_metabolite_sources=True,
        add_degradation=True,
        mark_mode="ancestor_all",
        prod_m_max=STD_PROD_M_MAX,
        prod_alpha_transcription=STD_ALPHA_TX,
        prod_alpha_translation=1.0,
        translation_output_scale=STD_TRANSLATION_SCALE,
        protein_production_output_scale=STD_PROTEIN_PROD_SCALE,
        promoter_incoming_edge_scale=STD_PROMOTER_INCOMING_SCALE,
        strong_e3_inhibition_weight=STD_STRONG_E3,
        weak_e3_edge_scale=STD_WEAK_E3,
        tf_dna_edge_weight=STD_TF_DNA_EDGE_W_BASE,
        tf_dna_route_independent=True,
        histone_to_promoter_state_scale=STD_HISTONE_TO_PROMOTER_STATE_SCALE,
        drop_hme_to_promoter_state=True,
        sign_repressive_chromatin_to_state=True,
        ups_degradation_policy=True,
        scale_non_e3_activity_inhibitors=True,
        activity_update_mode=STD_ACTIVITY_UPDATE_MODE,
        activity_relu_drive_scale=STD_ACTIVITY_RELU_DRIVE_SCALE,
        activity_relu_scale_mode=STD_ACTIVITY_RELU_SCALE_MODE,
        activity_sigmoid_m_max=STD_ACTIVITY_SIGMOID_M_MAX,
        chromhmm_emission_tsv=chromhmm_emission_tsv,
        interaction_weight_pulse=interaction_weight_pulse,
        **hl,
    )
    base.update(overrides)
    return SpnSimPolicy(**base)


def standard_b0_pulse_step(
    *,
    onset: int = STD_B0_ONSET,
    hold: int = STD_B0_HOLD,
) -> dict:
    """MDM2.protein ⊣ TP53.protein inhibition pulse schedule."""
    return {
        "source": "MDM2.protein",
        "target": "TP53.protein",
        "kind": "inhibition",
        "onset": int(onset),
        "hold": int(hold),
        "low_factor": STD_B0_PULSE_LOW_FACTOR,
        "restore": "step",
        "ramp_blocks": STD_B0_PULSE_RAMP_BLOCKS,
        "base_factor": 1.0,
    }


def standard_h1_b0_simulate_kwargs(
    *,
    chromhmm_emission_tsv: str,
    interaction_weight_pulse: dict | None = None,
    blocks: int = STD_B0_BLOCKS,
    runs: int = STD_B0_RUNS,
    n_jobs: int = 20,
    seed: int = 42,
    **overrides: Any,
) -> dict[str, Any]:
    """Keyword args for ``export_and_simulate_subgraph_dir`` under the standard policy."""
    pol = standard_h1_b0_spn_policy(
        chromhmm_emission_tsv=chromhmm_emission_tsv,
        interaction_weight_pulse=interaction_weight_pulse,
        blocks=blocks,
        runs=runs,
    )
    hl = standard_h1_b0_half_lives()
    kw: dict[str, Any] = dict(
        policy=pol,
        layer_rules=pol.layer_rules,
        edge_balance="metabolite",
        maintain_seeds=pol.maintain_seeds,
        maintain_metabolite_sources=pol.maintain_metabolite_sources,
        out_route_half_life_blocks=hl["out_route_half_life_blocks"],
        rna_half_life_blocks=hl["rna_half_life_blocks"],
        protein_half_life_blocks=hl["protein_half_life_blocks"],
        promoter_active_half_life_blocks=hl["promoter_active_half_life_blocks"],
        prod_alpha_transcription=pol.prod_alpha_transcription,
        prod_m_max=pol.prod_m_max,
        translation_output_scale=pol.translation_output_scale,
        protein_production_output_scale=pol.protein_production_output_scale,
        promoter_incoming_edge_scale=pol.promoter_incoming_edge_scale,
        strong_e3_inhibition_weight=pol.strong_e3_inhibition_weight,
        weak_e3_edge_scale=pol.weak_e3_edge_scale,
        tf_dna_edge_weight=pol.tf_dna_edge_weight,
        chromhmm_emission_tsv=chromhmm_emission_tsv,
        activity_update_mode=pol.activity_update_mode,
        activity_sigmoid_m_max=pol.activity_sigmoid_m_max,
        activity_relu_drive_scale=pol.activity_relu_drive_scale,
        activity_relu_scale_mode=pol.activity_relu_scale_mode,
        histone_to_promoter_state_scale=pol.histone_to_promoter_state_scale,
        drop_chromatin_enzyme_to_promoter_state=pol.drop_hme_to_promoter_state,
        sign_repressive_chromatin_to_state=pol.sign_repressive_chromatin_to_state,
        ups_degradation_policy=pol.ups_degradation_policy,
        scale_non_e3_activity_inhibitors=pol.scale_non_e3_activity_inhibitors,
        blocks=blocks,
        runs=runs,
        n_jobs=n_jobs,
        seed=seed,
        interaction_weight_pulse=interaction_weight_pulse,
    )
    kw.update(overrides)
    return kw


def is_standard_h1_b0_knobs(
    *,
    prod_m_max: float,
    prod_alpha_transcription: float,
    promoter_incoming_edge_scale: float,
    strong_e3_inhibition_weight: float,
    weak_e3_edge_scale: float,
    protein_half_life_blocks: float,
    rtol: float = 0.02,
) -> bool:
    """True when dogma/E3/promoter knobs match the standard H1 B0 preset."""
    def _close(a: float, b: float) -> bool:
        return abs(float(a) - float(b)) <= rtol * max(abs(float(b)), 1e-9)

    return (
        float(prod_m_max) <= 0.0
        and _close(prod_alpha_transcription, STD_ALPHA_TX)
        and _close(promoter_incoming_edge_scale, STD_PROMOTER_INCOMING_SCALE)
        and _close(strong_e3_inhibition_weight, STD_STRONG_E3)
        and _close(weak_e3_edge_scale, STD_WEAK_E3)
        and _close(protein_half_life_blocks, fast_e_half_life(STD_V15C_PROTEIN_HL))
    )

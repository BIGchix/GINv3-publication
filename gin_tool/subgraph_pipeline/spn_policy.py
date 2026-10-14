"""Default SPN simulation policy for subgraph_pipeline (versioned)."""

from __future__ import annotations

import math
from dataclasses import dataclass

# Baseline weak UPS E3 edge scale (was 0.3; fastE baseline uses 1/e).
WEAK_E3_EDGE_SCALE_BASELINE = 1.0 / math.e

from gin_tool.subgraph_pipeline.spn_layer_rules import LAYER_RULES_OFF, LAYER_RULES_V2

# v5: degradation on; maintain tokens on 祖先-only sources (Petri in-degree 0).
# v6: optional layer_rules=v2 (dogma production / regulatory / conserving transfer + sinks).
# v7: out-route f(m), miRNA translation repress, met in-degree-0 source seeding helpers.
# v8: v7 + chromHMM emission → signed global mark weights on hist→promoter.active
#     (w_j = mean_A e_ij − mean_I e_ij; promoter states only; not locus-specific).
SPN_PIPELINE_VERSION = "spn_v5_degrade_anc_source"
SPN_PIPELINE_VERSION_V2 = "spn_v7_outroute_mir_metsrc"
SPN_PIPELINE_VERSION_V8 = "spn_v8_chromhmm_emw"
# v9: default UPS strong/weak ESI policy + non-E3 activity-inhibitor gene-list scaling.
SPN_PIPELINE_VERSION_V9 = "spn_v9_ups_non_e3"
# v10: v9 + drop HME.protein→promoter;state (HME only via histone marks).
SPN_PIPELINE_VERSION_V10 = "spn_v10_hme_state_drop"
# v11: v10 + hist×0.1 + repressive sign + TF→DNA fixed per-block drive (superseded).
SPN_PIPELINE_VERSION_V11 = "spn_v11_tf_dna_fixed"
# v12: v10 + hist×0.1 + repressive sign + TF→DNA fixed edge weight w (δ ∝ m_TF×w, no out-route).
SPN_PIPELINE_VERSION_V12 = "spn_v12_tf_dna_edge"
# v13: v12 + uncapped dogma (prod_m_max≤0) + lower α_tx + translation output scale.
SPN_PIPELINE_VERSION_V13 = "spn_v13_uncapped_dogma"
# v14: v12 + uncapped dogma only (α_tx / translation scale unchanged).
SPN_PIPELINE_VERSION_V14 = "spn_v14_uncapped_only"
# v15: v12 + promoter incoming edge ×scale + tuned dogma (α_tx / translation).
SPN_PIPELINE_VERSION_V15 = "spn_v15_promoter_scaled_dogma"
# v16: v15 stack + global protein_production_output_scale on PRODUCE_TRANSLATION.
SPN_PIPELINE_VERSION_V16 = "spn_v16_protein_prod_scaled_dogma"
# std: validated H1 B0 (fast×e + e3ewe + uncapped + met-source maintain + o150 pulse).
SPN_PIPELINE_VERSION_STD = "spn_std_h1_b0"


def resolve_v2_pipeline_version(
    *,
    prod_m_max: float,
    prod_alpha_transcription: float,
    translation_output_scale: float,
    promoter_incoming_edge_scale: float = 1.0,
    protein_production_output_scale: float = 1.0,
    protein_half_life_blocks: float | None = None,
    strong_e3_inhibition_weight: float = 1.0,
    weak_e3_edge_scale: float = WEAK_E3_EDGE_SCALE_BASELINE,
) -> str:
    """Pick v12–std label from dogma / promoter knobs (layer_rules=v2 + TF→DNA edge)."""
    from gin_tool.subgraph_pipeline.spn_standard_policy import is_standard_h1_b0_knobs

    if protein_half_life_blocks is not None and is_standard_h1_b0_knobs(
        prod_m_max=prod_m_max,
        prod_alpha_transcription=prod_alpha_transcription,
        promoter_incoming_edge_scale=promoter_incoming_edge_scale,
        strong_e3_inhibition_weight=strong_e3_inhibition_weight,
        weak_e3_edge_scale=weak_e3_edge_scale,
        protein_half_life_blocks=protein_half_life_blocks,
    ):
        return SPN_PIPELINE_VERSION_STD
    if float(protein_production_output_scale) != 1.0:
        return SPN_PIPELINE_VERSION_V16
    if float(promoter_incoming_edge_scale) != 1.0:
        return SPN_PIPELINE_VERSION_V15
    if float(prod_m_max) > 0.0:
        return SPN_PIPELINE_VERSION_V12
    if float(translation_output_scale) != 1.0 or float(prod_alpha_transcription) <= 0.01:
        return SPN_PIPELINE_VERSION_V13
    return SPN_PIPELINE_VERSION_V14


@dataclass(frozen=True)
class SpnSimPolicy:
    version: str = SPN_PIPELINE_VERSION
    blocks: int = 200
    runs: int = 200
    maintain_seeds: bool = True
    maintain_every_n_blocks: int = 1
    """Re-inject / clamp seed places every N blocks (1=every block; 2=half rate)."""
    maintain_metabolite_sources: bool = True
    """When maintain_seeds=False, still clamp met↔met in-degree-0 source places each block."""
    add_degradation: bool = True

    add_self_inhibition: bool = False
    rna_half_life_blocks: float = 5.0
    mirna_half_life_blocks: float | None = None
    """miRNA degradation half-life; None = same as rna_half_life_blocks (required outlet: repression does not consume miRNA)."""
    protein_half_life_blocks: float = 20.0
    high_tokens: int = 10
    seed: int = 42
    n_jobs: int = 20
    """Parallel workers for independent MC runs (multiprocess pool)."""
    mark_mode: str = "ancestor_all"
    """Token injection: ancestor_all (=in-degree 0) | all_places | ancestor_protein | source_protein | seed_protein."""
    layer_rules: str = LAYER_RULES_OFF
    """off = legacy gin_consume_v1; v2 = layered semantics."""
    stoich: str = "off"
    """off | one_to_one | kegg — metabolic stoichiometry mode (kegg wired in phase C)."""
    prod_alpha_transcription: float = 1.0
    prod_alpha_translation: float = 1.0
    prod_m_max: float = 50.0
    """Transcription/translation sample ceiling; ≤0 = uncapped (use full level)."""
    translation_output_scale: float = 1.0
    """Legacy alias; prefer ``protein_production_output_scale`` (product applied at translation)."""
    protein_production_output_scale: float = 1.0
    """Global scale on activating PRODUCE_TRANSLATION mint (all RNA→protein edges; v2 rule)."""
    tf_proteolysis_eps: float = 0.0
    histone_half_life_blocks: float = 30.0
    promoter_active_half_life_blocks: float = 20.0
    metabolite_half_life_blocks: float | None = None
    """If set, degrade ALL metabolite places. Default None = only DEFAULT_EXTRA whitelist."""
    activity_sigmoid_m_max: float = 20.0
    """Capacity for promoter.active / promoter;state / histone marks under sigmoid updates."""
    activity_sigmoid_k: float = 1.0
    """Steepness of logit-drive and transcription readout sigmoid."""
    activity_update_mode: str = "sigmoid"
    """REGULATORY / activity write: ``sigmoid`` | ``relu`` | ``relu_active``."""
    activity_relu_drive_scale: float = 1.0
    """Scale factor α for promoter.active under ReLU (see activity_relu_scale_mode)."""
    activity_relu_scale_mode: str = "drive"
    """``drive``: multiply REGULATORY δ by α only (degradation unscaled).
    ``net``: leave δ unscaled; after each block apply ``m ← m0 + α·(m−m0)`` on
    promoter.active/;state (scales net input−output together)."""
    out_route_half_life_blocks: float = 20.0
    """First-order out-route turnover (blocks); δ_tot=f(m) then categorical edge pick. ≤0 disables routing."""
    hme_regulatory_budget_scale: float = 1.0
    """Multiply REGULATORY drive δ when source is an HME.protein (enzyme not consumed)."""
    chromhmm_emission_tsv: str | None = None
    """Optional chromHMM emission TSV → signed mark weights (enables spn_v8 when set with v2)."""
    chromhmm_mark_weights_tsv: str | None = None
    """Optional precomputed mark\\tweight TSV (overrides emission derivation)."""
    drop_protein_self_loops: bool = False
    """Remove homodimer R-V-R protein self-loops before Petri export."""
    fire_propensity_mode: str = "off"
    """off | geom_mean — Bernoulli fire skip on CATALYST/REGULATORY (v2)."""
    fire_propensity_scale: float = 1.0
    fire_propensity_ref: float = 0.0
    """Reference level for geom-mean propensity; 0 → high_tokens at simulate time."""
    interaction_weight_pulse: dict | None = None
    """Optional MDM2⊣TP53-style interaction weight pulse schedule."""
    regulatory_norm: str = "off"
    """off | target_softmax | hill_target_softmax."""
    regulatory_hill_k: float = 5.0
    regulatory_target_budget: float = 10.0
    ups_degradation_policy: bool = True
    """Scale weak E3 UPS edges; skip global protein degradation for strong ESI substrates."""
    ups_e3_members_tsv: str | None = None
    ups_esi_evidence_tsv: str | None = None
    ups_e3_pair_override_tsv: str | None = None
    weak_e3_edge_scale: float = WEAK_E3_EDGE_SCALE_BASELINE
    strong_e3_inhibition_weight: float = 1.0
    """Effective CATALYST inhibition w for strong ESI pairs (|w_in×w_out|; default 1.0)."""
    scale_non_e3_activity_inhibitors: bool = True
    """Gene-list rule: non-E3 activity-repressor enzymes scale inhibitory protein drains."""
    non_e3_activity_inhibitors_tsv: str | None = None
    non_e3_activity_inhibitor_scale: float = 0.1
    drop_hme_to_promoter_state: bool = True
    """Drop HME.protein→*.promoter;state; HME acts only via histone marks (nodeList_HME)."""
    histone_to_promoter_state_scale: float = 0.1
    """Scale histone_mark→*.promoter;state (1.0 = unchanged)."""
    promoter_incoming_edge_scale: float = 1.0
    """Multiply all GIN/Petri weights on edges targeting *.promoter;state / *.promoter.active."""
    sign_repressive_chromatin_to_state: bool = True
    """DNA-methylation-axis repressors → ;state weight sign from seed list (typically −1)."""
    repressive_chromatin_enzymes_tsv: str | None = None
    tf_dna_route_independent: bool = True
    """TF.protein→promoter hub: independent REGULATORY fire, exempt from out-route fan-out."""
    tf_dna_edge_weight: float = 0.5
    """Fixed multiplicative edge weight w for TF→DNA REGULATORY (δ ~ sample(m_TF × w))."""


DEFAULT_SPN_POLICY = SpnSimPolicy()


def policy_for_layer_rules(
    layer_rules: str,
    base: SpnSimPolicy | None = None,
    *,
    chromhmm_emission_tsv: str | None = None,
    chromhmm_mark_weights_tsv: str | None = None,
) -> SpnSimPolicy:
    """Return a policy with version/layer_rules set; v2 enables kind-wise degradation."""
    b = base or DEFAULT_SPN_POLICY
    em = (
        chromhmm_emission_tsv
        if chromhmm_emission_tsv is not None
        else getattr(b, "chromhmm_emission_tsv", None)
    )
    wt = (
        chromhmm_mark_weights_tsv
        if chromhmm_mark_weights_tsv is not None
        else getattr(b, "chromhmm_mark_weights_tsv", None)
    )
    kwargs = dict(
        blocks=b.blocks,
        runs=b.runs,
        maintain_seeds=b.maintain_seeds,
        maintain_every_n_blocks=getattr(b, "maintain_every_n_blocks", 1),
        maintain_metabolite_sources=getattr(b, "maintain_metabolite_sources", True),
        add_degradation=b.add_degradation,
        add_self_inhibition=b.add_self_inhibition,
        rna_half_life_blocks=b.rna_half_life_blocks,
        mirna_half_life_blocks=getattr(b, "mirna_half_life_blocks", None),
        protein_half_life_blocks=b.protein_half_life_blocks,
        high_tokens=b.high_tokens,
        seed=b.seed,
        n_jobs=b.n_jobs,
        mark_mode=b.mark_mode,
        prod_alpha_transcription=b.prod_alpha_transcription,
        prod_alpha_translation=b.prod_alpha_translation,
        prod_m_max=b.prod_m_max,
        translation_output_scale=float(getattr(b, "translation_output_scale", 1.0)),
        protein_production_output_scale=float(
            getattr(b, "protein_production_output_scale", 1.0)
        ),
        tf_proteolysis_eps=b.tf_proteolysis_eps,
        histone_half_life_blocks=b.histone_half_life_blocks,
        promoter_active_half_life_blocks=b.promoter_active_half_life_blocks,
        metabolite_half_life_blocks=b.metabolite_half_life_blocks,
        activity_sigmoid_m_max=b.activity_sigmoid_m_max,
        activity_sigmoid_k=b.activity_sigmoid_k,
        activity_update_mode=getattr(b, "activity_update_mode", "sigmoid"),
        activity_relu_drive_scale=getattr(b, "activity_relu_drive_scale", 1.0),
        activity_relu_scale_mode=getattr(b, "activity_relu_scale_mode", "drive"),
        out_route_half_life_blocks=getattr(b, "out_route_half_life_blocks", 20.0),
        interaction_weight_pulse=getattr(b, "interaction_weight_pulse", None),
        drop_protein_self_loops=getattr(b, "drop_protein_self_loops", False),
        fire_propensity_mode=getattr(b, "fire_propensity_mode", "off"),
        fire_propensity_scale=getattr(b, "fire_propensity_scale", 1.0),
        fire_propensity_ref=getattr(b, "fire_propensity_ref", 0.0),
        hme_regulatory_budget_scale=getattr(b, "hme_regulatory_budget_scale", 1.0),
        regulatory_norm=getattr(b, "regulatory_norm", "off"),
        regulatory_hill_k=getattr(b, "regulatory_hill_k", 5.0),
        regulatory_target_budget=getattr(b, "regulatory_target_budget", 10.0),
        ups_degradation_policy=getattr(b, "ups_degradation_policy", True),
        ups_e3_members_tsv=getattr(b, "ups_e3_members_tsv", None),
        ups_esi_evidence_tsv=getattr(b, "ups_esi_evidence_tsv", None),
        ups_e3_pair_override_tsv=getattr(b, "ups_e3_pair_override_tsv", None),
        weak_e3_edge_scale=getattr(
            b, "weak_e3_edge_scale", WEAK_E3_EDGE_SCALE_BASELINE
        ),
        strong_e3_inhibition_weight=float(
            getattr(b, "strong_e3_inhibition_weight", 1.0)
        ),
        scale_non_e3_activity_inhibitors=getattr(
            b, "scale_non_e3_activity_inhibitors", True
        ),
        non_e3_activity_inhibitors_tsv=getattr(
            b, "non_e3_activity_inhibitors_tsv", None
        ),
        non_e3_activity_inhibitor_scale=getattr(
            b, "non_e3_activity_inhibitor_scale", 0.1
        ),
        drop_hme_to_promoter_state=getattr(b, "drop_hme_to_promoter_state", True),
        histone_to_promoter_state_scale=getattr(
            b, "histone_to_promoter_state_scale", 0.1
        ),
        promoter_incoming_edge_scale=float(
            getattr(b, "promoter_incoming_edge_scale", 1.0)
        ),
        sign_repressive_chromatin_to_state=getattr(
            b, "sign_repressive_chromatin_to_state", True
        ),
        repressive_chromatin_enzymes_tsv=getattr(
            b, "repressive_chromatin_enzymes_tsv", None
        ),
        tf_dna_route_independent=getattr(b, "tf_dna_route_independent", True),
        tf_dna_edge_weight=float(
            getattr(b, "tf_dna_edge_weight", getattr(b, "tf_dna_fixed_drive", 0.5))
        ),
        chromhmm_emission_tsv=em,
        chromhmm_mark_weights_tsv=wt,
    )
    if layer_rules == LAYER_RULES_V2:
        ver = resolve_v2_pipeline_version(
            prod_m_max=float(kwargs["prod_m_max"]),
            prod_alpha_transcription=float(kwargs["prod_alpha_transcription"]),
            translation_output_scale=float(kwargs["translation_output_scale"]),
            promoter_incoming_edge_scale=float(kwargs["promoter_incoming_edge_scale"]),
            protein_production_output_scale=float(
                kwargs["protein_production_output_scale"]
            ),
            protein_half_life_blocks=float(kwargs["protein_half_life_blocks"]),
            strong_e3_inhibition_weight=float(kwargs["strong_e3_inhibition_weight"]),
            weak_e3_edge_scale=float(kwargs["weak_e3_edge_scale"]),
        )
        return SpnSimPolicy(
            version=ver,
            layer_rules=LAYER_RULES_V2,
            stoich=b.stoich if b.stoich != "off" else "one_to_one",
            **kwargs,
        )
    return SpnSimPolicy(
        version=b.version,
        layer_rules=LAYER_RULES_OFF,
        stoich=b.stoich,
        **kwargs,
    )

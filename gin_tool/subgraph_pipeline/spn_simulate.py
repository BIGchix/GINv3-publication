"""Export bipartite GIN subgraph to mEPN and run SPN simulation (scBionet engine)."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Sequence

from gin_tool.mapping import is_virtual_gin_node
from gin_tool.spn_export import (
    SpnPetriNet,
    SpnTransitionKind,
    allow_petri_arc,
    apply_place_tokens,
    build_petri_net,
    classify_virtual_transition,
    classify_place,
    write_mepn_graphml,
    is_mirna_place,
)
from gin_tool.subgraph_pipeline.spn_degradation import (
    SINK_CHROMATIN,
    SINK_HISTONE,
    SINK_METABOLITE,
    SINK_PROTEIN,
    SINK_RNA,
    add_molecule_degradation_layer,
)
from gin_tool.subgraph_pipeline.spn_layer_rules import LAYER_RULES_V2
from gin_tool.subgraph_pipeline.spn_policy import (
    DEFAULT_SPN_POLICY,
    SPN_PIPELINE_VERSION_V2,
    SPN_PIPELINE_VERSION_V8,
    SPN_PIPELINE_VERSION_V9,
    SPN_PIPELINE_VERSION_V10,
    SPN_PIPELINE_VERSION_V11,
    WEAK_E3_EDGE_SCALE_BASELINE,
    resolve_v2_pipeline_version,
    SpnSimPolicy,
)
from gin_tool.subgraph_pipeline.spn_simulator_gin import GINSPNSimulator, GIN_FIRING_POLICY
from gin_tool.subgraph_pipeline.edge_weight_balance import (
    balance_gin_edge_weights,
    parse_edge_balance_layers,
)
from gin_tool.spn_rehydrate import RehydratedEdge, RehydratedSubgraph


def load_gin_edges_tsv(path: str | Path) -> list[tuple[str, str, float]]:
    rows: list[tuple[str, str, float]] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        a, b, w = line.split("\t", 2)
        rows.append((a, b, float(w)))
    return rows


def build_petri_net_from_bipartite_edges(
    gin_edges: list[tuple[str, str, float]],
    *,
    seed_gin_names: set[str] | None = None,
    mark_mode: str = "seed_protein",
) -> SpnPetriNet:
    """Build place/transition net from full GIN bipartite edge list.

    For pure R↔V edge lists (the context export), arcs are materialized directly
    from ``gin_edges`` in O(E). Logical R→V→R expansion is skipped — it duplicated
    the same arcs and previously hit O(E²) via ``list.__contains__`` on ``net.arcs``.
    Place–place edges (rare) still get a synthetic transition.

    Token marking runs **after** arcs exist so ``source_protein`` can use in-degree.
    """
    real_nodes: set[str] = set()
    virtual_nodes: set[str] = set()
    place_place: list[tuple[str, str, float]] = []
    for a, b, w in gin_edges:
        a_virt = is_virtual_gin_node(a)
        b_virt = is_virtual_gin_node(b)
        if a_virt:
            virtual_nodes.add(a)
        else:
            real_nodes.add(a)
        if b_virt:
            virtual_nodes.add(b)
        else:
            real_nodes.add(b)
        if not a_virt and not b_virt:
            place_place.append((a, b, w))

    logical: list[RehydratedEdge] = [
        RehydratedEdge(source=a, target=b, sign=w, via_virtual=None)
        for a, b, w in place_place
    ]
    rehydrated = RehydratedSubgraph(
        real_nodes=real_nodes,
        virtual_nodes=virtual_nodes,
        gin_edges=list(gin_edges),
        logical_edges=logical,
        rehydrated_count=0,
        direct_count=len(logical),
    )

    # Defer marking until arcs are known (source_protein needs in-degree).
    net = build_petri_net(rehydrated, seed_gin_names=set())

    # Direct arcs from gin_edges (place→transition or transition→place).
    # Preserve signed weights (e.g. repressive histone → promoter;state = -1).
    arc_set = set(net.arcs)
    for a, b, w in gin_edges:
        if is_virtual_gin_node(a):
            if a not in net.transitions:
                net.transitions[a] = classify_virtual_transition(a)
        else:
            if a not in net.places:
                net.places[a] = classify_place(a)
        if is_virtual_gin_node(b):
            if b not in net.transitions:
                net.transitions[b] = classify_virtual_transition(b)
        else:
            if b not in net.places:
                net.places[b] = classify_place(b)
        if not allow_petri_arc(a, b):
            continue
        wf = float(w)
        if (a, b) not in arc_set:
            arc_set.add((a, b))
            net.arcs.append((a, b))
            net.arc_weights[(a, b)] = wf
        else:
            prev = net.arc_weights.get((a, b), 1.0)
            # Prefer an inhibitory sign if either copy is negative.
            if wf < 0 or prev >= 0:
                net.arc_weights[(a, b)] = wf

    apply_place_tokens(net, seed_gin_names=seed_gin_names or set(), mark_mode=mark_mode)
    return net


def _candidate_scbionet_roots() -> list[Path]:
    """Locate vendored SigPN engine (no torch import via scbionet.__init__)."""
    import os

    repo = Path(__file__).resolve().parents[2]
    roots: list[Path] = []
    env = os.environ.get("GIN_SCBIONET_ROOT")
    if env:
        roots.append(Path(env))
    roots.extend(
        [
            repo / "scbionet",
            repo.parent / "xda" / "scBionet" / "scbionet",
            repo.parent / "xda" / "scBionet",
        ]
    )
    return roots


def _load_scb_spn_modules():
    """Import scBionet SPN modules without pulling torch via scbionet.__init__."""
    import importlib.util
    import types

    model_path = None
    sim_path = None
    tried: list[str] = []
    for root in _candidate_scbionet_roots():
        for base in (root, root / "scbionet"):
            m = base / "spn" / "model.py"
            s = base / "spn" / "simulator.py"
            tried.append(str(m))
            if m.is_file() and s.is_file():
                model_path, sim_path = m, s
                break
        if model_path is not None:
            break
    if model_path is None or sim_path is None:
        raise FileNotFoundError(
            "missing scBionet SPN modules; set GIN_SCBIONET_ROOT or keep "
            f"scbionet/spn/ in the repo. tried={tried[:8]}"
        )

    spec_m = importlib.util.spec_from_file_location("scb_spn_model", model_path)
    mod_m = importlib.util.module_from_spec(spec_m)
    assert spec_m.loader is not None
    sys.modules["scb_spn_model"] = mod_m
    spec_m.loader.exec_module(mod_m)

    sim_code = sim_path.read_text(encoding="utf-8").replace(
        "from scbionet.spn.model import SPNInteraction, SPNModel",
        "from scb_spn_model import SPNInteraction, SPNModel",
    )
    mod_s = types.ModuleType("scb_spn_simulator")
    sys.modules["scb_spn_simulator"] = mod_s
    exec(compile(sim_code, str(sim_path), "exec"), mod_s.__dict__)
    return mod_m, mod_s


def petri_net_to_spn_model(net: SpnPetriNet, *, add_self_inhibition: bool = False):
    """Convert SpnPetriNet to scBionet SPNModel."""
    mod_m, _mod_s = _load_scb_spn_modules()
    SPNInteraction = mod_m.SPNInteraction
    SPNModel = mod_m.SPNModel

    places = sorted(net.places.keys())
    place_set = set(places)
    transitions = set(net.transitions.keys())

    arc_w = getattr(net, "arc_weights", None) or {}
    trans_in: dict[str, list[tuple[str, float]]] = {t: [] for t in transitions}
    trans_out: dict[str, list[tuple[str, float]]] = {t: [] for t in transitions}
    for src, tgt in net.arcs:
        w = float(arc_w.get((src, tgt), 1.0))
        if src in place_set and tgt in transitions:
            trans_in[tgt].append((src, w))
        elif src in transitions and tgt in place_set:
            trans_out[src].append((tgt, w))

    interactions: list[SPNInteraction] = []
    degradation_weights: dict[str, float] = getattr(net, "degradation_weights", {}) or {}
    sink_places = {SINK_RNA, SINK_PROTEIN, SINK_METABOLITE, SINK_HISTONE, SINK_CHROMATIN}

    for tid in transitions:
        if net.transitions.get(tid) == SpnTransitionKind.DEGRADATION:
            ins = trans_in.get(tid, [])
            if len(ins) == 1:
                src, _w_in = ins[0]
                w = float(degradation_weights.get(src, _w_in))
                interactions.append(
                    SPNInteraction(
                        source=src,
                        target=src,
                        kind="inhibition",
                        transition_id=tid,
                        weight=w,
                    )
                )
            continue

        t_kind = net.transitions.get(tid)
        for src, w_in in trans_in.get(tid, []):
            for tgt, w_out in trans_out.get(tid, []):
                if tgt in sink_places:
                    continue
                # Transcription produces RNA only (never promoter/DNA self-feed).
                if t_kind == SpnTransitionKind.TRANSCRIPTION and not allow_petri_arc(tid, tgt):
                    continue
                # GIN signs: histone often (-1,+1); protein inhibition often (-1,-1).
                # Treat any negative endpoint weight as inhibition (not only XOR/product<0).
                kind = "inhibition" if (w_in < 0 or w_out < 0) else "activation"
                # miRNA → translation is biologically repressive; GIN often stores +1.
                tid_l = str(tid).lower()
                if kind == "activation" and "translation" in tid_l and is_mirna_place(src):
                    kind = "inhibition"
                interactions.append(
                    SPNInteraction(
                        source=src,
                        target=tgt,
                        kind=kind,
                        transition_id=tid,
                        weight=abs(w_in * w_out) if (w_in * w_out) != 0 else abs(w_in or w_out),
                    )
                )

    model = SPNModel(places=places, interactions=interactions, place_labels={p: p for p in places})
    for p, tok in net.place_tokens.items():
        if tok and tok > 0:
            model.fixed_high.add(p)
    if add_self_inhibition:
        model.add_self_inhibition_loops()
    return model


def seed_gin_names_from_manifest(manifest_path: Path) -> set[str]:
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    if data.get("gin_anchor_nodes"):
        return set(data["gin_anchor_nodes"])
    genes = data.get("union_genes") or []
    kinds = data.get("seed_kinds") or {}
    names: set[str] = set()
    for g in genes:
        if kinds.get(g) == "metabolite":
            names.add(g)
            continue
        names.add(g)
        names.add(f"{g}.RNA")
        names.add(f"{g}.protein")
        names.add(f"{g}.promoter.active")
    return names


def apply_place_token_levels(net: SpnPetriNet, levels: dict[str, float]) -> set[str]:
    """Mark only places present in ``levels`` (value > 0); clear others."""
    marked: set[str] = set()
    for p in net.places:
        val = float(levels.get(p, 0.0) or 0.0)
        if val > 0.0:
            net.place_tokens[p] = 1
            marked.add(p)
        else:
            net.place_tokens[p] = 0
    return marked


def load_place_levels_tsv(path: str | Path) -> dict[str, float]:
    """Load ``place\\tlevel`` TSV (header optional)."""
    out: dict[str, float] = {}
    for i, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines()):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        place, raw = parts[0].strip(), parts[1].strip()
        if i == 0 and place.lower() in {"place", "name", "node"}:
            continue
        try:
            out[place] = float(raw)
        except ValueError:
            continue
    return out


def run_spn_simulation(
    model,
    *,
    out_dir: Path,
    stem: str,
    policy: SpnSimPolicy | None = None,
    blocks: int | None = None,
    runs: int | None = None,
    seed: int | None = None,
    maintain_seeds: bool | None = None,
    maintain_every_n_blocks: int | None = None,
    maintain_places: frozenset[str] | None = None,
    maintain_levels: dict[str, float] | None = None,
    keep_trajectories: bool = True,
    n_jobs: int | None = None,
    chromhmm_mark_weights: dict[str, float] | None = None,
    chromhmm_weight_meta: dict | None = None,
    interaction_weight_pulse: dict | None = None,
    nf_qc_pair: tuple[str, str] | None = None,
    nf_warm: int = 20,
    nf_max_lag: int = 40,
    regulatory_norm: str | None = None,
    regulatory_hill_k: float | None = None,
    regulatory_target_budget: float | None = None,
) -> dict:
    pol = policy or DEFAULT_SPN_POLICY
    blocks = pol.blocks if blocks is None else blocks
    runs = pol.runs if runs is None else runs
    seed = pol.seed if seed is None else seed
    maintain_seeds = pol.maintain_seeds if maintain_seeds is None else maintain_seeds
    maintain_every_n_blocks = (
        getattr(pol, "maintain_every_n_blocks", 1)
        if maintain_every_n_blocks is None
        else int(maintain_every_n_blocks)
    )
    n_jobs = pol.n_jobs if n_jobs is None else n_jobs
    out_dir.mkdir(parents=True, exist_ok=True)
    pulse = interaction_weight_pulse
    if pulse is None:
        pulse = getattr(pol, "interaction_weight_pulse", None)
    reg_norm = (
        str(getattr(pol, "regulatory_norm", "off") or "off")
        if regulatory_norm is None
        else str(regulatory_norm)
    )
    reg_k = (
        float(getattr(pol, "regulatory_hill_k", 5.0))
        if regulatory_hill_k is None
        else float(regulatory_hill_k)
    )
    reg_bud = (
        float(getattr(pol, "regulatory_target_budget", 10.0))
        if regulatory_target_budget is None
        else float(regulatory_target_budget)
    )
    sim = GINSPNSimulator(
        n_blocks=blocks,
        n_runs=runs,
        seed=seed,
        high_tokens=pol.high_tokens,
        maintain_seeds=maintain_seeds,
        maintain_every_n_blocks=maintain_every_n_blocks,
        maintain_places=maintain_places,
        maintain_levels=maintain_levels,
        n_jobs=n_jobs,
        layer_rules=getattr(pol, "layer_rules", "off"),
        prod_alpha_transcription=getattr(pol, "prod_alpha_transcription", 1.0),
        prod_alpha_translation=getattr(pol, "prod_alpha_translation", 1.0),
        prod_m_max=getattr(pol, "prod_m_max", 50.0),
        translation_output_scale=float(getattr(pol, "translation_output_scale", 1.0)),
        protein_production_output_scale=float(
            getattr(pol, "protein_production_output_scale", 1.0)
        ),
        tf_proteolysis_eps=getattr(pol, "tf_proteolysis_eps", 0.0),
        activity_sigmoid_m_max=getattr(pol, "activity_sigmoid_m_max", 20.0),
        activity_sigmoid_k=getattr(pol, "activity_sigmoid_k", 1.0),
        activity_update_mode=getattr(pol, "activity_update_mode", "sigmoid"),
        activity_relu_drive_scale=getattr(pol, "activity_relu_drive_scale", 1.0),
        activity_relu_scale_mode=getattr(pol, "activity_relu_scale_mode", "drive"),
        out_route_half_life_blocks=getattr(pol, "out_route_half_life_blocks", 20.0),
        hme_regulatory_budget_scale=getattr(pol, "hme_regulatory_budget_scale", 1.0),
        chromhmm_mark_weights=chromhmm_mark_weights,
        fire_propensity_mode=str(getattr(pol, "fire_propensity_mode", "off") or "off"),
        fire_propensity_scale=float(getattr(pol, "fire_propensity_scale", 1.0)),
        fire_propensity_ref=float(
            getattr(pol, "fire_propensity_ref", 0.0) or float(pol.high_tokens)
        ),
        interaction_weight_pulse=pulse,
        nf_qc_pair=nf_qc_pair,
        nf_warm=int(nf_warm),
        nf_max_lag=int(nf_max_lag),
        regulatory_norm=reg_norm,
        regulatory_hill_k=reg_k,
        regulatory_target_budget=reg_bud,
        tf_dna_route_independent=bool(getattr(pol, "tf_dna_route_independent", True)),
        tf_dna_edge_weight=float(
            getattr(pol, "tf_dna_edge_weight", getattr(pol, "tf_dna_fixed_drive", 0.5))
        ),
    )
    result = sim.run(model)
    csv_path = out_dir / f"{stem}_trajectory.csv"
    json_path = out_dir / f"{stem}_trajectory.json"
    if keep_trajectories:
        result.save_csv(csv_path)
        result.save_json(json_path)

    if result.nf_passrate is not None and result.nf_per_run is not None:
        nf_path = out_dir / f"{stem}_nf_passrate.json"
        nf_path.write_text(
            json.dumps(
                {
                    "pair": list(nf_qc_pair) if nf_qc_pair else None,
                    "warm": int(nf_warm),
                    "max_lag": int(nf_max_lag),
                    "passrate": result.nf_passrate,
                    "per_run": result.nf_per_run,
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )
        if result.pair_block_samples is not None and result.pair_place_names is not None:
            import numpy as np

            npz_path = out_dir / f"{stem}_pair_runs.npz"
            a_name, b_name = result.pair_place_names
            np.savez_compressed(
                npz_path,
                places=np.array([a_name, b_name]),
                block_samples=result.pair_block_samples,
            )

    # top variable places at final block
    final = result.mean_markings[-1]
    order = sorted(range(len(final)), key=lambda i: -final[i])
    top = [(result.place_names[i], float(final[i])) for i in order[:15]]

    summary = {
        "spn_pipeline_version": pol.version,
        "n_places": len(model.places),
        "n_interactions": len(model.interactions),
        "n_marked_seeds": len(model.fixed_high),
        "marked_seed_names": sorted(getattr(model, "fixed_high", set()) or []),
        "firing_policy": result.firing_policy or GIN_FIRING_POLICY,
        "layer_rules": getattr(pol, "layer_rules", "off"),
        "stoich": getattr(pol, "stoich", "off"),
        "activity_sigmoid_m_max": getattr(pol, "activity_sigmoid_m_max", 20.0),
        "activity_sigmoid_k": getattr(pol, "activity_sigmoid_k", 1.0),
        "activity_update_mode": getattr(pol, "activity_update_mode", "sigmoid"),
        "activity_relu_drive_scale": getattr(pol, "activity_relu_drive_scale", 1.0),
        "activity_relu_scale_mode": getattr(pol, "activity_relu_scale_mode", "drive"),
        "fire_propensity_mode": getattr(pol, "fire_propensity_mode", "off"),
        "fire_propensity_scale": getattr(pol, "fire_propensity_scale", 1.0),
        "fire_propensity_ref": float(
            getattr(pol, "fire_propensity_ref", 0.0) or pol.high_tokens
        ),
        "promoter_active_half_life_blocks": getattr(
            pol, "promoter_active_half_life_blocks", 20.0
        ),
        "hme_regulatory_budget_scale": getattr(pol, "hme_regulatory_budget_scale", 1.0),
        "chromhmm_emw": chromhmm_weight_meta,
        "chromhmm_n_marks": len(chromhmm_mark_weights or {}),
        "ledger_totals": result.ledger_totals,
        "maintain_seeds": maintain_seeds,
        "maintain_every_n_blocks": maintain_every_n_blocks,
        "n_maintain_levels": len(maintain_levels or {}),
        "prod_alpha_transcription": getattr(pol, "prod_alpha_transcription", 1.0),
        "prod_alpha_translation": getattr(pol, "prod_alpha_translation", 1.0),
        "blocks": blocks,
        "runs": runs,
        "n_jobs": n_jobs,
        "keep_trajectories": keep_trajectories,
        "trajectory_csv": str(csv_path) if keep_trajectories else None,
        "top_places_final_block": top,
        "regulatory_norm": reg_norm,
        "regulatory_hill_k": reg_k,
        "regulatory_target_budget": reg_bud,
    }
    if result.nf_passrate is not None:
        summary["nf_passrate"] = result.nf_passrate
        summary["nf_passrate_json"] = str(out_dir / f"{stem}_nf_passrate.json")
        summary["nf_qc_pair"] = list(nf_qc_pair) if nf_qc_pair else None
        if result.pair_block_samples is not None:
            summary["pair_runs_npz"] = str(out_dir / f"{stem}_pair_runs.npz")
    (out_dir / f"{stem}_sim_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def _run_reasonableness_qc_after_sim(
    *,
    sim_dir: Path,
    stem: str,
    subgraph_dir: Path,
    model,
    summary: dict,
) -> dict:
    from gin_tool.subgraph_pipeline.spn_reasonableness_qc import (
        emit_spn_reasonableness_qc,
        qc_spn_reasonableness_from_sim_dir,
    )

    # Persist mark_mode / degradation before QC so disk summary matches in-memory.
    (sim_dir / f"{stem}_sim_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    try:
        qc = qc_spn_reasonableness_from_sim_dir(
            sim_dir,
            stem,
            subgraph_dir=subgraph_dir,
            marked_seed_names=sorted(getattr(model, "fixed_high", set()) or []),
        )
    except FileNotFoundError as exc:
        print(f"[SPN-QC] skip reasonableness (missing traj): {exc}", flush=True)
        return summary
    out_json = sim_dir / f"{stem}_reasonableness_qc.json"
    emit_spn_reasonableness_qc(qc, out_json=out_json, print_summary=True)
    summary["reasonableness_qc"] = {
        "ok": qc.ok,
        "overall": qc.report.get("verdict", {}).get("overall"),
        "n_warnings": len(qc.warnings),
        "n_failures": len(qc.failures),
        "path": str(out_json),
    }
    (sim_dir / f"{stem}_sim_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def export_and_simulate_subgraph_dir(
    subgraph_dir: str | Path,
    *,
    policy: SpnSimPolicy | None = None,
    blocks: int | None = None,
    runs: int | None = None,
    add_degradation: bool | None = None,
    rna_half_life_blocks: float | None = None,
    mirna_half_life_blocks: float | None = None,
    protein_half_life_blocks: float | None = None,
    add_self_inhibition: bool | None = None,
    maintain_seeds: bool | None = None,
    maintain_every_n_blocks: int | None = None,
    maintain_metabolite_sources: bool | None = None,
    mark_mode: str | None = None,
    place_levels: dict[str, float] | None = None,
    keep_trajectories: bool = True,
    sim_subdir: str = "spn_sim",
    stem: str = "p53",
    n_jobs: int | None = None,
    layer_rules: str | None = None,
    stoich: str | None = None,
    tf_proteolysis_eps: float | None = None,
    activity_sigmoid_m_max: float | None = None,
    activity_sigmoid_k: float | None = None,
    activity_update_mode: str | None = None,
    activity_relu_drive_scale: float | None = None,
    activity_relu_scale_mode: str | None = None,
    promoter_active_half_life_blocks: float | None = None,
    out_route_half_life_blocks: float | None = None,
    edge_balance: str | list[str] | None = None,
    prod_alpha_transcription: float | None = None,
    prod_alpha_translation: float | None = None,
    prod_m_max: float | None = None,
    translation_output_scale: float | None = None,
    protein_production_output_scale: float | None = None,
    hme_regulatory_budget_scale: float | None = None,
    chromhmm_emission_tsv: str | Path | None = None,
    chromhmm_mark_weights_tsv: str | Path | None = None,
    drop_protein_self_loops: bool | None = None,
    fire_propensity_mode: str | None = None,
    fire_propensity_scale: float | None = None,
    fire_propensity_ref: float | None = None,
    interaction_weight_pulse: dict | None = None,
    seed: int | None = None,
    nf_qc_pair: tuple[str, str] | None = None,
    nf_warm: int = 20,
    nf_max_lag: int = 40,
    regulatory_norm: str | None = None,
    regulatory_hill_k: float | None = None,
    regulatory_target_budget: float | None = None,
    drop_chromatin_enzyme_to_promoter_state: bool | None = None,
    histone_to_promoter_state_scale: float | None = None,
    promoter_incoming_edge_scale: float | None = None,
    sign_repressive_chromatin_to_state: bool | None = None,
    repressive_chromatin_enzymes_tsv: str | Path | None = None,
    tf_dna_route_independent: bool | None = None,
    tf_dna_edge_weight: float | None = None,
    tf_dna_fixed_drive: float | None = None,
    ups_degradation_policy: bool | None = None,
    weak_e3_edge_scale: float | None = None,
    strong_e3_inhibition_weight: float | None = None,
    ups_e3_members_tsv: str | Path | None = None,
    ups_esi_evidence_tsv: str | Path | None = None,
    ups_e3_pair_override_tsv: str | Path | None = None,
    drop_protein_inhibition_pairs: Sequence[tuple[str, str]] | bool | None = None,
    force_skip_protein_degradation: Sequence[str] | None = None,
    scale_non_e3_activity_inhibitors: bool | None = None,
    non_e3_activity_inhibitors_tsv: str | Path | None = None,
    non_e3_activity_inhibitor_scale: float | None = None,
) -> dict:
    pol = policy or DEFAULT_SPN_POLICY
    blocks = pol.blocks if blocks is None else blocks
    runs = pol.runs if runs is None else runs
    n_jobs = pol.n_jobs if n_jobs is None else n_jobs
    add_degradation = pol.add_degradation if add_degradation is None else add_degradation
    rna_half_life_blocks = pol.rna_half_life_blocks if rna_half_life_blocks is None else rna_half_life_blocks
    mirna_half_life_blocks = (
        getattr(pol, "mirna_half_life_blocks", None)
        if mirna_half_life_blocks is None
        else float(mirna_half_life_blocks)
    )
    protein_half_life_blocks = (
        pol.protein_half_life_blocks if protein_half_life_blocks is None else protein_half_life_blocks
    )
    add_self_inhibition = pol.add_self_inhibition if add_self_inhibition is None else add_self_inhibition
    maintain_seeds = pol.maintain_seeds if maintain_seeds is None else maintain_seeds
    maintain_every_n_blocks = (
        getattr(pol, "maintain_every_n_blocks", 1)
        if maintain_every_n_blocks is None
        else int(maintain_every_n_blocks)
    )
    maintain_metabolite_sources = (
        getattr(pol, "maintain_metabolite_sources", True)
        if maintain_metabolite_sources is None
        else bool(maintain_metabolite_sources)
    )
    mark_mode = pol.mark_mode if mark_mode is None else mark_mode
    layer_rules = getattr(pol, "layer_rules", "off") if layer_rules is None else layer_rules
    stoich = getattr(pol, "stoich", "off") if stoich is None else stoich
    tf_proteolysis_eps = (
        getattr(pol, "tf_proteolysis_eps", 0.0)
        if tf_proteolysis_eps is None
        else float(tf_proteolysis_eps)
    )
    activity_sigmoid_m_max = (
        getattr(pol, "activity_sigmoid_m_max", 20.0)
        if activity_sigmoid_m_max is None
        else float(activity_sigmoid_m_max)
    )
    activity_sigmoid_k = (
        getattr(pol, "activity_sigmoid_k", 1.0)
        if activity_sigmoid_k is None
        else float(activity_sigmoid_k)
    )
    activity_update_mode = (
        getattr(pol, "activity_update_mode", "sigmoid")
        if activity_update_mode is None
        else str(activity_update_mode)
    ).strip().lower() or "sigmoid"
    activity_relu_drive_scale = (
        getattr(pol, "activity_relu_drive_scale", 1.0)
        if activity_relu_drive_scale is None
        else float(activity_relu_drive_scale)
    )
    activity_relu_scale_mode = (
        str(getattr(pol, "activity_relu_scale_mode", "drive") or "drive")
        if activity_relu_scale_mode is None
        else str(activity_relu_scale_mode)
    ).strip().lower() or "drive"
    if activity_relu_scale_mode not in ("drive", "net"):
        raise ValueError(
            f"activity_relu_scale_mode must be 'drive' or 'net', got {activity_relu_scale_mode!r}"
        )
    promoter_active_half_life_blocks = (
        getattr(pol, "promoter_active_half_life_blocks", 20.0)
        if promoter_active_half_life_blocks is None
        else float(promoter_active_half_life_blocks)
    )
    out_route_half_life_blocks = (
        getattr(pol, "out_route_half_life_blocks", 20.0)
        if out_route_half_life_blocks is None
        else float(out_route_half_life_blocks)
    )
    prod_alpha_transcription = (
        getattr(pol, "prod_alpha_transcription", 1.0)
        if prod_alpha_transcription is None
        else float(prod_alpha_transcription)
    )
    prod_alpha_translation = (
        getattr(pol, "prod_alpha_translation", 1.0)
        if prod_alpha_translation is None
        else float(prod_alpha_translation)
    )
    prod_m_max_val = (
        float(getattr(pol, "prod_m_max", 50.0))
        if prod_m_max is None
        else float(prod_m_max)
    )
    translation_output_scale_val = (
        float(getattr(pol, "translation_output_scale", 1.0))
        if translation_output_scale is None
        else float(translation_output_scale)
    )
    protein_production_output_scale_val = (
        float(getattr(pol, "protein_production_output_scale", 1.0))
        if protein_production_output_scale is None
        else float(protein_production_output_scale)
    )
    promoter_incoming_edge_scale_val = (
        float(getattr(pol, "promoter_incoming_edge_scale", 1.0))
        if promoter_incoming_edge_scale is None
        else float(promoter_incoming_edge_scale)
    )
    weak_e3_edge_scale_val = (
        float(getattr(pol, "weak_e3_edge_scale", WEAK_E3_EDGE_SCALE_BASELINE))
        if weak_e3_edge_scale is None
        else float(weak_e3_edge_scale)
    )
    strong_e3_inhibition_weight_val = (
        float(getattr(pol, "strong_e3_inhibition_weight", 1.0))
        if strong_e3_inhibition_weight is None
        else float(strong_e3_inhibition_weight)
    )
    tf_dna_edge_weight_raw = (
        float(tf_dna_edge_weight)
        if tf_dna_edge_weight is not None
        else (
            float(tf_dna_fixed_drive)
            if tf_dna_fixed_drive is not None
            else float(
                getattr(
                    pol,
                    "tf_dna_edge_weight",
                    getattr(pol, "tf_dna_fixed_drive", 0.5),
                )
            )
        )
    )
    tf_dna_edge_weight_effective = tf_dna_edge_weight_raw * promoter_incoming_edge_scale_val
    hme_regulatory_budget_scale = (
        getattr(pol, "hme_regulatory_budget_scale", 1.0)
        if hme_regulatory_budget_scale is None
        else float(hme_regulatory_budget_scale)
    )
    chromhmm_emission_tsv = (
        getattr(pol, "chromhmm_emission_tsv", None)
        if chromhmm_emission_tsv is None
        else chromhmm_emission_tsv
    )
    chromhmm_mark_weights_tsv = (
        getattr(pol, "chromhmm_mark_weights_tsv", None)
        if chromhmm_mark_weights_tsv is None
        else chromhmm_mark_weights_tsv
    )
    drop_protein_self_loops = (
        bool(getattr(pol, "drop_protein_self_loops", False))
        if drop_protein_self_loops is None
        else bool(drop_protein_self_loops)
    )
    fire_propensity_mode = (
        str(getattr(pol, "fire_propensity_mode", "off") or "off")
        if fire_propensity_mode is None
        else str(fire_propensity_mode)
    )
    fire_propensity_scale = (
        float(getattr(pol, "fire_propensity_scale", 1.0))
        if fire_propensity_scale is None
        else float(fire_propensity_scale)
    )
    fire_propensity_ref = (
        float(getattr(pol, "fire_propensity_ref", 0.0))
        if fire_propensity_ref is None
        else float(fire_propensity_ref)
    )
    seed = pol.seed if seed is None else int(seed)
    if layer_rules == LAYER_RULES_V2 and stoich == "off":
        stoich = "one_to_one"
    use_v8 = layer_rules == LAYER_RULES_V2 and bool(
        chromhmm_emission_tsv or chromhmm_mark_weights_tsv
    )
    # Rebuild policy with CLI overrides so simulator sees layer_rules / alphas.
    pol = SpnSimPolicy(
        version=(
            resolve_v2_pipeline_version(
                prod_m_max=prod_m_max_val,
                prod_alpha_transcription=prod_alpha_transcription,
                translation_output_scale=translation_output_scale_val,
                promoter_incoming_edge_scale=promoter_incoming_edge_scale_val,
                protein_production_output_scale=protein_production_output_scale_val,
                protein_half_life_blocks=protein_half_life_blocks,
                strong_e3_inhibition_weight=strong_e3_inhibition_weight_val,
                weak_e3_edge_scale=weak_e3_edge_scale_val,
            )
            if layer_rules == LAYER_RULES_V2
            else pol.version
        ),
        blocks=pol.blocks,
        runs=pol.runs,
        maintain_seeds=maintain_seeds,
        maintain_every_n_blocks=maintain_every_n_blocks,
        maintain_metabolite_sources=maintain_metabolite_sources,
        add_degradation=pol.add_degradation,
        add_self_inhibition=pol.add_self_inhibition,
        rna_half_life_blocks=rna_half_life_blocks,
        mirna_half_life_blocks=mirna_half_life_blocks,
        protein_half_life_blocks=protein_half_life_blocks,
        high_tokens=pol.high_tokens,
        seed=seed,
        n_jobs=pol.n_jobs,
        mark_mode=pol.mark_mode,
        layer_rules=layer_rules,
        stoich=stoich,
        prod_alpha_transcription=prod_alpha_transcription,
        prod_alpha_translation=prod_alpha_translation,
        prod_m_max=prod_m_max_val,
        translation_output_scale=translation_output_scale_val,
        protein_production_output_scale=protein_production_output_scale_val,
        tf_proteolysis_eps=tf_proteolysis_eps,
        histone_half_life_blocks=getattr(pol, "histone_half_life_blocks", 30.0),
        promoter_active_half_life_blocks=promoter_active_half_life_blocks,
        metabolite_half_life_blocks=getattr(pol, "metabolite_half_life_blocks", None),
        activity_sigmoid_m_max=activity_sigmoid_m_max,
        activity_sigmoid_k=activity_sigmoid_k,
        activity_update_mode=activity_update_mode,
        activity_relu_drive_scale=activity_relu_drive_scale,
        activity_relu_scale_mode=activity_relu_scale_mode,
        out_route_half_life_blocks=out_route_half_life_blocks,
        hme_regulatory_budget_scale=hme_regulatory_budget_scale,
        chromhmm_emission_tsv=str(chromhmm_emission_tsv) if chromhmm_emission_tsv else None,
        chromhmm_mark_weights_tsv=(
            str(chromhmm_mark_weights_tsv) if chromhmm_mark_weights_tsv else None
        ),
        drop_protein_self_loops=drop_protein_self_loops,
        fire_propensity_mode=fire_propensity_mode,
        fire_propensity_scale=fire_propensity_scale,
        fire_propensity_ref=fire_propensity_ref,
        interaction_weight_pulse=(
            interaction_weight_pulse
            if interaction_weight_pulse is not None
            else getattr(pol, "interaction_weight_pulse", None)
        ),
        regulatory_norm=(
            regulatory_norm
            if regulatory_norm is not None
            else getattr(pol, "regulatory_norm", "off")
        ),
        regulatory_hill_k=(
            float(regulatory_hill_k)
            if regulatory_hill_k is not None
            else float(getattr(pol, "regulatory_hill_k", 5.0))
        ),
        regulatory_target_budget=(
            float(regulatory_target_budget)
            if regulatory_target_budget is not None
            else float(getattr(pol, "regulatory_target_budget", 10.0))
        ),
        ups_degradation_policy=getattr(pol, "ups_degradation_policy", True),
        ups_e3_members_tsv=getattr(pol, "ups_e3_members_tsv", None),
        ups_esi_evidence_tsv=getattr(pol, "ups_esi_evidence_tsv", None),
        ups_e3_pair_override_tsv=getattr(pol, "ups_e3_pair_override_tsv", None),
        weak_e3_edge_scale=getattr(pol, "weak_e3_edge_scale", 0.3),
        strong_e3_inhibition_weight=float(
            getattr(pol, "strong_e3_inhibition_weight", 1.0)
        ),
        scale_non_e3_activity_inhibitors=getattr(
            pol, "scale_non_e3_activity_inhibitors", True
        ),
        non_e3_activity_inhibitors_tsv=getattr(
            pol, "non_e3_activity_inhibitors_tsv", None
        ),
        non_e3_activity_inhibitor_scale=getattr(
            pol, "non_e3_activity_inhibitor_scale", 0.1
        ),
        drop_hme_to_promoter_state=(
            bool(drop_chromatin_enzyme_to_promoter_state)
            if drop_chromatin_enzyme_to_promoter_state is not None
            else bool(getattr(pol, "drop_hme_to_promoter_state", True))
        ),
        histone_to_promoter_state_scale=(
            float(histone_to_promoter_state_scale)
            if histone_to_promoter_state_scale is not None
            else float(getattr(pol, "histone_to_promoter_state_scale", 0.1))
        ),
        promoter_incoming_edge_scale=promoter_incoming_edge_scale_val,
        sign_repressive_chromatin_to_state=(
            bool(sign_repressive_chromatin_to_state)
            if sign_repressive_chromatin_to_state is not None
            else bool(getattr(pol, "sign_repressive_chromatin_to_state", True))
        ),
        repressive_chromatin_enzymes_tsv=(
            str(repressive_chromatin_enzymes_tsv)
            if repressive_chromatin_enzymes_tsv is not None
            else getattr(pol, "repressive_chromatin_enzymes_tsv", None)
        ),
        tf_dna_route_independent=(
            bool(tf_dna_route_independent)
            if tf_dna_route_independent is not None
            else bool(getattr(pol, "tf_dna_route_independent", True))
        ),
        tf_dna_edge_weight=tf_dna_edge_weight_effective,
    )
    # Resolve promoter-edge kwargs after policy rebuild.
    drop_hme = bool(getattr(pol, "drop_hme_to_promoter_state", True))
    hist_state_scale = float(getattr(pol, "histone_to_promoter_state_scale", 0.1))
    incoming_prom_scale = float(getattr(pol, "promoter_incoming_edge_scale", 1.0))
    sign_repr = bool(getattr(pol, "sign_repressive_chromatin_to_state", True))
    repr_tsv = repressive_chromatin_enzymes_tsv or getattr(
        pol, "repressive_chromatin_enzymes_tsv", None
    )
    subgraph_dir = Path(subgraph_dir)
    edges_path = subgraph_dir / "gin_edges.tsv"
    if not edges_path.is_file():
        alt = subgraph_dir / "subgraph_edges.tsv"
        if alt.is_file():
            edges_path = alt
        else:
            raise FileNotFoundError(edges_path)
    manifest_path = subgraph_dir / "seed_manifest.json"

    gin_edges = load_gin_edges_tsv(edges_path)
    promoter_edge_stats = None
    if drop_hme or hist_state_scale != 1.0 or incoming_prom_scale != 1.0 or sign_repr:
        from gin_tool.subgraph_pipeline.promoter_edge_policy import apply_promoter_edge_policy

        gin_edges, promoter_edge_stats = apply_promoter_edge_policy(
            gin_edges,
            drop_chromatin_enzyme_to_state=drop_hme,
            histone_to_state_scale=hist_state_scale,
            promoter_incoming_edge_scale=incoming_prom_scale,
            sign_repressive_to_state=sign_repr,
            repressive_genes_tsv=repr_tsv,
        )
        print(
            f"[hme-promoter-state] edges {promoter_edge_stats.n_in}->{promoter_edge_stats.n_out} "
            f"drop_HME_to_state={promoter_edge_stats.n_dropped_enzyme_to_state} "
            f"genes={len(promoter_edge_stats.dropped_enzyme_genes or [])} "
            f"list={promoter_edge_stats.enzyme_list_source} "
            f"hist_to_state_scaled={promoter_edge_stats.n_scaled_hist_to_state} "
            f"hist_scale={promoter_edge_stats.hist_scale} "
            f"incoming_promoter_scaled={promoter_edge_stats.n_scaled_incoming_promoter} "
            f"incoming_scale={promoter_edge_stats.incoming_promoter_scale} "
            f"repressive_signed={promoter_edge_stats.n_repressive_signed} "
            f"repr_genes={promoter_edge_stats.repressive_signed_genes}",
            flush=True,
        )
    ups_policy_stats = None
    skip_protein_degradation: set[str] = set()
    use_ups = (
        bool(getattr(pol, "ups_degradation_policy", True))
        if ups_degradation_policy is None
        else bool(ups_degradation_policy)
    )
    if use_ups:
        from gin_tool.subgraph_pipeline.ups_degradation_policy import (
            DEFAULT_E3_MEMBERS,
            DEFAULT_ESI_EVIDENCE,
            apply_ups_edge_scale,
        )

        weak_scale = weak_e3_edge_scale_val
        strong_inhib_w = strong_e3_inhibition_weight_val
        members_tsv = ups_e3_members_tsv or getattr(pol, "ups_e3_members_tsv", None)
        evidence_tsv = ups_esi_evidence_tsv or getattr(pol, "ups_esi_evidence_tsv", None)
        override_tsv = ups_e3_pair_override_tsv or getattr(
            pol, "ups_e3_pair_override_tsv", None
        )
        # Soft-skip if seed tables missing (non-p53 projects).
        mem_path = Path(members_tsv) if members_tsv else DEFAULT_E3_MEMBERS
        ev_path = Path(evidence_tsv) if evidence_tsv else DEFAULT_ESI_EVIDENCE
        if not mem_path.is_file() or not ev_path.is_file():
            print(
                f"[ups-degradation-policy] SKIP missing tables "
                f"members={mem_path.is_file()} evidence={ev_path.is_file()}",
                flush=True,
            )
        else:
            gin_edges, ups_policy_stats, skip_protein_degradation = apply_ups_edge_scale(
                gin_edges,
                weak_scale=weak_scale,
                strong_e3_inhibition_weight=strong_inhib_w,
                e3_members_tsv=members_tsv,
                esi_evidence_tsv=evidence_tsv,
                override_tsv=override_tsv,
            )
            print(
                f"[ups-degradation-policy] pairs={ups_policy_stats.n_ups_pairs} "
                f"strong={ups_policy_stats.n_strong_pairs} weak={ups_policy_stats.n_weak_pairs} "
                f"arcs_scaled={ups_policy_stats.n_arcs_scaled} weak_scale={weak_scale} "
                f"strong_inhib_w={strong_inhib_w} "
                f"strong_arcs_scaled={ups_policy_stats.n_strong_arcs_scaled} "
                f"skip_global_protein={len(skip_protein_degradation)}",
                flush=True,
            )
    non_e3_stats = None
    use_non_e3 = (
        bool(getattr(pol, "scale_non_e3_activity_inhibitors", True))
        if scale_non_e3_activity_inhibitors is None
        else bool(scale_non_e3_activity_inhibitors)
    )
    if use_non_e3:
        from gin_tool.subgraph_pipeline.ups_degradation_policy import (
            apply_non_e3_activity_inhibitor_scale,
            load_esi_evidence_tables,
        )

        ne_path = non_e3_activity_inhibitors_tsv or getattr(
            pol, "non_e3_activity_inhibitors_tsv", None
        )
        ne_scale = (
            float(getattr(pol, "non_e3_activity_inhibitor_scale", 0.1))
            if non_e3_activity_inhibitor_scale is None
            else float(non_e3_activity_inhibitor_scale)
        )
        e3_mem = None
        try:
            e3_mem = load_esi_evidence_tables().e3_members
        except Exception:  # noqa: BLE001
            e3_mem = frozenset()
        gin_edges, non_e3_stats = apply_non_e3_activity_inhibitor_scale(
            gin_edges,
            gene_list_tsv=ne_path,
            default_scale=ne_scale,
            e3_members=e3_mem,
        )
        print(
            f"[non-e3-activity-inhibitor] arcs_scaled={non_e3_stats.n_arcs_scaled} "
            f"scale={ne_scale} genes_hit={non_e3_stats.genes_hit}",
            flush=True,
        )
    drop_inhib_stats = None
    if drop_protein_inhibition_pairs:
        from gin_tool.subgraph_pipeline.ups_degradation_policy import (
            DROP_TP53_INHIBITOR_PAIRS,
            drop_protein_inhibition_pairs as _drop_inh,
        )

        pairs = (
            DROP_TP53_INHIBITOR_PAIRS
            if drop_protein_inhibition_pairs is True
            else list(drop_protein_inhibition_pairs)
        )
        gin_edges, drop_inhib_stats = _drop_inh(gin_edges, pairs=pairs)
        print(
            f"[drop-inhibition] pairs={drop_inhib_stats.pairs} "
            f"dropped={drop_inhib_stats.n_dropped} "
            f"edges {drop_inhib_stats.n_in}->{drop_inhib_stats.n_out}",
            flush=True,
        )
    if force_skip_protein_degradation:
        skip_protein_degradation |= {str(p) for p in force_skip_protein_degradation}
    self_loop_stats = None
    if drop_protein_self_loops:
        from gin_tool.subgraph_pipeline.protein_self_loop_prune import (
            drop_protein_self_loops as _drop_protein_self_loops,
        )

        gin_edges, self_loop_stats = _drop_protein_self_loops(gin_edges)
        print(
            f"[protein-self-loop] removed={self_loop_stats.n_removed} "
            f"edges {self_loop_stats.n_in}->{self_loop_stats.n_out}",
            flush=True,
        )
    bal_layers: list[str]
    if isinstance(edge_balance, (list, tuple, set)):
        bal_layers = parse_edge_balance_layers(",".join(str(x) for x in edge_balance))
    else:
        bal_layers = parse_edge_balance_layers(edge_balance)
    bal_stats = None
    if bal_layers:
        gin_edges, bal_stats = balance_gin_edge_weights(gin_edges, layers=bal_layers)
        print(
            f"[edge-balance] layers={bal_layers} places={bal_stats.n_places_considered} "
            f"balanced={bal_stats.n_balanced} missing_side={bal_stats.n_skipped_missing_side} "
            f"edges_scaled={bal_stats.n_edges_scaled}",
            flush=True,
        )
    seeds = seed_gin_names_from_manifest(manifest_path) if manifest_path.is_file() else set()
    try:
        from gin_tool.histone_promoter_sign_qc import (
            audit_histone_promoter_signs,
            write_histone_promoter_sign_qc,
        )

        hist_qc = audit_histone_promoter_signs(gin_edges)
        write_histone_promoter_sign_qc(subgraph_dir / "histone_promoter_sign_qc.json", hist_qc)
        print(
            f"[histone-qc] hist→promoter;state inhibitory={hist_qc.n_inhibitory} "
            f"activating={hist_qc.n_activating} paths_to_active={hist_qc.n_paths_to_promoter_active} "
            f"(inhibitory_paths={hist_qc.n_inhibitory_paths_to_active})",
            flush=True,
        )
    except Exception as exc:  # noqa: BLE001 — best-effort QC
        print(f"[histone-qc] skip: {exc}", flush=True)
        hist_qc = None
    print(
        f"[spn] build petri net edges={len(gin_edges)} seeds={len(seeds)} "
        f"mark_mode={mark_mode} degradation={add_degradation} (serial; "
        f"ProcessPool n_jobs starts after export)",
        flush=True,
    )
    net = build_petri_net_from_bipartite_edges(
        gin_edges, seed_gin_names=seeds, mark_mode=mark_mode
    )
    levels_in_net: dict[str, float] | None = None
    if place_levels:
        levels_in_net = {p: float(v) for p, v in place_levels.items() if p in net.places and float(v) > 0}
        marked = apply_place_token_levels(net, levels_in_net)
        mark_mode = f"place_levels({mark_mode})"
        print(
            f"[spn] place_levels override: requested={len(place_levels)} "
            f"in_net={len(levels_in_net)} marked={len(marked)}",
            flush=True,
        )
    print(
        f"[spn] petri net places={len(net.places)} transitions={len(net.transitions)} "
        f"arcs={len(net.arcs)} marked={sum(1 for v in net.place_tokens.values() if v > 0)}",
        flush=True,
    )

    maintain_met_places: frozenset[str] = frozenset()
    if maintain_metabolite_sources and not maintain_seeds:
        from gin_tool.subgraph_pipeline.met_source_seed import (
            metabolite_source_places_in_net,
        )

        maintain_met_places = frozenset(
            metabolite_source_places_in_net(gin_edges, set(net.places))
        )
        print(
            f"[met-maintain] in-degree-0 met sources in net: {len(maintain_met_places)}",
            flush=True,
        )

    deg_stats = None
    if add_degradation:
        v2 = layer_rules == LAYER_RULES_V2
        deg_stats = add_molecule_degradation_layer(
            net,
            rna_half_life_blocks=rna_half_life_blocks,
            mirna_half_life_blocks=mirna_half_life_blocks,
            protein_half_life_blocks=protein_half_life_blocks,
            histone_half_life_blocks=(
                getattr(pol, "histone_half_life_blocks", 30.0) if v2 else None
            ),
            promoter_active_half_life_blocks=(
                getattr(pol, "promoter_active_half_life_blocks", 20.0) if v2 else None
            ),
            # Never blanket-degrade all metabolites unless policy explicitly sets a half-life.
            metabolite_half_life_blocks=getattr(pol, "metabolite_half_life_blocks", None),
            skip_protein_places=skip_protein_degradation or None,
        )

    graphml_path = subgraph_dir / "spn.graphml"
    write_mepn_graphml(graphml_path, net)
    qc = {
        "spn_pipeline_version": pol.version,
        "layer_rules": layer_rules,
        "stoich": stoich,
        "petri_net": net.to_dict(),
        "graphml": str(graphml_path),
        "maintain_seeds": maintain_seeds,
        "maintain_metabolite_sources": maintain_metabolite_sources,
        "n_maintain_met_sources": len(maintain_met_places),
        "maintain_met_sources": sorted(maintain_met_places),
        "mark_mode": mark_mode,
        "add_degradation": add_degradation,
        "keep_trajectories": keep_trajectories,
        "n_inhibitory_arc_weights": sum(
            1 for w in (getattr(net, "arc_weights", {}) or {}).values() if w < 0
        ),
    }
    if deg_stats is not None:
        qc["degradation"] = deg_stats.to_dict()
    if bal_stats is not None:
        qc["edge_balance"] = bal_stats.to_dict()
    if layer_rules == LAYER_RULES_V2 or stoich in ("one_to_one", "kegg"):
        try:
            from gin_tool.spn_export import SpnPlaceKind as _SPK
            from gin_tool.subgraph_pipeline.stoichiometry_resolve import report_stoich_coverage

            place_is_met = {
                p
                for p, k in net.places.items()
                if k == _SPK.METABOLITE
                or str(p).startswith(("CHEBI", "SID"))
                or str(p).endswith(".met")
            }
            met_tids: list[str] = []
            for src, tgt in net.arcs:
                if src in place_is_met and tgt in net.transitions:
                    met_tids.append(tgt)
                elif tgt in place_is_met and src in net.transitions:
                    met_tids.append(src)
            met_tids = sorted(set(met_tids))
            cov = report_stoich_coverage(met_tids, met_transition_ids=met_tids)
            qc["stoich_coverage"] = {**cov.to_dict(), "mode": stoich}
            print(
                f"[stoich] mode={stoich} met_trans={cov.n_met_edges} "
                f"keggr={cov.n_keggr_tagged} resolved={cov.n_resolved} "
                f"fallback_1to1={cov.n_fallback_one_to_one}",
                flush=True,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"[stoich] skip coverage: {exc}", flush=True)
    if hist_qc is not None:
        qc["histone_promoter_signs"] = hist_qc.to_dict()
    (subgraph_dir / "spn_export_qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")

    model = petri_net_to_spn_model(net, add_self_inhibition=add_self_inhibition)
    try:
        from gin_tool.signed_edge_spn_qc import audit_signed_edges_in_spn, write_signed_edge_spn_qc

        signed_qc = audit_signed_edges_in_spn(
            gin_edges,
            arc_weights=getattr(net, "arc_weights", {}) or {},
            arcs=net.arcs,
            interactions=model.interactions,
        )
        write_signed_edge_spn_qc(subgraph_dir / "signed_edge_spn_qc.json", signed_qc)
        print(
            f"[signed-qc] neg_gin={signed_qc.n_neg_gin_edges} neg_arcs={signed_qc.n_neg_arcs_in_petri} "
            f"dropped={signed_qc.n_neg_dropped} rvr_inh={signed_qc.n_got_rvr_inhibition}/"
            f"{signed_qc.n_expected_rvr_inhibition} both_neg={signed_qc.n_both_negative_rvr} "
            f"ok={signed_qc.to_dict()['ok']}",
            flush=True,
        )
        qc["signed_edge_spn"] = signed_qc.to_dict()
        (subgraph_dir / "spn_export_qc.json").write_text(
            json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8"
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[signed-qc] skip: {exc}", flush=True)
    print(
        f"[spn] model places={len(model.places)} interactions={len(model.interactions)} "
        f"fixed_high={len(model.fixed_high)}; starting MC runs={runs} n_jobs={n_jobs}",
        flush=True,
    )
    sim_dir = subgraph_dir / sim_subdir
    sim_dir.mkdir(parents=True, exist_ok=True)
    chromhmm_weights: dict[str, float] | None = None
    chromhmm_meta: dict | None = None
    if use_v8:
        from gin_tool.subgraph_pipeline.chromhmm_emission_weights import (
            resolve_chromhmm_mark_weights,
            write_mark_weights_json,
            write_mark_weights_tsv,
        )

        emw = resolve_chromhmm_mark_weights(
            emission_tsv=chromhmm_emission_tsv,
            weights_tsv=chromhmm_mark_weights_tsv,
        )
        if emw is None:
            raise ValueError("spn_v8 requested but no chromHMM emission/weights resolved")
        chromhmm_weights = dict(emw.weights)
        chromhmm_meta = emw.to_dict()
        write_mark_weights_json(emw, sim_dir / f"{stem}_chromhmm_mark_weights.json")
        write_mark_weights_tsv(emw.weights, sim_dir / f"{stem}_chromhmm_mark_weights.tsv")
        print(
            f"[chromhmm-emw] marks={len(emw.weights)} "
            f"promoter_states={emw.n_promoter_states} "
            f"A={emw.n_active_states} I={emw.n_inactive_states} "
            f"source={emw.source_path}",
            flush=True,
        )
        qc["chromhmm_emw"] = chromhmm_meta
        (subgraph_dir / "spn_export_qc.json").write_text(
            json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8"
        )
    summary = run_spn_simulation(
        model,
        out_dir=sim_dir,
        stem=stem,
        policy=pol,
        blocks=blocks,
        runs=runs,
        seed=seed,
        maintain_seeds=maintain_seeds,
        maintain_every_n_blocks=maintain_every_n_blocks,
        maintain_places=maintain_met_places or None,
        maintain_levels=levels_in_net,
        keep_trajectories=keep_trajectories,
        n_jobs=n_jobs,
        chromhmm_mark_weights=chromhmm_weights,
        chromhmm_weight_meta=chromhmm_meta,
        interaction_weight_pulse=getattr(pol, "interaction_weight_pulse", None),
        nf_qc_pair=nf_qc_pair,
        nf_warm=nf_warm,
        nf_max_lag=nf_max_lag,
        regulatory_norm=getattr(pol, "regulatory_norm", None),
        regulatory_hill_k=getattr(pol, "regulatory_hill_k", None),
        regulatory_target_budget=getattr(pol, "regulatory_target_budget", None),
    )
    summary["add_degradation"] = add_degradation
    summary["mark_mode"] = mark_mode
    summary["graphml"] = str(graphml_path)
    if bal_stats is not None:
        summary["edge_balance"] = bal_stats.to_dict()
    if self_loop_stats is not None:
        summary["protein_self_loop_prune"] = self_loop_stats.to_dict()
    if promoter_edge_stats is not None:
        summary["promoter_edge_policy"] = promoter_edge_stats.to_dict()
    if ups_policy_stats is not None:
        summary["ups_degradation_policy"] = ups_policy_stats.to_dict()
    if non_e3_stats is not None:
        summary["non_e3_activity_inhibitor_scale"] = non_e3_stats.to_dict()
    if drop_inhib_stats is not None:
        summary["drop_protein_inhibition"] = drop_inhib_stats.to_dict()
    if deg_stats is not None:
        summary["degradation"] = deg_stats.to_dict()
    summary["fire_propensity_mode"] = getattr(pol, "fire_propensity_mode", "off")
    summary["fire_propensity_scale"] = getattr(pol, "fire_propensity_scale", 1.0)
    summary["fire_propensity_ref"] = (
        float(getattr(pol, "fire_propensity_ref", 0.0) or pol.high_tokens)
    )
    if interaction_weight_pulse:
        summary["interaction_weight_pulse"] = dict(interaction_weight_pulse)
    summary["tf_dna_route_independent"] = bool(getattr(pol, "tf_dna_route_independent", True))
    summary["tf_dna_edge_weight"] = float(
        getattr(pol, "tf_dna_edge_weight", getattr(pol, "tf_dna_fixed_drive", 0.5))
    )
    summary["promoter_incoming_edge_scale"] = float(
        getattr(pol, "promoter_incoming_edge_scale", 1.0)
    )
    summary["prod_m_max"] = float(getattr(pol, "prod_m_max", 50.0))
    summary["prod_alpha_transcription"] = float(getattr(pol, "prod_alpha_transcription", 1.0))
    summary["translation_output_scale"] = float(getattr(pol, "translation_output_scale", 1.0))
    summary["protein_production_output_scale"] = float(
        getattr(pol, "protein_production_output_scale", 1.0)
    )
    if levels_in_net is not None:
        summary["place_levels_n"] = len(levels_in_net)
    summary["maintain_metabolite_sources"] = maintain_metabolite_sources
    if maintain_met_places:
        summary["maintain_met_sources_n"] = len(maintain_met_places)
    if keep_trajectories:
        summary = _run_reasonableness_qc_after_sim(
            sim_dir=sim_dir,
            stem=stem,
            subgraph_dir=subgraph_dir,
            model=model,
            summary=summary,
        )
    return summary

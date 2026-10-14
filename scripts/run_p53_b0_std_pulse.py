#!/usr/bin/env python3
"""H1 B0 standard Pulse — publication SPECIAL CASE (MDM2⊣TP53 on H1 SigPN).

SigPN = signaling Petri net (not stochastic PN). This script hard-codes the
manuscript B0 protocol; other interventions need a custom driver that reuses
``export_and_simulate_subgraph_dir`` / ``spn_standard_policy``.
See docs/SIGPN_CUSTOM_SIMULATIONS.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from gin_tool.subgraph_pipeline.spn_policy import SPN_PIPELINE_VERSION_STD
from gin_tool.subgraph_pipeline.spn_standard_policy import (
    STD_B0_BLOCKS,
    STD_B0_HOLD,
    STD_B0_ONSET,
    STD_B0_RUNS,
    STD_FAST_E_TIMESCALE,
    STD_PROMOTER_INCOMING_SCALE,
    STD_PROD_M_MAX,
    STD_ALPHA_TX,
    STD_STRONG_E3,
    STD_TF_DNA_EDGE_W_BASE,
    STD_TRANSLATION_SCALE,
    STD_V15C_OUT_ROUTE_HL,
    STD_V15C_PROMOTER_ACTIVE_HL,
    STD_V15C_PROTEIN_HL,
    STD_V15C_RNA_HL,
    STD_WEAK_E3,
    fast_e_half_life,
    standard_b0_pulse_step,
    standard_h1_b0_simulate_kwargs,
)

ROOT = Path(__file__).resolve().parents[1]
PRUNE = ROOT / "gin_out/p53_cell_grn_prune_v4_seedhist/h1_wa01"
LEV = ROOT / "runtime/omics_seed/cellline_levels/h1_wa01/h1_wa01_rna_plus_met_sources_place_levels.tsv"
EMW = ROOT / "runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv"
OUT = ROOT / "runtime/spn_phys5"

TAG = "std"
PROTEIN_HL = fast_e_half_life(STD_V15C_PROTEIN_HL)
RNA_HL = fast_e_half_life(STD_V15C_RNA_HL)
PROMOTER_ACTIVE_HL = fast_e_half_life(STD_V15C_PROMOTER_ACTIVE_HL)
OUT_ROUTE_HL = fast_e_half_life(STD_V15C_OUT_ROUTE_HL)
TRACK = (
    "MDM2.promoter.active",
    "MDM2.RNA",
    "MDM2.protein",
    "TP53.protein",
)
NF_PAIR = ("TP53.protein", "MDM2.protein")


def run_arm(
    *,
    arm: str,
    runs: int,
    jobs: int,
    force: bool,
    hold: int,
    onset: int,
    blocks: int,
    seed: int = 42,
) -> dict:
    from gin_tool.subgraph_pipeline.spn_simulate import (
        export_and_simulate_subgraph_dir,
        load_place_levels_tsv,
    )

    sub = f"spn_sim_std_B0ind_{TAG}_h{hold}_{arm}_oneshot"
    stem = f"p53_B0ind_{TAG}_h{hold}_{arm}_once_b{blocks}_n{runs}"
    sim_dir = PRUNE / sub
    traj = sim_dir / f"{stem}_trajectory.json"
    if not force and traj.is_file():
        return {
            "arm": arm,
            "skipped": True,
            "sim_subdir": sub,
            "stem": stem,
            "trajectory_json": str(traj),
            "pair_runs_npz": str(sim_dir / f"{stem}_pair_runs.npz"),
        }

    print(
        f"[std {arm}] onset={onset} hold={hold} blocks={blocks} n={runs} "
        f"pipeline={SPN_PIPELINE_VERSION_STD} fast×e α_tx={STD_ALPHA_TX} "
        f"strong_e3={STD_STRONG_E3} weak_e3={STD_WEAK_E3:.4f} …",
        flush=True,
    )
    kw = standard_h1_b0_simulate_kwargs(
        chromhmm_emission_tsv=str(EMW),
        blocks=blocks,
        runs=runs,
        n_jobs=jobs,
        seed=int(seed),
    )
    kw.update(
        place_levels=load_place_levels_tsv(LEV),
        sim_subdir=sub,
        stem=stem,
        keep_trajectories=True,
        nf_qc_pair=NF_PAIR,
    )
    if arm == "step":
        kw["interaction_weight_pulse"] = standard_b0_pulse_step(onset=onset, hold=hold)
    summary = export_and_simulate_subgraph_dir(PRUNE, **kw)
    return {
        "arm": arm,
        "skipped": False,
        "sim_subdir": sub,
        "stem": stem,
        "trajectory_json": str(traj),
        "pair_runs_npz": summary.get(
            "pair_runs_npz", str(sim_dir / f"{stem}_pair_runs.npz")
        ),
        "spn_pipeline_version": summary.get("spn_pipeline_version"),
        "promoter_incoming_edge_scale": summary.get("promoter_incoming_edge_scale"),
        "prod_m_max": summary.get("prod_m_max"),
        "prod_alpha_transcription": summary.get("prod_alpha_transcription"),
        "translation_output_scale": summary.get("translation_output_scale"),
        "tf_dna_edge_weight": summary.get("tf_dna_edge_weight"),
        "strong_e3_inhibition_weight": STD_STRONG_E3,
        "weak_e3_edge_scale": STD_WEAK_E3,
        "tf_dna_edge_weight_base": STD_TF_DNA_EDGE_W_BASE,
    }


def follow_vs_control(ctrl_npz: Path, step_npz: Path, *, onset: int, hold: int) -> dict:
    zc = np.load(ctrl_npz)
    zs = np.load(step_npz)
    places = [str(x) for x in zc["places"]]
    ia, ib = places.index("TP53.protein"), places.index("MDM2.protein")
    C, S = zc["block_samples"], zs["block_samples"]
    hold_sl = slice(onset, onset + hold)
    post_sl = slice(onset + hold, min(onset + hold + 40, S.shape[1]))
    d_tp = S[:, hold_sl, ia].mean(1) - C[:, hold_sl, ia].mean(1)
    d_md = S[:, post_sl, ib].mean(1) - C[:, post_sl, ib].mean(1)
    follow = (d_tp > 0) & (d_md > 0)
    return {
        "mean_dTP53_hold_vs_ctrl": float(d_tp.mean()),
        "mean_dMDM2_post_vs_ctrl": float(d_md.mean()),
        "median_dTP53_hold": float(np.median(d_tp)),
        "median_dMDM2_post": float(np.median(d_md)),
        "rate_follow_vs_ctrl": float(follow.mean()),
        "rate_MDM2_post_up": float((d_md > 0).mean()),
        "rate_TP53_hold_up": float((d_tp > 0).mean()),
    }


def extract_track_means(traj_path: Path) -> dict[str, np.ndarray]:
    d = json.loads(traj_path.read_text(encoding="utf-8"))
    places = list(d["places"])
    mean = np.asarray(d["mean_markings"], float)
    return {name: mean[:, places.index(name)].copy() for name in TRACK}


def control_dynamics_summary(traj_path: Path, *, pre_end: int) -> dict:
    means = extract_track_means(traj_path)
    out: dict[str, dict] = {}
    pre = slice(0, pre_end)
    for name in ("TP53.protein", "MDM2.protein"):
        y = means[name]
        y_pre = y[pre]
        peak_i = int(np.argmax(y_pre))
        out[name] = {
            "pre_mean": float(y_pre.mean()),
            "pre_end_vs_start": float(y[pre_end - 1] / max(y[0], 1e-9)),
            "pre_peak_block": peak_i,
            "pre_peak_value": float(y[peak_i]),
            "late_end_vs_pre_mean": float(y[-1] / max(y_pre.mean(), 1e-9)),
            "monotone_pre": bool(np.all(np.diff(y_pre) >= -1e-9)),
        }
    return out


def layer_coupling_step_minus_ctrl(ctrl_traj: Path, step_traj: Path, *, onset: int, hold: int) -> dict:
    c = extract_track_means(ctrl_traj)
    s = extract_track_means(step_traj)
    hold_sl = slice(onset, onset + hold)
    post_sl = slice(onset + hold, onset + hold + 40)
    return {
        name: {
            "hold_delta": float(s[name][hold_sl].mean() - c[name][hold_sl].mean()),
            "post_delta": float(s[name][post_sl].mean() - c[name][post_sl].mean()),
        }
        for name in TRACK
    }


def save_tracks(by_arm: dict, hold: int, onset: int, blocks: int, follow: dict) -> dict[str, str]:
    OUT.mkdir(parents=True, exist_ok=True)
    prefix = OUT / f"p53_B0_{TAG}_h{hold}"
    arrays: dict[str, np.ndarray] = {}
    series: dict = {
        "tag": TAG,
        "spn_pipeline_version": SPN_PIPELINE_VERSION_STD,
        "timescale_factor": STD_FAST_E_TIMESCALE,
        "hold": hold,
        "onset": onset,
        "half_lives": {
            "protein": PROTEIN_HL,
            "rna": RNA_HL,
            "promoter_active": PROMOTER_ACTIVE_HL,
            "out_route": OUT_ROUTE_HL,
        },
        "prod_alpha_transcription": STD_ALPHA_TX,
        "promoter_incoming_edge_scale": STD_PROMOTER_INCOMING_SCALE,
        "prod_m_max": STD_PROD_M_MAX,
        "places": list(TRACK),
        "arms": {},
        "follow_vs_control": follow,
    }
    n_t = None
    for arm, meta in by_arm.items():
        means = extract_track_means(Path(meta["trajectory_json"]))
        series["arms"][arm] = {
            "summary": {
                k: {
                    "hold_mean": float(v[onset : onset + hold].mean()),
                    "post_40": float(v[onset + hold : onset + hold + 40].mean()),
                    "late_post_to_end": float(
                        v[min(onset + hold + 40, len(v) - 1) :].mean()
                    ),
                }
                for k, v in means.items()
            }
        }
        for k, v in means.items():
            arrays[f"{arm}__{k}"] = v
            n_t = len(v) if n_t is None else n_t
    block = np.arange(n_t or (blocks + 1), dtype=np.int32)
    arrays["block"] = block
    npz = Path(str(prefix) + "_track_means.npz")
    np.savez_compressed(npz, **arrays)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), dpi=140, sharex=True)
    order = ("TP53.protein", "MDM2.protein", "MDM2.RNA", "MDM2.promoter.active")
    colors = {"control": "#4c78a8", "step": "#e45756"}
    for ax, name in zip(axes.ravel(), order):
        for arm in by_arm:
            ax.plot(block, arrays[f"{arm}__{name}"], color=colors[arm], lw=1.8, label=arm)
        ax.axvspan(onset, onset + hold, color="0.85", alpha=0.55, zorder=0)
        ax.set_title(name)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle(
        f"std onset={onset} blocks={blocks} ΔTP53={follow['mean_dTP53_hold_vs_ctrl']:+.2f} "
        f"ΔMDM2_post={follow['mean_dMDM2_post_vs_ctrl']:+.3f} "
        f"follow={follow['rate_follow_vs_ctrl']:.0%}",
        fontsize=10,
    )
    fig.tight_layout()
    png = Path(str(prefix) + "_track_means.png")
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    jpath = Path(str(prefix) + "_track_means.json")
    jpath.write_text(json.dumps(series, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"npz": str(npz), "png": str(png), "json": str(jpath)}


def main() -> int:
    global PRUNE, LEV, EMW, OUT, TAG

    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, default=STD_B0_RUNS)
    ap.add_argument("--n-jobs", type=int, default=20)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--onset", type=int, default=STD_B0_ONSET)
    ap.add_argument("--blocks", type=int, default=STD_B0_BLOCKS)
    ap.add_argument("--hold", type=int, default=STD_B0_HOLD)
    ap.add_argument(
        "--subgraph-dir",
        type=Path,
        default=None,
        help="H1 pruned cell dir with subgraph_edges.tsv (default: historical seedhist h1)",
    )
    ap.add_argument(
        "--tag",
        default=None,
        help="Report/sim tag suffix (default: std). Use a distinct tag for Stage A repro.",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Report output dir (default: runtime/spn_phys5)",
    )
    ap.add_argument("--levels", type=Path, default=None, help="Place levels TSV")
    ap.add_argument("--emw", type=Path, default=None, help="chromHMM emission TSV")
    ap.add_argument(
        "--seed",
        type=int,
        default=42,
        help="RNG seed for SPN stochastic runs (default 42; prior gin11a pulse used 42)",
    )
    args = ap.parse_args()
    if args.subgraph_dir is not None:
        PRUNE = args.subgraph_dir.resolve()
    if args.levels is not None:
        LEV = args.levels.resolve()
    if args.emw is not None:
        EMW = args.emw.resolve()
    if args.out_dir is not None:
        OUT = args.out_dir.resolve()
    if args.tag is not None:
        TAG = str(args.tag).strip() or TAG
    hold = int(args.hold)
    onset = int(args.onset)
    blocks = int(args.blocks)
    OUT.mkdir(parents=True, exist_ok=True)
    if not (PRUNE / "subgraph_edges.tsv").is_file():
        raise SystemExit(f"missing subgraph_edges.tsv under {PRUNE}")

    by = {
        arm: run_arm(
            arm=arm,
            runs=args.runs,
            jobs=args.n_jobs,
            force=args.force,
            hold=hold,
            onset=onset,
            blocks=blocks,
            seed=int(args.seed),
        )
        for arm in ("control", "step")
    }

    follow = follow_vs_control(
        Path(by["control"]["pair_runs_npz"]),
        Path(by["step"]["pair_runs_npz"]),
        onset=onset,
        hold=hold,
    )
    coupling = layer_coupling_step_minus_ctrl(
        Path(by["control"]["trajectory_json"]),
        Path(by["step"]["trajectory_json"]),
        onset=onset,
        hold=hold,
    )
    ctrl_dyn = control_dynamics_summary(
        Path(by["control"]["trajectory_json"]), pre_end=onset
    )
    paths = save_tracks(by, hold, onset, blocks, follow)
    dmd = follow["mean_dMDM2_post_vs_ctrl"]
    dtp = follow["mean_dTP53_hold_vs_ctrl"]
    rate = follow["rate_follow_vs_ctrl"]

    report = {
        "tag": TAG,
        "pipeline": SPN_PIPELINE_VERSION_STD,
        "legacy_alias": "v15cfea10e3ewe_o150b350_metsrc",
        "subgraph_dir": str(PRUNE),
        "place_levels": str(LEV),
        "chromhmm_emission_tsv": str(EMW),
        "blocks": blocks,
        "strong_e3_inhibition_weight": STD_STRONG_E3,
        "weak_e3_edge_scale": STD_WEAK_E3,
        "timescale_factor": STD_FAST_E_TIMESCALE,
        "half_lives": {
            "protein_half_life_blocks": PROTEIN_HL,
            "rna_half_life_blocks": RNA_HL,
            "promoter_active_half_life_blocks": PROMOTER_ACTIVE_HL,
            "out_route_half_life_blocks": OUT_ROUTE_HL,
        },
        "promoter_incoming_edge_scale": STD_PROMOTER_INCOMING_SCALE,
        "prod_m_max": STD_PROD_M_MAX,
        "prod_alpha_transcription": STD_ALPHA_TX,
        "translation_output_scale": STD_TRANSLATION_SCALE,
        "hold": hold,
        "onset": onset,
        "runs": args.runs,
        "seed": int(args.seed),
        "arms": by,
        "follow_vs_control": follow,
        "control_dynamics_pre_pulse": ctrl_dyn,
        "layer_coupling_step_minus_ctrl": coupling,
        "track_outputs": paths,
        "verdict_plain": (
            f"std onset={onset} blocks={blocks} hold={hold}: "
            f"ΔTP53={dtp:+.3f} follow={rate:.1%} ΔMDM2_post={dmd:+.3f}"
        ),
    }
    rp = OUT / f"p53_B0_{TAG}_h{hold}_report.json"
    rp.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(report["verdict_plain"], flush=True)
    print(f"[report] {rp}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

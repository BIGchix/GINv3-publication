#!/usr/bin/env python3
"""Fig4 C/D only: old fig4-2 noprune pipeline × frozen Stage A gin11a context.

- Shared Stage A edges (NO cell GRN prune)
- spn_v5_degrade_anc_source / DEFAULT_SPN_POLICY
- 200 runs × 200 blocks, maintain_seeds=True
- Unified RNA(+met-source) place-levels per cell (same file family for all 6 lines)
- Heatmaps: C=histone_mark (wine), D=promoter.active (red)
- Does not touch Fig4B or pulse artifacts.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from gin_tool.subgraph_pipeline.spn_policy import DEFAULT_SPN_POLICY, SPN_PIPELINE_VERSION
from gin_tool.subgraph_pipeline.spn_simulate import (
    export_and_simulate_subgraph_dir,
    load_place_levels_tsv,
)
from run_fig4bcd_gin11a_from_stageA import (  # noqa: E402
    CELLS,
    _ensure_colormaps,
    _load_traj_plateau,
    _pearson_matrix,
    _place_layer,
    plot_cell_corr_heatmap,
)

# Same file family for every cell (avoid H1 seed-only rna_place_levels asymmetry).
LEVELS_SUFFIX = "rna_plus_met_sources_place_levels.tsv"
STEM_TAG = "metsrc"


def _prepare_work_dir(*, stage_a: Path, work: Path) -> Path:
    work.mkdir(parents=True, exist_ok=True)
    src_e = stage_a / "subgraph_edges.tsv"
    src_m = stage_a / "seed_manifest.json"
    if not src_e.is_file():
        raise FileNotFoundError(src_e)
    dst_e = work / "subgraph_edges.tsv"
    if not dst_e.is_file() or dst_e.stat().st_size != src_e.stat().st_size:
        shutil.copy2(src_e, dst_e)
    if src_m.is_file() and not (work / "seed_manifest.json").is_file():
        shutil.copy2(src_m, work / "seed_manifest.json")
    meta = {
        "pipeline": "fig4cd_noprune_gin11a_stageA",
        "spn_pipeline_version": SPN_PIPELINE_VERSION,
        "stage_a": str(stage_a.relative_to(ROOT)),
        "note": "Shared context; no GRN prune; matches fig4-2 C/D method with new Stage A input",
    }
    (work / "RUN_META.json").write_text(json.dumps(meta, indent=2) + "\n")
    return work


def _levels_path(cell_id: str) -> Path:
    p = ROOT / "runtime/omics_seed/cellline_levels" / cell_id / f"{cell_id}_{LEVELS_SUFFIX}"
    if not p.is_file():
        raise FileNotFoundError(p)
    return p


def run_cell(
    *,
    work: Path,
    cell_id: str,
    runs: int,
    blocks: int,
    n_jobs: int,
    seed: int,
    force: bool,
) -> Path:
    levels = _levels_path(cell_id)
    stem = f"p53_noprune_gin11a_{cell_id}_{STEM_TAG}_b{blocks}_n{runs}"
    traj = work / "spn_sim" / f"{stem}_trajectory.json"
    if traj.is_file() and not force:
        print(f"[sim] skip existing {traj}", flush=True)
        return traj
    print(
        f"[sim] {cell_id} spn={SPN_PIPELINE_VERSION} blocks={blocks} runs={runs} "
        f"maintain_seeds={DEFAULT_SPN_POLICY.maintain_seeds} levels={levels.name}",
        flush=True,
    )
    export_and_simulate_subgraph_dir(
        work,
        policy=DEFAULT_SPN_POLICY,
        blocks=blocks,
        runs=runs,
        n_jobs=n_jobs,
        seed=seed,
        maintain_seeds=True,
        mark_mode="ancestor_all",
        place_levels=load_place_levels_tsv(levels),
        keep_trajectories=True,
        sim_subdir="spn_sim",
        stem=stem,
        # Explicitly off: do not inherit std B0 / chromHMM knobs
        layer_rules="off",
        chromhmm_emission_tsv=None,
    )
    if not traj.is_file():
        raise FileNotFoundError(traj)
    return traj


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stage-a",
        type=Path,
        default=ROOT / "gin_out/p53_shared_context_v4_rh1_gin11a",
    )
    ap.add_argument(
        "--work-dir",
        type=Path,
        default=ROOT / "gin_out/p53_fig4_gin11a_noprune_cd",
    )
    ap.add_argument("--fig-dir", type=Path, default=ROOT / "manuscript/figures/fig4")
    ap.add_argument("--blocks", type=int, default=200)
    ap.add_argument("--runs", type=int, default=200)
    ap.add_argument("--n-jobs", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--plateau-blocks", type=int, default=100)
    ap.add_argument("--skip-sim", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    stage_a = args.stage_a if args.stage_a.is_absolute() else ROOT / args.stage_a
    work = args.work_dir if args.work_dir.is_absolute() else ROOT / args.work_dir
    fig_dir = args.fig_dir if args.fig_dir.is_absolute() else ROOT / args.fig_dir
    work = _prepare_work_dir(stage_a=stage_a, work=work)

    labels = [lab for _, lab, _ in CELLS]
    trajs: dict[str, Path] = {}
    if not args.skip_sim:
        for cell_id, lab, _ in CELLS:
            trajs[lab] = run_cell(
                work=work,
                cell_id=cell_id,
                runs=args.runs,
                blocks=args.blocks,
                n_jobs=args.n_jobs,
                seed=args.seed,
                force=args.force,
            )
    else:
        for cell_id, lab, _ in CELLS:
            stem = f"p53_noprune_gin11a_{cell_id}_{STEM_TAG}_b{args.blocks}_n{args.runs}"
            traj = work / "spn_sim" / f"{stem}_trajectory.json"
            if not traj.is_file():
                raise FileNotFoundError(traj)
            trajs[lab] = traj

    traj_series = {
        lab: _load_traj_plateau(trajs[lab], n_blocks=args.plateau_blocks) for lab in labels
    }

    # Layer corr via existing helper (also writes md companion)
    corr_json = work / "corr_by_layer_noprune_gin11a_plateau100.json"
    corr_md = work / "corr_by_layer_noprune_gin11a_plateau100.md"
    cmd = [
        sys.executable,
        "-u",
        str(ROOT / "scripts/compare_spn_corr_by_layer.py"),
        "--out-json",
        str(corr_json),
        "--out-md",
        str(corr_md),
        "--traj-aggregate",
        "plateau_mean",
        "--plateau-blocks",
        str(args.plateau_blocks),
    ]
    for lab in labels:
        cmd += ["--traj", lab, str(trajs[lab])]
        cell_id = next(c for c, l, _ in CELLS if l == lab)
        cmd += ["--levels", lab, str(_levels_path(cell_id))]
    print("[corr]", " ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=ROOT)

    _ensure_colormaps()
    layer_data = {
        layer: {
            lab: {p: v for p, v in traj_series[lab].items() if _place_layer(p) == layer}
            for lab in labels
        }
        for layer in ("histone_mark", "promoter.active")
    }
    # Color scale from actual offdiag range (nice floor → 1.0); dendrogram aligned in plotter.
    panel_specs = (
        ("C", "histone_mark", "wine", "Simulated histones", "fig4_C_gin11a_noprune_histone"),
        ("D", "promoter.active", "vivid_red", "Simulated promoters", "fig4_D_gin11a_noprune_promoter"),
    )
    metas = {}
    for letter, layer, cmap, title, stem_name in panel_specs:
        corr, keys = _pearson_matrix(layer_data[layer], labels)
        stem = fig_dir / stem_name
        meta = plot_cell_corr_heatmap(
            corr,
            labels,
            out_stem=stem,
            title=title,
            cmap=cmap,
            vmin=None,
            vmax=1.0,
            cbar_label="Pearson r",
            cbar_ticks=None,
            linewidths=1.0,
            show_title=True,
        )
        meta.update({"panel": letter, "layer": layer, "n_places": len(keys)})
        metas[letter] = meta
        print(
            f"[plot] {letter} {layer} n={len(keys)} offdiag={meta['offdiag_range']} vmin={meta['vmin']}",
            flush=True,
        )

    summary = {
        "pipeline": "fig4cd_noprune_gin11a_stageA",
        "spn_pipeline_version": SPN_PIPELINE_VERSION,
        "stage_a": str(stage_a.relative_to(ROOT)),
        "work_dir": str(work.relative_to(ROOT)),
        "blocks": args.blocks,
        "runs": args.runs,
        "plateau_blocks": args.plateau_blocks,
        "grn_prune": False,
        "place_levels": f"runtime/omics_seed/cellline_levels/{{cell}}/{{cell}}_{LEVELS_SUFFIX}",
        "levels_note": "Unified file family for all 6 cells (fixes H1 seed-only rna_place_levels asymmetry)",
        "trajs": {k: str(v.relative_to(ROOT)) for k, v in trajs.items()},
        "corr_json": str(corr_json.relative_to(ROOT)),
        "panels": metas,
    }
    (work / "fig4cd_noprune_gin11a_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

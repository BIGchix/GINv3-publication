#!/usr/bin/env python3
"""Fig4 B–D from frozen Stage A (gin11a): prune → std control sims → layer PCC → heatmaps.

Does not overwrite publication H1 pulse dir; writes under:
  gin_out/p53_cell_grn_prune_v4_gin11a/
  gin_out/p53_fig4_gin11a/
  manuscript/figures/fig4/fig4_{B,C,D}_gin11a_*
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from gin_tool.subgraph_pipeline.spn_policy import SPN_PIPELINE_VERSION_STD
from gin_tool.subgraph_pipeline.spn_standard_policy import (
    STD_B0_BLOCKS,
    standard_h1_b0_simulate_kwargs,
)
from gin_tool.subgraph_pipeline.spn_simulate import (
    export_and_simulate_subgraph_dir,
    load_place_levels_tsv,
)

CELLS = (
    ("h1_wa01", "H1", "h1_wa01_E003_promoter_emission.tsv"),
    ("h9_wa09", "H9", "h9_wa09_E008_promoter_emission.tsv"),
    ("huvec", "HUVEC", "huvec_E122_promoter_emission.tsv"),
    ("imr90", "IMR90", "imr90_E017_promoter_emission.tsv"),
    ("keratinocyte", "NHEK", "keratinocyte_E127_promoter_emission.tsv"),
    ("tcell_cd4", "CD4", "tcell_cd4_E043_promoter_emission.tsv"),
)

# Board fig4-2 B–D (SI): Input RNA / Simulated histones / Simulated promoters
BOARD_LAYERS = (
    ("B", "input_RNA", "Input RNA"),
    ("C", "histone_mark", "Simulated histone marks"),
    ("D", "promoter.active", "Simulated promoter.active"),
)


def _place_layer(name: str) -> str:
    low = name.lower()
    if name.endswith(".protein"):
        return "protein"
    if ".promoter.active" in low:
        return "promoter.active"
    if ".promoter." in low:
        mark = name.split(".promoter.", 1)[1]
        if mark.startswith(("H2", "H3", "H4")):
            return "histone_mark"
        return "promoter.other"
    if low.endswith(".rna"):
        return "RNA"
    if name.startswith(("SID", "CID", "CHEBI", "HMDB", "CAS")):
        return "metabolite"
    return "other"


def prune_cells(*, stage_a: Path, out_root: Path, seed_file: Path) -> None:
    cmd = [
        sys.executable,
        "-u",
        str(ROOT / "scripts/build_pathway_cell_grn_prune.py"),
        "--context-edges",
        str(stage_a / "subgraph_edges.tsv"),
        "--seed-file",
        str(seed_file),
        "--out-root",
        str(out_root),
        "--cells",
        ",".join(c[0] for c in CELLS),
    ]
    print("[prune]", " ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=ROOT)
    # copy seed_manifest into each cell dir for SPN export seeds
    man = stage_a / "seed_manifest.json"
    if man.is_file():
        for cell_id, _, _ in CELLS:
            dst = out_root / cell_id / "seed_manifest.json"
            if not dst.is_file():
                shutil.copy2(man, dst)


def run_control_sim(
    *,
    prune_root: Path,
    cell_id: str,
    emw: str,
    runs: int,
    blocks: int,
    n_jobs: int,
    force: bool,
    seed: int,
) -> Path:
    cell_dir = prune_root / cell_id
    levels = ROOT / "runtime/omics_seed/cellline_levels" / cell_id / f"{cell_id}_rna_plus_met_sources_place_levels.tsv"
    emission = ROOT / "runtime/chromhmm/panel_emissions" / emw
    sub = f"spn_sim_fig4bcd_gin11a_ctrl_b{blocks}_n{runs}"
    stem = f"p53_fig4bcd_gin11a_{cell_id}_ctrl_b{blocks}_n{runs}"
    traj = cell_dir / sub / f"{stem}_trajectory.json"
    if traj.is_file() and not force:
        print(f"[sim] skip existing {traj}", flush=True)
        return traj
    for p in (cell_dir / "subgraph_edges.tsv", levels, emission):
        if not p.is_file():
            raise FileNotFoundError(p)
    print(f"[sim] {cell_id} blocks={blocks} runs={runs} jobs={n_jobs}", flush=True)
    kw = standard_h1_b0_simulate_kwargs(
        chromhmm_emission_tsv=str(emission),
        blocks=blocks,
        runs=runs,
        n_jobs=n_jobs,
        seed=seed,
    )
    kw.update(
        place_levels=load_place_levels_tsv(levels),
        sim_subdir=sub,
        stem=stem,
        keep_trajectories=True,
    )
    export_and_simulate_subgraph_dir(cell_dir, **kw)
    if not traj.is_file():
        raise FileNotFoundError(traj)
    return traj


def _load_traj_plateau(path: Path, *, n_blocks: int = 100) -> dict[str, float]:
    d = json.loads(path.read_text(encoding="utf-8"))
    mean = np.asarray(d["mean_markings"], float)
    places = list(d["places"])
    tail = mean[-n_blocks:] if mean.shape[0] >= n_blocks else mean
    return {p: float(tail[:, i].mean()) for i, p in enumerate(places)}


def _load_levels_rna(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            name, raw = parts[0], parts[1]
            if name.lower() in {"place", "name", "node"} or raw.lower() in {"level", "value", "token"}:
                continue
            try:
                val = float(raw)
            except ValueError:
                continue
            if name.endswith(".RNA") or name.lower().endswith(".rna"):
                out[name] = val
    return out


def _pearson_matrix(
    series: dict[str, dict[str, float]], labels: list[str]
) -> tuple[np.ndarray, list[str]]:
    # common places across cells
    keys = set(series[labels[0]])
    for lab in labels[1:]:
        keys &= set(series[lab])
    keys = sorted(keys)
    if len(keys) < 5:
        raise RuntimeError(f"too few shared places: {len(keys)}")
    X = np.array([[series[lab][k] for k in keys] for lab in labels], dtype=float)
    # row-wise standardize then corrcoef on rows
    Xn = X - X.mean(1, keepdims=True)
    denom = np.linalg.norm(Xn, axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    Xn = Xn / denom
    return Xn @ Xn.T, keys


# Board / WeChat fig4 panel style (Aug 2026): white→vivid hue, top dendrogram only,
# colorbar right, ylabels left, xtick 45°, no cell borders. See gin_out/p53_fig4/*_meta.json.
_CMAP_STOPS: dict[str, list[str]] = {
    "vivid_green": ["#ffffff", "#d1ebd2", "#90cd93", "#00c853", "#42a047"],
    "wine": ["#ffffff", "#fce2eb", "#f595b5", "#bf1659", "#9c0048"],
    "vivid_red": ["#ffffff", "#ffe5e8", "#f47c7b", "#e8433f", "#d50000"],
    "vivid_orange": ["#ffffff", "#fff4e4", "#ffce87", "#ff9d0e", "#ff6d00"],
    "vivid_blue": ["#ffffff", "#e3f2fd", "#4baaf5", "#2196f3", "#00b0ff"],
    "vivid_teal": ["#ffffff", "#e0f7fa", "#7edde9", "#00e5ff", "#00acc1"],
    "vivid_purple": ["#ffffff", "#f3e5f5", "#ce93d8", "#8e24aa", "#6a1b9a"],
}


def _ensure_colormaps() -> None:
    """Register fig4 board colormaps sampled from gin_out/p53_fig4 reference panels."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    for name, stops in _CMAP_STOPS.items():
        if name in plt.colormaps():
            continue
        cmap = LinearSegmentedColormap.from_list(name, stops, N=256)
        try:
            mpl.colormaps.register(cmap, name=name)
        except ValueError:
            pass


def _nice_vmin(offdiag_min: float, *, fixed: float | None = None, step: float | None = None) -> float:
    if fixed is not None:
        return float(fixed)
    if not np.isfinite(offdiag_min):
        return 0.0
    # Tight floor near unity; coarser when dynamic range is wide.
    if step is None:
        step = 0.01 if offdiag_min >= 0.9 else 0.05
    return float(max(0.0, np.floor(offdiag_min / step) * step))


def _nice_cbar_ticks(vmin: float, vmax: float, *, max_ticks: int = 6) -> list[float]:
    """Even ticks spanning [vmin, vmax], inclusive ends."""
    span = float(vmax) - float(vmin)
    if not np.isfinite(span) or span <= 0:
        return [float(vmin), float(vmax)]
    # Prefer 0.01 / 0.02 / 0.05 / 0.1 steps
    for step in (0.01, 0.02, 0.025, 0.05, 0.1):
        n = int(round(span / step))
        if 2 <= n <= max_ticks - 1:
            ticks = [float(vmin) + i * step for i in range(n + 1)]
            ticks[-1] = float(vmax)
            return [round(t, 6) for t in ticks]
    return [float(vmin), float(vmax)]


def plot_cell_corr_heatmap(
    corr: np.ndarray,
    labels: list[str],
    *,
    out_stem: Path,
    title: str,
    cmap: str,
    vmin: float | None,
    vmax: float,
    cbar_label: str | None = "Pearson r",
    cbar_ticks: list[float] | None = None,
    linewidths: float = 1.0,
    show_title: bool = True,
) -> dict:
    """Board-style C/D heatmap: data-scaled colorbar; dendrogram locked to heatmap columns."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import squareform

    _ensure_colormaps()
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "axes.linewidth": 0.6,
        }
    )
    n = corr.shape[0]
    dist = 1.0 - np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(dist, 0.0)
    dist = np.maximum(dist, 0.0)
    condensed = squareform(dist, checks=False)
    Z = hierarchy.linkage(condensed, method="average")
    order = hierarchy.leaves_list(Z)
    labels_ord = [labels[i] for i in order]
    corr_ord = corr[np.ix_(order, order)]
    off = corr[~np.eye(n, dtype=bool)]
    off_min = float(off.min()) if off.size else 0.0
    vmin_use = _nice_vmin(off_min, fixed=vmin)
    vmax_use = float(vmax)
    ticks_use = cbar_ticks if cbar_ticks is not None else _nice_cbar_ticks(vmin_use, vmax_use)

    # Heatmap is the anchor; dendrogram + cbar are attached so columns stay aligned.
    fig_h = 5.6 if show_title else 5.2
    fig, ax = plt.subplots(figsize=(5.4, fig_h))
    divider = make_axes_locatable(ax)
    ax_d = divider.append_axes("top", size="14%", pad=0.06)
    cax = divider.append_axes("right", size="5%", pad=0.12)

    # Draw dendrogram in heatmap x-units (0..n): scipy icoord leaf marks 5,15,...,10n-5
    # map to cell centers 0.5,1.5,...,n-0.5 — keeps stems on column centers.
    ddata = hierarchy.dendrogram(Z, no_plot=True, color_threshold=0)
    for xs, ys in zip(ddata["icoord"], ddata["dcoord"]):
        ax_d.plot(
            [x / 10.0 for x in xs],
            ys,
            color="#111111",
            solid_capstyle="round",
            lw=1.1,
        )
    ymax = max(max(ys) for ys in ddata["dcoord"]) if ddata["dcoord"] else 1.0
    ax_d.set_xlim(0.0, float(n))
    ax_d.set_ylim(0.0, ymax * 1.05)
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    for sp in ax_d.spines.values():
        sp.set_visible(False)
    ax_d.set_facecolor("white")
    ax_d.tick_params(length=0)
    ax_d.set_xlabel("")
    ax_d.set_ylabel("")

    sns.heatmap(
        corr_ord,
        ax=ax,
        cmap=cmap,
        vmin=vmin_use,
        vmax=vmax_use,
        xticklabels=labels_ord,
        yticklabels=labels_ord,
        square=True,
        linewidths=float(linewidths),
        linecolor="white",
        cbar_ax=cax,
        cbar_kws={"ticks": ticks_use},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xlim(0.0, float(n))
    ax_d.set_xlim(0.0, float(n))
    ax.tick_params(axis="x", labelrotation=45, labelsize=12, pad=2)
    ax.tick_params(axis="y", labelrotation=0, labelsize=12, pad=2)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontweight("bold")
        lab.set_color("#111111")
    for lab in ax.get_xticklabels():
        lab.set_ha("right")
    cax.tick_params(labelsize=10)
    if cbar_label:
        cax.set_ylabel(cbar_label, fontsize=11, rotation=270, labelpad=14)

    fig.canvas.draw()
    pos = ax.get_position()
    dpos = ax_d.get_position()
    # Re-lock dendrogram width to heatmap after square=True may shrink the parent.
    ax_d.set_position([pos.x0, dpos.y0, pos.width, dpos.height])
    ax_d.set_xlim(ax.get_xlim())

    title_artist = None
    if show_title:
        dpos = ax_d.get_position()
        x_mid = pos.x0 + pos.width / 2.0
        title_artist = fig.text(
            x_mid,
            min(dpos.y1 + 0.02, 0.98),
            title,
            ha="center",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            color="#111111",
        )

    # Leave room for y labels / rotated x labels / title; avoid bbox_inches='tight'
    # which reflows axes and breaks dendrogram↔heatmap column alignment.
    fig.subplots_adjust(left=0.18, right=0.88, bottom=0.16, top=0.88 if show_title else 0.92)
    fig.canvas.draw()
    pos = ax.get_position()
    dpos = ax_d.get_position()
    ax_d.set_position([pos.x0, dpos.y0, pos.width, dpos.height])
    ax_d.set_xlim(ax.get_xlim())
    if title_artist is not None:
        dpos = ax_d.get_position()
        title_artist.set_position((pos.x0 + pos.width / 2.0, min(dpos.y1 + 0.02, 0.98)))

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(
            out_stem.with_suffix(f".{ext}"),
            facecolor="white",
            edgecolor="none",
            bbox_inches=None,
            dpi=300,
        )
    plt.close(fig)
    meta = {
        "title": title,
        "leaf_order": labels_ord,
        "offdiag_range": [float(off.min()), float(off.max())] if off.size else [None, None],
        "metric": "pearson",
        "linkage": "average on 1-r",
        "vmin": vmin_use,
        "vmax": vmax_use,
        "cmap": cmap,
        "cbar_label": cbar_label,
        "cbar_ticks": ticks_use,
        "linewidths": float(linewidths),
        "scale": "data-driven nice floor of offdiag min" if vmin is None else "fixed vmin",
        "layout": "dendrogram in heatmap x-units via make_axes_locatable (column-aligned)",
        "png": str(out_stem.with_suffix(".png").relative_to(ROOT)),
        "pdf": str(out_stem.with_suffix(".pdf").relative_to(ROOT)),
    }
    out_stem.with_name(out_stem.name + "_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stage-a",
        type=Path,
        default=ROOT / "gin_out/p53_shared_context_v4_rh1_gin11a",
    )
    ap.add_argument(
        "--prune-root",
        type=Path,
        default=ROOT / "gin_out/p53_cell_grn_prune_v4_gin11a",
    )
    ap.add_argument(
        "--seed-file",
        type=Path,
        default=ROOT / "seedLists/p53_pi3k_curated/genes.txt",
    )
    ap.add_argument("--out-dir", type=Path, default=ROOT / "gin_out/p53_fig4_gin11a")
    ap.add_argument("--fig-dir", type=Path, default=ROOT / "manuscript/figures/fig4")
    ap.add_argument("--blocks", type=int, default=STD_B0_BLOCKS)
    ap.add_argument("--runs", type=int, default=200)
    ap.add_argument("--n-jobs", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--plateau-blocks", type=int, default=100)
    ap.add_argument("--skip-prune", action="store_true")
    ap.add_argument("--skip-sim", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    stage_a = args.stage_a if args.stage_a.is_absolute() else ROOT / args.stage_a
    prune_root = args.prune_root if args.prune_root.is_absolute() else ROOT / args.prune_root
    out_dir = args.out_dir if args.out_dir.is_absolute() else ROOT / args.out_dir
    fig_dir = args.fig_dir if args.fig_dir.is_absolute() else ROOT / args.fig_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not (stage_a / "subgraph_edges.tsv").is_file():
        raise SystemExit(f"missing Stage A edges: {stage_a}")

    if not args.skip_prune:
        prune_cells(stage_a=stage_a, out_root=prune_root, seed_file=args.seed_file)

    trajs: dict[str, Path] = {}
    levels: dict[str, Path] = {}
    for cell_id, label, emw in CELLS:
        levels[label] = (
            ROOT
            / "runtime/omics_seed/cellline_levels"
            / cell_id
            / f"{cell_id}_rna_plus_met_sources_place_levels.tsv"
        )
        if args.skip_sim:
            # find existing
            matches = sorted((prune_root / cell_id).glob("spn_sim_fig4bcd_gin11a_ctrl_*/p53_fig4bcd_gin11a_*_trajectory.json"))
            if not matches:
                raise SystemExit(f"no traj for {cell_id}; run without --skip-sim")
            trajs[label] = matches[-1]
        else:
            trajs[label] = run_control_sim(
                prune_root=prune_root,
                cell_id=cell_id,
                emw=emw,
                runs=args.runs,
                blocks=args.blocks,
                n_jobs=args.n_jobs,
                force=args.force,
                seed=args.seed,
            )

    labels = [lab for _, lab, _ in CELLS]
    traj_series = {lab: _load_traj_plateau(trajs[lab], n_blocks=args.plateau_blocks) for lab in labels}
    level_series = {lab: _load_levels_rna(levels[lab]) for lab in labels}

    # Full layer summary via compare_spn_corr_by_layer
    corr_json = out_dir / "corr_by_layer_gin11a_plateau100.json"
    corr_md = out_dir / "corr_by_layer_gin11a_plateau100.md"
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
        cmd += ["--levels", lab, str(levels[lab])]
    print("[corr]", " ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=ROOT)

    # Board heatmaps B–D + extras
    layer_data: dict[str, dict[str, dict[str, float]]] = {"input_RNA": level_series}
    for layer in ("RNA", "protein", "promoter.active", "histone_mark", "metabolite"):
        layer_data[layer] = {
            lab: {p: v for p, v in traj_series[lab].items() if _place_layer(p) == layer}
            for lab in labels
        }

    # Board colors (gin_out/p53_fig4): B green 0.2–1; C wine; D red.
    # Fixed vmin only where historical board locked it; else nice-floor offdiag min.
    style = {
        "input_RNA": ("vivid_green", 0.20, 1.0),
        "histone_mark": ("wine", None, 1.0),
        "promoter.active": ("vivid_red", None, 1.0),
        "RNA": ("vivid_teal", None, 1.0),
        "protein": ("vivid_orange", None, 1.0),
        "metabolite": ("vivid_blue", None, 1.0),
    }
    panel_map = {
        "B": ("input_RNA", "Fig. 4B  Input RNA (ENCODE place levels)"),
        "C": ("histone_mark", "Fig. 4C  Simulated histone marks"),
        "D": ("promoter.active", "Fig. 4D  Simulated promoter.active"),
    }
    metas = {}
    for letter, (layer, title) in panel_map.items():
        cmap, vmin, vmax = style[layer]
        corr, keys = _pearson_matrix(layer_data[layer], labels)
        stem = fig_dir / f"fig4_{letter}_gin11a_{layer}"
        meta = plot_cell_corr_heatmap(
            corr,
            labels,
            out_stem=stem,
            title=title,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        meta["n_places"] = len(keys)
        meta["layer"] = layer
        meta["panel"] = letter
        metas[letter] = meta
        print(f"[plot] {letter} {layer} offdiag={meta['offdiag_range']} n={len(keys)}", flush=True)

    # Also dump companion layers under gin_out
    for layer in ("RNA", "protein", "metabolite"):
        cmap, vmin, vmax = style[layer]
        corr, keys = _pearson_matrix(layer_data[layer], labels)
        stem = out_dir / f"fig4_companion_gin11a_{layer}"
        plot_cell_corr_heatmap(
            corr, labels, out_stem=stem, title=f"Simulated {layer}", cmap=cmap, vmin=vmin, vmax=vmax
        )

    summary = {
        "pipeline": SPN_PIPELINE_VERSION_STD,
        "stage_a": str(stage_a.relative_to(ROOT)),
        "prune_root": str(prune_root.relative_to(ROOT)),
        "blocks": args.blocks,
        "runs": args.runs,
        "plateau_blocks": args.plateau_blocks,
        "trajs": {k: str(v.relative_to(ROOT)) for k, v in trajs.items()},
        "corr_json": str(corr_json.relative_to(ROOT)),
        "panels": metas,
        "board_note": "fig4-2 B–D = Input RNA / histone / promoter.active (SI)",
    }
    (out_dir / "fig4bcd_gin11a_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

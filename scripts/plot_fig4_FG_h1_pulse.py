#!/usr/bin/env python3
"""Figure 4 panels F–G: H1 standard B0 pulse trajectories (publication style).

F: TP53.protein (control vs step), blocks 100–350, peak annotated.
G: MDM2 promoter.active / RNA / protein (absolute tokens), small multiples.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "runtime/spn_phys5/p53_B0_v15cfea10e3ewe_o150b350_metsrc_h40_report.json"
OUT_DIR = ROOT / "manuscript/figures/fig4"

ONSET = 150
HOLD = 40
X0, X1 = 100, 350

# Colorblind-friendly (Okabe–Ito inspired)
C_TP53 = "#0072B2"
C_PROM = "#E69F00"
C_RNA = "#009E73"
C_PROT = "#CC79A7"
C_CTRL = "#666666"


def _setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _load_arm_series(traj_path: Path, places: list[str]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    d = json.loads(traj_path.read_text(encoding="utf-8"))
    all_places = list(d["places"])
    mean = np.asarray(d["mean_markings"], dtype=float)
    std = np.asarray(d["std_markings"], dtype=float)
    n = int(d["n_runs"])
    sem = std / np.sqrt(max(n, 1))
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name in places:
        j = all_places.index(name)
        out[name] = (mean[:, j].copy(), sem[:, j].copy())
    return out


def _blocks(n: int) -> np.ndarray:
    return np.arange(n, dtype=int)


def _mask(blocks: np.ndarray) -> np.ndarray:
    return (blocks >= X0) & (blocks <= X1)


def _peak_in_window(
    blocks: np.ndarray, y: np.ndarray, *, start: int, end: int
) -> tuple[int, float]:
    idx = np.where((blocks >= start) & (blocks <= end))[0]
    i = int(idx[int(np.argmax(y[idx]))])
    return int(blocks[i]), float(y[i])


def _shade_pulse(ax: plt.Axes) -> None:
    ax.axvspan(ONSET, ONSET + HOLD, color="#B0B0B0", alpha=0.22, lw=0, zorder=0)
    ymin, ymax = ax.get_ylim()
    ax.text(
        ONSET + HOLD / 2,
        ymax - 0.02 * (ymax - ymin),
        "pulse",
        ha="center",
        va="top",
        fontsize=6.5,
        color="#555555",
        style="italic",
        zorder=1,
    )


def _annotate_peak(
    ax: plt.Axes,
    *,
    block: int,
    token: float,
    color: str,
    xytext: tuple[float, float],
) -> None:
    ax.plot(block, token, "o", ms=4.5, color=color, zorder=5, markeredgecolor="white", markeredgewidth=0.6)
    ax.annotate(
        f"{token:.1f} @ {block}",
        xy=(block, token),
        xytext=xytext,
        textcoords="offset points",
        fontsize=6.5,
        color=color,
        ha="left",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=0, shrinkB=2),
        zorder=6,
    )


def _plot_mean_band(
    ax: plt.Axes,
    blocks: np.ndarray,
    mean: np.ndarray,
    sem: np.ndarray,
    *,
    color: str,
    ls: str,
    label: str,
    lw: float = 1.4,
) -> None:
    m = _mask(blocks)
    x = blocks[m]
    y = mean[m]
    e = sem[m]
    ax.fill_between(x, y - e, y + e, color=color, alpha=0.18, lw=0, zorder=2)
    ax.plot(x, y, ls=ls, color=color, lw=lw, label=label, zorder=3)


def plot_panel_f(
    ctrl: dict[str, tuple[np.ndarray, np.ndarray]],
    step: dict[str, tuple[np.ndarray, np.ndarray]],
    blocks: np.ndarray,
    out: Path,
) -> Path:
    fig, ax = plt.subplots(figsize=(3.2, 2.4))
    c_m, c_e = ctrl["TP53.protein"]
    s_m, s_e = step["TP53.protein"]
    _plot_mean_band(ax, blocks, c_m, c_e, color=C_CTRL, ls="--", label="control", lw=1.2)
    _plot_mean_band(ax, blocks, s_m, s_e, color=C_TP53, ls="-", label="pulse", lw=1.5)

    ax.set_xlim(X0, X1)
    ax.set_xlabel("Block")
    ax.set_ylabel("TP53 protein (tokens)")
    ax.set_title("F  Inhibition arm", loc="left", fontweight="bold", pad=6)
    _shade_pulse(ax)

    pk_b, pk_v = _peak_in_window(blocks, s_m, start=ONSET, end=ONSET + HOLD + 20)
    # place annotation to the right of peak if room
    xoff = 12 if pk_b < 280 else -70
    yoff = 8
    _annotate_peak(ax, block=pk_b, token=pk_v, color=C_TP53, xytext=(xoff, yoff))

    ax.legend(frameon=False, loc="upper right", handlelength=2.2)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out.with_suffix(f".{ext}"))
    plt.close(fig)
    return out


def plot_panel_g(
    ctrl: dict[str, tuple[np.ndarray, np.ndarray]],
    step: dict[str, tuple[np.ndarray, np.ndarray]],
    blocks: np.ndarray,
    out: Path,
) -> Path:
    layers = [
        ("MDM2.promoter.active", "MDM2 promoter.active", C_PROM),
        ("MDM2.RNA", "MDM2 RNA", C_RNA),
        ("MDM2.protein", "MDM2 protein", C_PROT),
    ]
    fig, axes = plt.subplots(
        3, 1, figsize=(3.2, 4.6), sharex=True, gridspec_kw={"hspace": 0.12}
    )
    ann_offsets = [(14, 6), (14, 6), (14, 6)]

    for ax, (place, title, color), xytext in zip(axes, layers, ann_offsets):
        c_m, c_e = ctrl[place]
        s_m, s_e = step[place]
        _plot_mean_band(ax, blocks, c_m, c_e, color=C_CTRL, ls="--", label="control", lw=1.1)
        _plot_mean_band(ax, blocks, s_m, s_e, color=color, ls="-", label="pulse", lw=1.4)
        ax.set_xlim(X0, X1)
        ax.set_ylabel("Tokens")
        ax.text(0.02, 0.92, title, transform=ax.transAxes, ha="left", va="top", fontsize=7.5, fontweight="semibold")
        _shade_pulse(ax)
        pk_b, pk_v = _peak_in_window(blocks, s_m, start=ONSET, end=ONSET + HOLD + 80)
        xoff = 12 if pk_b < 280 else -75
        _annotate_peak(ax, block=pk_b, token=pk_v, color=color, xytext=(xoff, xytext[1]))

    axes[0].set_title("G  Induction arm (MDM2 layers)", loc="left", fontweight="bold", pad=6)
    axes[-1].set_xlabel("Block")
    axes[0].legend(frameon=False, loc="upper right", ncol=2, handlelength=2.0)

    fig.subplots_adjust(hspace=0.12, left=0.16, right=0.98, top=0.93, bottom=0.08)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out.with_suffix(f".{ext}"))
    plt.close(fig)
    return out


def plot_fg_row(
    ctrl: dict[str, tuple[np.ndarray, np.ndarray]],
    step: dict[str, tuple[np.ndarray, np.ndarray]],
    blocks: np.ndarray,
    out: Path,
) -> Path:
    """Optional side-by-side F|G composite for layout preview."""
    fig = plt.figure(figsize=(7.0, 4.6))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.05], wspace=0.32, hspace=0.14)

    ax_f = fig.add_subplot(gs[:, 0])
    c_m, c_e = ctrl["TP53.protein"]
    s_m, s_e = step["TP53.protein"]
    _plot_mean_band(ax_f, blocks, c_m, c_e, color=C_CTRL, ls="--", label="control", lw=1.2)
    _plot_mean_band(ax_f, blocks, s_m, s_e, color=C_TP53, ls="-", label="pulse", lw=1.5)
    ax_f.set_xlim(X0, X1)
    ax_f.set_xlabel("Block")
    ax_f.set_ylabel("TP53 protein (tokens)")
    ax_f.set_title("F  Inhibition arm", loc="left", fontweight="bold", pad=6)
    _shade_pulse(ax_f)
    pk_b, pk_v = _peak_in_window(blocks, s_m, start=ONSET, end=ONSET + HOLD + 20)
    _annotate_peak(ax_f, block=pk_b, token=pk_v, color=C_TP53, xytext=(12, 8))
    ax_f.legend(frameon=False, loc="upper right")

    layers = [
        ("MDM2.promoter.active", "MDM2 promoter.active", C_PROM),
        ("MDM2.RNA", "MDM2 RNA", C_RNA),
        ("MDM2.protein", "MDM2 protein", C_PROT),
    ]
    for i, (place, title, color) in enumerate(layers):
        ax = fig.add_subplot(gs[i, 1], sharex=ax_f if i else None)
        c_m, c_e = ctrl[place]
        s_m, s_e = step[place]
        _plot_mean_band(ax, blocks, c_m, c_e, color=C_CTRL, ls="--", label="control", lw=1.1)
        _plot_mean_band(ax, blocks, s_m, s_e, color=color, ls="-", label="pulse", lw=1.4)
        ax.set_xlim(X0, X1)
        ax.set_ylabel("Tokens")
        ax.text(0.02, 0.92, title, transform=ax.transAxes, ha="left", va="top", fontsize=7.5)
        _shade_pulse(ax)
        pk_b, pk_v = _peak_in_window(blocks, s_m, start=ONSET, end=ONSET + HOLD + 80)
        _annotate_peak(ax, block=pk_b, token=pk_v, color=color, xytext=(12, 6))
        if i == 0:
            ax.set_title("G  Induction arm (MDM2 layers)", loc="left", fontweight="bold", pad=6)
            ax.legend(frameon=False, loc="upper right", ncol=2, handlelength=2.0)
        if i == 2:
            ax.set_xlabel("Block")
        else:
            ax.tick_params(labelbottom=False)

    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", type=Path, default=REPORT)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    _setup_style()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    ctrl_path = Path(report["arms"]["control"]["trajectory_json"])
    step_path = Path(report["arms"]["step"]["trajectory_json"])
    places = [
        "TP53.protein",
        "MDM2.promoter.active",
        "MDM2.RNA",
        "MDM2.protein",
    ]
    ctrl = _load_arm_series(ctrl_path, places)
    step = _load_arm_series(step_path, places)
    n = next(iter(ctrl.values()))[0].shape[0]
    blocks = _blocks(n)

    out_dir = args.out_dir
    meta = {
        "source_report": str(args.report),
        "control_traj": str(ctrl_path),
        "step_traj": str(step_path),
        "onset": ONSET,
        "hold": HOLD,
        "x_range": [X0, X1],
        "n_runs": report.get("runs"),
        "peaks_step": {},
    }
    for place in places:
        b, v = _peak_in_window(blocks, step[place][0], start=ONSET, end=ONSET + HOLD + 80)
        meta["peaks_step"][place] = {"block": b, "tokens": round(v, 3)}

    f_path = plot_panel_f(ctrl, step, blocks, out_dir / "fig4_F_TP53_protein")
    g_path = plot_panel_g(ctrl, step, blocks, out_dir / "fig4_G_MDM2_layers")
    row_path = plot_fg_row(ctrl, step, blocks, out_dir / "fig4_FG_row")
    (out_dir / "fig4_FG_meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {f_path.with_suffix('.pdf')}")
    print(f"wrote {g_path.with_suffix('.pdf')}")
    print(f"wrote {row_path.with_suffix('.pdf')}")
    print(json.dumps(meta["peaks_step"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

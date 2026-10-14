#!/usr/bin/env python3
"""Figure 4F–I: 2×2 pulse panels (publication style, matching fig4-2.pdf).

F TP53.protein · G MDM2.promoter.active · H MDM2.RNA · I MDM2.protein
Style aligned with scripts/plot_fig4_FG_h1_pulse.py / manuscript fig4-2.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

ONSET = 150
HOLD = 40
X0, X1 = 100, 350

C_TP53 = "#0072B2"
C_PROM = "#E69F00"
C_RNA = "#009E73"
C_PROT = "#CC79A7"
C_CTRL = "#666666"

PANELS = [
    ("F", "TP53.protein", "TP53 protein", C_TP53),
    ("G", "MDM2.promoter.active", "MDM2 promoter", C_PROM),
    ("H", "MDM2.RNA", "MDM2 RNA", C_RNA),
    ("I", "MDM2.protein", "MDM2 protein", C_PROT),
]


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
        "Pulse",
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


def plot_fhi_grid(
    ctrl: dict[str, tuple[np.ndarray, np.ndarray]],
    step: dict[str, tuple[np.ndarray, np.ndarray]],
    blocks: np.ndarray,
    out: Path,
    *,
    meta_extra: dict | None = None,
) -> dict:
    _setup_style()
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.8), sharex=True)
    axes_flat = axes.ravel()
    peaks: dict[str, dict] = {}

    for ax, (letter, place, title, color) in zip(axes_flat, PANELS):
        c_m, c_e = ctrl[place]
        s_m, s_e = step[place]
        _plot_mean_band(ax, blocks, c_m, c_e, color=C_CTRL, ls="--", label="Control", lw=1.15)
        _plot_mean_band(ax, blocks, s_m, s_e, color=color, ls="-", label="Pulse", lw=1.5)
        ax.set_xlim(X0, X1)
        ax.set_ylabel("Tokens")
        ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", pad=4)
        end = ONSET + HOLD + (20 if place == "TP53.protein" else 80)
        pk_b, pk_v = _peak_in_window(blocks, s_m, start=ONSET, end=end)
        peaks[place] = {"block": pk_b, "tokens": round(pk_v, 3)}
        # set ylim before shade so Pulse label sits correctly
        m = _mask(blocks)
        ymax = float(max((s_m[m] + s_e[m]).max(), (c_m[m] + c_e[m]).max()))
        ymin = float(min((s_m[m] - s_e[m]).min(), (c_m[m] - c_e[m]).min()))
        pad = 0.08 * max(ymax - ymin, 1.0)
        ax.set_ylim(max(0.0, ymin - pad), ymax + pad)
        _shade_pulse(ax)
        xoff = 10 if pk_b < 280 else -72
        _annotate_peak(ax, block=pk_b, token=pk_v, color=color, xytext=(xoff, 7))
        ax.legend(frameon=False, loc="upper right", handlelength=2.0)

    for ax in axes[1, :]:
        ax.set_xlabel("Pseudo-time block")

    fig.tight_layout(w_pad=1.2, h_pad=1.0)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out.with_suffix(f".{ext}"))
    plt.close(fig)

    meta = {
        "panels": "F–I",
        "layout": "2x2",
        "onset": ONSET,
        "hold": HOLD,
        "x_range": [X0, X1],
        "peaks_step": peaks,
        "png": str(out.with_suffix(".png").relative_to(ROOT)) if out.is_absolute() else str(out.with_suffix(".png")),
        "pdf": str(out.with_suffix(".pdf").relative_to(ROOT)) if out.is_absolute() else str(out.with_suffix(".pdf")),
    }
    if meta_extra:
        meta.update(meta_extra)
    meta_path = out.with_name(out.name + "_meta.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--report",
        type=Path,
        default=ROOT / "runtime/spn_phys5/p53_B0_std_stageA_gin11a_s48_h40_report.json",
    )
    ap.add_argument(
        "--out-stem",
        type=Path,
        default=ROOT / "manuscript/figures/fig4/fig4_FHI_gin11a_s48",
    )
    args = ap.parse_args()

    report = json.loads(args.report.read_text(encoding="utf-8"))
    ctrl_path = Path(report["arms"]["control"]["trajectory_json"])
    step_path = Path(report["arms"]["step"]["trajectory_json"])
    places = [p for _, p, _, _ in PANELS]
    ctrl = _load_arm_series(ctrl_path, places)
    step = _load_arm_series(step_path, places)
    n = next(iter(ctrl.values()))[0].shape[0]
    blocks = np.arange(n, dtype=int)

    out = args.out_stem if args.out_stem.is_absolute() else ROOT / args.out_stem
    meta = plot_fhi_grid(
        ctrl,
        step,
        blocks,
        out,
        meta_extra={
            "source_report": str(args.report),
            "control_traj": str(ctrl_path),
            "step_traj": str(step_path),
            "seed": report.get("seed"),
            "tag": report.get("tag"),
            "n_runs": report.get("runs"),
            "follow_vs_control": report.get("follow_vs_control"),
        },
    )
    print(json.dumps(meta, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

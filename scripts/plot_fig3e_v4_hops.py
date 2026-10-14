#!/usr/bin/env python3
"""Figure 3E: context subgraph size vs real-hop (v4, no GRN prune / no chromHMM).

Grouped bars: real nodes + total nodes (L) and total edges (R twin axis).
Optional pie row: node-type composition (Protein/RNA/Promoter/Histone/…).
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

C_REAL = "#2A6F6F"
C_TOTAL = "#7A6A5A"
C_EDGES = "#C45C26"
KIND_ORDER = ("Protein", "RNA", "Promoter", "Histone", "Metabolite", "Complex", "Virtual")
KIND_COLORS = {
    "Protein": "#1B4F72",
    "RNA": "#1E8449",
    "Promoter": "#5D6D7E",
    "Histone": "#2874A6",
    "Metabolite": "#B7950B",
    "Complex": "#6C3483",
    "Virtual": "#AAB7B8",
}


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
        }
    )


def _classify_kinds(edges_path: Path) -> dict[str, int]:
    from gin_tool.mapping import is_virtual_gin_node
    from gin_tool.spn_export import SpnPlaceKind, classify_place

    nodes: set[str] = set()
    with edges_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            nodes.add(parts[0])
            nodes.add(parts[1])
    kinds: Counter[str] = Counter()
    for n in nodes:
        if is_virtual_gin_node(n) or ";" in n:
            kinds["Virtual"] += 1
            continue
        k = classify_place(n)
        if k == SpnPlaceKind.PROTEIN:
            kinds["Protein"] += 1
        elif k == SpnPlaceKind.RNA:
            kinds["RNA"] += 1
        elif k in (
            SpnPlaceKind.PROMOTER,
            SpnPlaceKind.CHROMATIN_LOCUS,
            SpnPlaceKind.DNA_REGION,
        ):
            kinds["Promoter"] += 1
        elif k == SpnPlaceKind.METABOLITE:
            kinds["Metabolite"] += 1
        elif k == SpnPlaceKind.PROTEIN_COMPLEX:
            kinds["Complex"] += 1
        elif k == SpnPlaceKind.HISTONE_MARK:
            kinds["Histone"] += 1
        else:
            kinds["Protein"] += 1
    return dict(kinds)


def _load_rows(stats_json: Path, roots: list[Path] | None) -> list[dict]:
    payload = json.loads(stats_json.read_text(encoding="utf-8"))
    rows = list(payload["rows"])
    rows.sort(key=lambda r: r["real_hop"])
    if roots is not None:
        by_rh = {r["real_hop"]: r for r in rows}
        rows = []
        for root in roots:
            # prefer matching by path suffix
            hit = next((r for r in payload["rows"] if Path(r["root"]).name == root.name), None)
            if hit is None:
                raise SystemExit(f"no stats row for {root}")
            rows.append(hit)
        rows.sort(key=lambda r: r["real_hop"])
        _ = by_rh
    for r in rows:
        root = Path(r["root"])
        if not root.is_absolute():
            root = ROOT / root
        edges = root / "subgraph_edges.tsv"
        r["kinds"] = _classify_kinds(edges) if edges.is_file() else {}
        r["root_path"] = str(root)
    return rows


def _fmt(n: int) -> str:
    return f"{n:,}"


def plot_fig3e(rows: list[dict], out_stem: Path, *, with_pies: bool) -> dict:
    _setup_style()
    x = np.arange(len(rows))
    width = 0.26
    real = [int(r["n_real_nodes"]) for r in rows]
    total = [int(r["n_total_nodes"]) for r in rows]
    edges = [int(r["n_edges"]) for r in rows]
    labels = [
        f"Real-hop={r['real_hop']}\n(Graph-hop={r['graph_hop']})" for r in rows
    ]

    if with_pies:
        fig = plt.figure(figsize=(7.2, 5.6))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.45, wspace=0.25)
        ax = fig.add_subplot(gs[0, :])
        pie_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    else:
        fig, ax = plt.subplots(figsize=(5.2, 3.2))
        pie_axes = []

    ax2 = ax.twinx()
    b1 = ax.bar(x - width, real, width, color=C_REAL, label="Real nodes", zorder=3)
    b2 = ax.bar(x, total, width, color=C_TOTAL, label="Total nodes", zorder=3)
    b3 = ax2.bar(x + width, edges, width, color=C_EDGES, label="Total edges", zorder=3)

    for bars, vals, color in (
        (b1, real, C_REAL),
        (b2, total, C_TOTAL),
        (b3, edges, C_EDGES),
    ):
        for rect, v in zip(bars, vals):
            ax_ref = ax2 if bars is b3 else ax
            ax_ref.text(
                rect.get_x() + rect.get_width() / 2,
                rect.get_height(),
                _fmt(v),
                ha="center",
                va="bottom",
                fontsize=6.5,
                color=color,
                fontweight="bold",
                clip_on=False,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Number of nodes")
    ax2.set_ylabel("Number of edges", color=C_EDGES)
    ax2.tick_params(axis="y", colors=C_EDGES)
    ax2.spines["right"].set_color(C_EDGES)
    ax.set_ylim(0, max(total) * 1.22)
    ax2.set_ylim(0, max(edges) * 1.22)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax2.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)

    # default depth marker (real-hop=2)
    if any(r["real_hop"] == 2 for r in rows):
        i2 = next(i for i, r in enumerate(rows) if r["real_hop"] == 2)
        i1 = next((i for i, r in enumerate(rows) if r["real_hop"] == 1), None)
        if i1 is not None:
            ax.plot(
                [i1 - width, i2 - width],
                [real[i1], real[i2]],
                color=C_REAL,
                lw=1.2,
                zorder=4,
            )
            ax.annotate(
                "default",
                xy=(i2 - width, real[i2]),
                xytext=(8, 10),
                textcoords="offset points",
                fontsize=7,
                color=C_REAL,
                fontweight="bold",
            )

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=False, ncol=3)

    ver = rows[0].get("context_pipeline_version", "context_v4")
    ax.set_title(
        "Context subgraph size vs real-hop\n"
        f"KEGG p53 · RWR-C TFs · {ver}\n"
        "no GRN allowlist prune · no chromHMM",
        loc="left",
        fontsize=8,
        pad=8,
    )

    if with_pies:
        for axp, r in zip(pie_axes, rows):
            kinds = r.get("kinds") or {}
            sizes = [int(kinds.get(k, 0)) for k in KIND_ORDER]
            colors = [KIND_COLORS[k] for k in KIND_ORDER]
            # drop zero slices for cleaner pies
            nz = [(s, c, k) for s, c, k in zip(sizes, colors, KIND_ORDER) if s > 0]
            if not nz:
                axp.axis("off")
                continue
            sizes_nz, colors_nz, keys_nz = zip(*nz)
            axp.pie(
                sizes_nz,
                colors=colors_nz,
                startangle=90,
                wedgeprops={"linewidth": 0.4, "edgecolor": "white"},
            )
            axp.set_title(
                f"real-hop = {r['real_hop']}  (n = {_fmt(int(r['n_total_nodes']))})",
                fontsize=8,
            )
            # compact counts under pie
            bits = [f"{k[0] if k != 'Virtual' else 'V'}={kinds.get(k, 0)}" for k in KIND_ORDER if kinds.get(k, 0)]
            axp.text(
                0.5,
                -0.12,
                "  ".join(bits),
                transform=axp.transAxes,
                ha="center",
                va="top",
                fontsize=5.5,
                color="#444444",
            )
        # shared legend
        handles = [
            mpl.patches.Patch(facecolor=KIND_COLORS[k], edgecolor="none", label=k)
            for k in KIND_ORDER
        ]
        fig.legend(
            handles=handles,
            loc="lower center",
            ncol=7,
            frameon=False,
            bbox_to_anchor=(0.5, 0.01),
            fontsize=7,
        )

    out_stem = out_stem if out_stem.is_absolute() else (ROOT / out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    png = out_stem.with_suffix(".png")
    pdf = out_stem.with_suffix(".pdf")
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)

    def _rel(p: Path) -> str:
        try:
            return str(p.resolve().relative_to(ROOT))
        except ValueError:
            return str(p)

    meta = {
        "panel": "fig3E",
        "pipeline": ver,
        "rows": [
            {
                "real_hop": r["real_hop"],
                "graph_hop": r["graph_hop"],
                "n_real_nodes": r["n_real_nodes"],
                "n_total_nodes": r["n_total_nodes"],
                "n_edges": r["n_edges"],
                "kinds": r.get("kinds"),
            }
            for r in rows
        ],
        "png": _rel(png),
        "pdf": _rel(pdf),
    }
    meta_path = out_stem.with_name(out_stem.name + "_meta.json")
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stats",
        type=Path,
        default=ROOT / "benchmarks/subgraph-speed/fig3e_v4_context_hop_stats.json",
    )
    ap.add_argument(
        "--out-stem",
        type=Path,
        default=ROOT / "manuscript/figures/fig3/fig3_E_v4_hops",
    )
    ap.add_argument("--no-pies", action="store_true", help="Bars only (manuscript E strip)")
    args = ap.parse_args()
    rows = _load_rows(args.stats, None)
    meta = plot_fig3e(rows, args.out_stem, with_pies=not args.no_pies)
    print(json.dumps(meta, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

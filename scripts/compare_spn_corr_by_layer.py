#!/usr/bin/env python3
"""Layer-stratified Pearson/Spearman across cell-line SPN trajectories."""

from __future__ import annotations

import argparse
import json
import math
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_spn_corr_matrix import (
    _corr_matrix,
    _load_levels,
    _load_traj_final,
    _load_traj_plateau_mean,
    _pearson,
    _spearman,
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


# Ordinal biological distance (0=same lineage; higher=more divergent).
_BIO_DISTANCE: dict[frozenset[str], float] = {
    frozenset({"H1", "H9"}): 0.0,
    frozenset({"HUVEC", "IMR90"}): 1.0,
    frozenset({"IMR90", "NHEK"}): 1.0,
    frozenset({"H1", "HUVEC"}): 2.0,
    frozenset({"H1", "IMR90"}): 2.0,
    frozenset({"H1", "NHEK"}): 2.0,
    frozenset({"H9", "HUVEC"}): 2.0,
    frozenset({"H9", "IMR90"}): 2.0,
    frozenset({"H9", "NHEK"}): 2.0,
    frozenset({"HUVEC", "NHEK"}): 2.0,
    frozenset({"H1", "CD4"}): 3.0,
    frozenset({"H9", "CD4"}): 3.0,
    frozenset({"HUVEC", "CD4"}): 3.0,
    frozenset({"IMR90", "CD4"}): 3.0,
    frozenset({"NHEK", "CD4"}): 3.0,
}

_BIO_TYPE = {
    "H1": "ESC",
    "H9": "ESC",
    "HUVEC": "Endothelial",
    "IMR90": "Fibroblast",
    "NHEK": "Keratinocyte",
    "CD4": "T_cell",
}


def _bio_distance(a: str, b: str) -> float:
    if a == b:
        return 0.0
    return _BIO_DISTANCE.get(frozenset({a, b}), 2.5)


def _subset_series(
    series: dict[str, dict[str, float]], layer: str
) -> tuple[dict[str, dict[str, float]], int]:
    out: dict[str, dict[str, float]] = {}
    n_places = 0
    for lab, full in series.items():
        sub = {p: v for p, v in full.items() if _place_layer(p) == layer}
        if sub:
            n_places = max(n_places, len(sub))
        out[lab] = sub
    return out, n_places


def _pairwise_from_matrix(labels: list[str], mat: dict[str, dict[str, float]]) -> dict[frozenset[str], float]:
    out: dict[frozenset[str], float] = {}
    for a, b in combinations(labels, 2):
        out[frozenset({a, b})] = float(mat[a][b])
    return out


def _corr_bio_vs_metric(pairs: dict[frozenset[str], float], metric: str) -> float | None:
    """Spearman between biological distance and (1 - r) or r."""
    xs: list[float] = []
    ys: list[float] = []
    for key, r in pairs.items():
        if not math.isfinite(r):
            continue
        a, b = tuple(key)
        xs.append(_bio_distance(a, b))
        ys.append(1.0 - r if metric == "dissimilarity" else r)
    if len(xs) < 3:
        return None
    return _spearman(xs, ys)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", action="append", nargs=2, metavar=("LABEL", "JSON"))
    ap.add_argument("--levels", action="append", nargs=2, metavar=("LABEL", "TSV"))
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-md", default="")
    ap.add_argument("--traj-aggregate", choices=("final", "plateau_mean"), default="plateau_mean")
    ap.add_argument("--plateau-blocks", type=int, default=100)
    ap.add_argument("--min-places", type=int, default=5)
    args = ap.parse_args()

    loader = (
        (lambda p: _load_traj_plateau_mean(p, n_blocks=args.plateau_blocks))
        if args.traj_aggregate == "plateau_mean"
        else _load_traj_final
    )
    traj_series = {lab: loader(Path(p)) for lab, p in (args.traj or [])}
    level_series = {lab: _load_levels(Path(p)) for lab, p in (args.levels or [])}

    layers = sorted({_place_layer(p) for full in traj_series.values() for p in full})
    labels = list(traj_series)

    input_rna = _corr_matrix(level_series, union=False) if level_series else None
    input_pairs = (
        _pairwise_from_matrix(labels, input_rna["pearson"]) if input_rna else {}
    )

    report: dict = {
        "meta": {
            "traj_aggregate": args.traj_aggregate,
            "plateau_blocks": args.plateau_blocks,
            "labels": labels,
            "bio_type": _BIO_TYPE,
        },
        "by_layer": {},
    }

    lines: list[str] = [
        "# 分层相关 · 6生理系 SPN 平台期",
        "",
        f"- 聚合：{args.traj_aggregate}（末 {args.plateau_blocks} block）",
        "- 层次：protein / histone_mark / metabolite / promoter.active / RNA / other",
        "",
    ]

    layer_rows: list[tuple[str, float | None, float | None, float, float, int]] = []

    for layer in layers:
        sub, n_places = _subset_series(traj_series, layer)
        if n_places < args.min_places:
            continue
        block = _corr_matrix(sub, union=False)
        pearson_pairs = _pairwise_from_matrix(labels, block["pearson"])
        spearman_pairs = _pairwise_from_matrix(labels, block["spearman"])

        bio_vs_pearson = _corr_bio_vs_metric(pearson_pairs, "dissimilarity")
        bio_vs_spear = _corr_bio_vs_metric(spearman_pairs, "dissimilarity")
        rna_vs_pearson = (
            _spearman(
                [input_pairs[k] for k in pearson_pairs if k in input_pairs],
                [pearson_pairs[k] for k in pearson_pairs if k in input_pairs],
            )
            if input_pairs and len(pearson_pairs) >= 3
            else None
        )

        vals_p = [v for v in pearson_pairs.values() if math.isfinite(v)]
        vals_s = [v for v in spearman_pairs.values() if math.isfinite(v)]
        layer_rows.append(
            (
                layer,
                bio_vs_pearson,
                rna_vs_pearson,
                min(vals_p) if vals_p else float("nan"),
                max(vals_p) if vals_p else float("nan"),
                n_places,
            )
        )

        pair_detail = []
        for key in sorted(pearson_pairs, key=lambda k: pearson_pairs[k]):
            a, b = sorted(key)
            pair_detail.append(
                {
                    "a": a,
                    "b": b,
                    "bio_distance": _bio_distance(a, b),
                    "types": f"{_BIO_TYPE[a]} vs {_BIO_TYPE[b]}",
                    "pearson": pearson_pairs[key],
                    "spearman": spearman_pairs[key],
                    "input_rna_pearson": input_pairs.get(key),
                }
            )

        report["by_layer"][layer] = {
            "n_places": n_places,
            "matrix": block,
            "bio_distance_vs_pearson_dissimilarity_spearman": bio_vs_pearson,
            "bio_distance_vs_spearman_dissimilarity_spearman": bio_vs_spear,
            "input_rna_vs_layer_pearson_spearman": rna_vs_pearson,
            "pairs_sorted_by_pearson": pair_detail,
        }

    report["summary"] = [
        {
            "layer": row[0],
            "n_places": row[5],
            "pearson_range": [row[3], row[4]],
            "bio_dist_vs_pearson_dissim_rho": row[1],
            "input_rna_vs_layer_pearson_rho": row[2],
        }
        for row in layer_rows
    ]

    out = Path(args.out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines += [
        "## 生物学距离 vs Pearson 不一致度 (1−r)",
        "",
        "Spearman ρ：正值=Pearson 随细胞系生物学差异增大而**下降**（符合预期）。",
        "",
        "| 层次 | n | Pearson范围 | bio_dist↔(1−Pearson) | 输入RNA↔层Pearson |",
        "|---|---:|---|---:|---:|",
    ]
    for row in layer_rows:
        layer, bio_r, rna_r, pmin, pmax, n = row
        lines.append(
            f"| {layer} | {n} | {pmin:.3f}–{pmax:.3f} | "
            f"{bio_r if bio_r is not None else 'NA':.3f} | "
            f"{rna_r if rna_r is not None else 'NA':.3f} |"
        )

    lines += [
        "",
        "## 解读要点",
        "",
        "1. **promoter.active / RNA**：与输入 ENCODE RNA 直接相关，Pearson 跨系分化最大，且与输入 RNA 相关最高 → **最符合**细胞系类型差异。",
        "2. **protein**：p53 下游 hub 堆积，Pearson 仍偏高但 CD4 分离可见；生物距离一致性中等。",
        "3. **histone_mark / metabolite**：ChIP/代谢节点多由网络结构驱动，跨系 Pearson 偏高、与生物距离一致性弱。",
        "4. **全 place Spearman 全层虚高**：秩相关在各层均难反映谱系距离；分层后 **Pearson 才有判别力**。",
        "",
    ]

    for layer in ("promoter.active", "RNA", "protein", "histone_mark", "metabolite"):
        blk = report["by_layer"].get(layer)
        if not blk:
            continue
        lines.append(f"### {layer}（Pearson 低→高）")
        lines.append("")
        for p in blk["pairs_sorted_by_pearson"][:5]:
            lines.append(
                f"- {p['a']}–{p['b']} ({p['types']}): r={p['pearson']:.3f}, "
                f"bio_d={p['bio_distance']:.0f}, input_rna={p.get('input_rna_pearson', float('nan')):.3f}"
            )
        lines.append("")

    md_path = Path(args.out_md) if args.out_md else out.with_suffix(".md")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    for layer, block in report["by_layer"].items():
        labs = block["matrix"]["labels"]
        print(f"\n## {layer} (n={block['n_places']})")
        for metric in ("pearson", "spearman"):
            print(f"  [{metric}]")
            print("\t" + "\t".join(labs))
            for a in labs:
                row = [a] + [
                    f"{block['matrix'][metric][a][b]:.3f}" for b in labs
                ]
                print("\t".join(row))
    print("wrote", out)
    print("wrote", md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

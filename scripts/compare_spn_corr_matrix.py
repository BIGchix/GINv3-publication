#!/usr/bin/env python3
"""Spearman / Pearson correlation matrix across multiple SPN trajectory finals (+ optional level TSVs)."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def _spearman(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")

    def ranks(vals: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: vals[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    deny = math.sqrt(sum((b - my) ** 2 for b in ry))
    if denx == 0 or deny == 0:
        return float("nan")
    return num / (denx * deny)


def _load_traj_final(path: Path) -> dict[str, float]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {p: float(v) for p, v in zip(data["places"], data["mean_markings"][-1])}


def _load_traj_plateau_mean(path: Path, *, n_blocks: int) -> dict[str, float]:
    data = json.loads(path.read_text(encoding="utf-8"))
    places = data["places"]
    mean = data["mean_markings"]
    if not mean:
        return {}
    window = mean[-n_blocks:] if n_blocks > 0 else mean
    n = len(window)
    acc = [0.0] * len(places)
    for row in window:
        for i, v in enumerate(row):
            acc[i] += float(v)
    inv = 1.0 / n
    return {p: acc[i] * inv for i, p in enumerate(places)}


def _load_levels(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open(encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                out[row["place"]] = float(row["level"])
            except (KeyError, ValueError):
                continue
    return out


def _pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    denx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    deny = math.sqrt(sum((b - my) ** 2 for b in ys))
    if denx == 0 or deny == 0:
        return float("nan")
    return num / (denx * deny)


def _corr_matrix(
    series: dict[str, dict[str, float]], *, union: bool = False
) -> dict:
    labels = list(series)
    spearman = {a: {} for a in labels}
    pearson = {a: {} for a in labels}
    n_shared = {a: {} for a in labels}
    for a in labels:
        for b in labels:
            if union:
                keys = sorted(set(series[a]) | set(series[b]))
                xa = [series[a].get(k, 0.0) for k in keys]
                yb = [series[b].get(k, 0.0) for k in keys]
            else:
                keys = sorted(set(series[a]) & set(series[b]))
                xa = [series[a][k] for k in keys]
                yb = [series[b][k] for k in keys]
            spearman[a][b] = _spearman(xa, yb)
            pearson[a][b] = _pearson(xa, yb)
            n_shared[a][b] = len(keys)
    return {"labels": labels, "spearman": spearman, "pearson": pearson, "n_shared": n_shared}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--traj",
        action="append",
        nargs=2,
        metavar=("LABEL", "JSON"),
        help="SPN trajectory JSON (repeatable)",
    )
    ap.add_argument(
        "--levels",
        action="append",
        nargs=2,
        metavar=("LABEL", "TSV"),
        help="place levels TSV (repeatable)",
    )
    ap.add_argument("--out-json", required=True)
    ap.add_argument(
        "--traj-aggregate",
        choices=("final", "plateau_mean"),
        default="final",
        help="final=last block; plateau_mean=mean over last --plateau-blocks",
    )
    ap.add_argument(
        "--plateau-blocks",
        type=int,
        default=100,
        help="blocks averaged for plateau_mean (default: last 100 of 200)",
    )
    args = ap.parse_args()

    report: dict = {}
    if args.traj:
        loader = (
            (lambda p: _load_traj_plateau_mean(p, n_blocks=args.plateau_blocks))
            if args.traj_aggregate == "plateau_mean"
            else _load_traj_final
        )
        series = {lab: loader(Path(p)) for lab, p in args.traj}
        traj_key = (
            f"trajectory_plateau_mean_last{args.plateau_blocks}"
            if args.traj_aggregate == "plateau_mean"
            else "trajectory_final"
        )
        report[traj_key] = _corr_matrix(series, union=False)
        report["meta"] = {
            "traj_aggregate": args.traj_aggregate,
            "plateau_blocks": args.plateau_blocks,
        }
    if args.levels:
        series = {lab: _load_levels(Path(p)) for lab, p in args.levels}
        report["input_levels_intersect"] = _corr_matrix(series, union=False)
        report["input_levels_union0"] = _corr_matrix(series, union=True)

    out = Path(args.out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # compact print
    for block_name, block in report.items():
        if block_name == "meta" or "spearman" not in block:
            continue
        labs = block["labels"]
        for metric in ("spearman", "pearson"):
            print(f"## {block_name} ({metric})")
            print("\t" + "\t".join(labs))
            for a in labs:
                row = [a] + [f"{block[metric][a][b]:.3f}" for b in labs]
                print("\t".join(row))
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

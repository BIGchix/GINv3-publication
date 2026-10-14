"""Cross-cell-line layer correlation clustermaps (promoter.active / histone)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def _spearman_cols(X: np.ndarray) -> np.ndarray:
    """Spearman correlation between columns of X (n_samples × n_features)."""
    n, p = X.shape
    if p < 2 or n < 2:
        return np.eye(max(p, 1))
    ranks = np.zeros_like(X, dtype=float)
    for j in range(p):
        col = X[:, j]
        order = np.argsort(col, kind="mergesort")
        r = np.empty(n, dtype=float)
        i = 0
        while i < n:
            k = i
            while k + 1 < n and col[order[k + 1]] == col[order[i]]:
                k += 1
            avg = (i + k) / 2.0 + 1.0
            for t in range(i, k + 1):
                r[order[t]] = avg
            i = k + 1
        ranks[:, j] = r
    corr = np.corrcoef(ranks, rowvar=False)
    if corr.ndim == 0:
        return np.array([[float(corr)]])
    return corr


def plateau_mean_from_csv(
    path: Path,
    *,
    start: int,
    end: int,
) -> dict[str, float]:
    df = pd.read_csv(path)
    if "block" not in df.columns:
        raise ValueError(f"missing block column in {path}")
    sl = df["block"].between(start, end - 1)
    sub = df.loc[sl]
    if sub.empty:
        sub = df.tail(max(end - start, 1))
    out: dict[str, float] = {}
    for col in df.columns:
        if col == "block":
            continue
        out[col] = float(sub[col].mean())
    return out


def cell_place_matrix_from_trajectories(
    cell_csv: dict[str, Path],
    places: list[str],
    *,
    plateau_start: int = 50,
    plateau_end: int = 150,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Build matrix rows=cells, cols=places (plateau mean tokens)."""
    headers = [set(pd.read_csv(pth, nrows=0).columns) for pth in cell_csv.values()]
    common = set.intersection(*headers) - {"block"}
    used = [p for p in places if p in common]
    if not used:
        used = sorted(common)

    labels = list(cell_csv.keys())
    rows: list[list[float]] = []
    for lab in labels:
        vals = plateau_mean_from_csv(
            cell_csv[lab], start=plateau_start, end=plateau_end
        )
        rows.append([float(vals.get(p, 0.0)) for p in used])
    return np.asarray(rows, dtype=float), labels, used


def plot_layer_corr_from_cell_matrix(
    X: np.ndarray,
    labels: list[str],
    *,
    out_path: Path,
    layer: str,
    figsize: tuple[float, float] = (5.4, 5.0),
) -> dict:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import squareform

    corr = _spearman_cols(X)
    np.fill_diagonal(corr, 1.0)
    dist = 1.0 - np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    linkage = hierarchy.linkage(condensed, method="average")
    order = list(hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(linkage, condensed)))

    offdiag = corr[~np.eye(corr.shape[0], dtype=bool)]
    mean_offdiag = float(np.nanmean(offdiag)) if offdiag.size else float("nan")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cg = sns.clustermap(
        corr,
        row_linkage=linkage,
        col_linkage=linkage,
        cmap="vlag",
        vmin=-1,
        vmax=1,
        linewidths=0.2,
        figsize=figsize,
        cbar_kws={"label": "Spearman ρ"},
    )
    cg.ax_heatmap.set_title(f"{layer} · cross-place (cells={','.join(labels)})", fontsize=10)
    cg.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close("all")

    return {
        "path": str(out_path),
        "order": [int(x) for x in order],
        "mean_offdiag": mean_offdiag,
        "linkage": [[float(v) for v in row] for row in linkage.tolist()],
        "cell_labels": labels,
    }


def export_layer_corr_bundle(
    cell_csv: dict[str, Path],
    places: list[str],
    *,
    layer: str,
    out_dir: Path,
    prefix: str,
    plateau_start: int = 50,
    plateau_end: int = 150,
) -> dict:
    """Save correlation matrix + dendrogram order for Figure 4 replotting."""
    X, cell_labels, used = cell_place_matrix_from_trajectories(
        cell_csv,
        places,
        plateau_start=plateau_start,
        plateau_end=plateau_end,
    )
    corr = _spearman_cols(X)
    np.fill_diagonal(corr, 1.0)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{prefix}_corr_{layer.replace('.', '_')}"

    # cell × place plateau matrix
    mat_tsv = out_dir / f"{stem}_cell_place_plateau.tsv"
    with mat_tsv.open("w", encoding="utf-8") as fh:
        fh.write("cell\t" + "\t".join(used) + "\n")
        for i, lab in enumerate(cell_labels):
            fh.write(lab + "\t" + "\t".join(f"{X[i, j]:.6g}" for j in range(len(used))) + "\n")

    corr_tsv = out_dir / f"{stem}_matrix.tsv"
    with corr_tsv.open("w", encoding="utf-8") as fh:
        fh.write("place\t" + "\t".join(used) + "\n")
        for i, pi in enumerate(used):
            fh.write(pi + "\t" + "\t".join(f"{corr[i, j]:.6g}" for j in range(len(used))) + "\n")

    png = out_dir / f"{stem}_clustermap.png"
    plot_info = plot_layer_corr_from_cell_matrix(
        X, cell_labels, out_path=png, layer=layer
    )

    bundle = {
        "layer": layer,
        "plateau_blocks": [plateau_start, plateau_end],
        "cell_labels": cell_labels,
        "places": used,
        "corr_spearman": corr.tolist(),
        "cell_place_matrix": X.tolist(),
        "clustermap_order": plot_info["order"],
        "linkage": plot_info["linkage"],
        "mean_offdiag": plot_info["mean_offdiag"],
        "paths": {
            "cell_place_plateau_tsv": str(mat_tsv),
            "corr_matrix_tsv": str(corr_tsv),
            "clustermap_png": str(png),
        },
    }
    json_path = out_dir / f"{stem}_bundle.json"
    json_path.write_text(json.dumps(bundle, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez_compressed(
        out_dir / f"{stem}_bundle.npz",
        places=np.array(used),
        cell_labels=np.array(cell_labels),
        X=X,
        corr=corr,
        order=np.array(plot_info["order"], dtype=int),
    )
    bundle["paths"]["bundle_json"] = str(json_path)
    bundle["paths"]["bundle_npz"] = str(out_dir / f"{stem}_bundle.npz")
    return bundle

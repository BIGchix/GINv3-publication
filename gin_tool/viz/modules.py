"""Module community heatmap and MDS bubble chart."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def load_module_partition(modules_json: Path, names: list[str]) -> dict[int, int]:
    """Map subgraph node index -> module id."""
    data = json.loads(modules_json.read_text(encoding="utf-8"))
    name_to_idx = {n: i for i, n in enumerate(names)}
    mod_by_node: dict[int, int] = {}
    for mod in data:
        mid = int(mod["id"])
        for nm in mod["nodes"]:
            if nm in name_to_idx:
                mod_by_node[name_to_idx[nm]] = mid
    return mod_by_node


def render_module_heatmap(viz_dir: Path, distances_tsv: Path) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    df = pd.read_csv(distances_tsv, sep="\t")
    mods = sorted(set(df["module_i"]).union(df["module_j"]))
    m = len(mods)
    mat = np.full((m, m), np.nan)
    for _, row in df.iterrows():
        i, j = int(row["module_i"]), int(row["module_j"])
        d = row["min_geodesic"]
        if pd.isna(d) or row.get("unreachable", 0):
            val = np.nan
        else:
            val = float(d)
        mat[i, j] = mat[j, i] = val
    np.fill_diagonal(mat, 0.0)

    cg = sns.clustermap(
        mat,
        cmap="viridis_r",
        linewidths=0.5,
        figsize=(max(5, m * 0.6), max(4, m * 0.5)),
        cbar_kws={"label": "min geodesic"},
    )
    cg.savefig(viz_dir / "module_heatmap.png", dpi=120, bbox_inches="tight")
    plt.close("all")


def render_module_bubble(
    viz_dir: Path,
    modules_json: Path,
    distances_tsv: Path,
    seed_names: set[str],
) -> None:
    import json
    import matplotlib.pyplot as plt
    import pandas as pd
    from sklearn.manifold import MDS

    mods = json.loads(modules_json.read_text(encoding="utf-8"))
    n_mod = len(mods)
    if n_mod < 2 or n_mod > 20:
        return

    df = pd.read_csv(distances_tsv, sep="\t")
    mat = np.full((n_mod, n_mod), 0.0)
    for _, row in df.iterrows():
        i, j = int(row["module_i"]), int(row["module_j"])
        d = row["min_geodesic"]
        if pd.isna(d) or row.get("unreachable", 0):
            d = float(n_mod) + 1.0
        mat[i, j] = mat[j, i] = float(d)

    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    xy = mds.fit_transform(mat)

    sizes = []
    seed_fracs = []
    for mod in mods:
        nodes = mod["nodes"]
        sizes.append(len(nodes))
        n_seed = sum(1 for nm in nodes if nm in seed_names)
        seed_fracs.append(n_seed / max(len(nodes), 1))

    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(
        xy[:, 0],
        xy[:, 1],
        s=[40 + 12 * s for s in sizes],
        c=seed_fracs,
        cmap="YlOrRd",
        alpha=0.85,
        edgecolors="#333",
    )
    for i, mod in enumerate(mods):
        ax.annotate(f"M{mod['id']}", (xy[i, 0], xy[i, 1]), fontsize=9, ha="center")
    plt.colorbar(sc, ax=ax, label="seed fraction in module")
    ax.set_title("Module MDS (size ~ |module|)")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(viz_dir / "module_bubble.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

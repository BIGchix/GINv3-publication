"""Null debiased context bullet-style charts."""

from __future__ import annotations

import json
from pathlib import Path


def render_null_summary_bullet(viz_dir: Path, null_json: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    data = json.loads(null_json.read_text(encoding="utf-8"))
    block = data.get("subgraph_size_vs_null", {})
    if not block:
        return

    metrics = []
    obs, lo, hi = [], [], []
    for key in ("subgraph_nodes", "subgraph_edges_sym_nnz", "sym_mean_deg_proxy"):
        if key not in block:
            continue
        q = block[key]
        metrics.append(key.replace("_", "\n"))
        obs.append(q["observed"])
        lo.append(q["null_p05"])
        hi.append(q["null_p95"])

    if not metrics:
        return

    fig, ax = plt.subplots(figsize=(8, max(2, len(metrics) * 0.8)))
    y = np.arange(len(metrics))
    ax.hlines(y, lo, hi, color="#3498db", lw=6, alpha=0.5, label="null p05–p95")
    ax.scatter(obs, y, color="#c0392b", zorder=3, s=80, label="observed")
    ax.set_yticks(y)
    ax.set_yticklabels(metrics)
    ax.set_xlabel("value")
    ax.set_title("Observed vs random-seed null")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(viz_dir / "null_summary_bullet.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

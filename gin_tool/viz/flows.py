"""Cross-layer Sankey from flows.json."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from gin_tool.mapping import gin_entity_layer, is_virtual_gin_node


def _compress_layers(path_nodes: list[str]) -> list[str]:
    layers: list[str] = []
    for nm in path_nodes:
        if is_virtual_gin_node(nm):
            continue
        ly = gin_entity_layer(nm)
        if ly == "virtual":
            continue
        if not layers or layers[-1] != ly:
            layers.append(ly)
    return layers


def render_layer_sankey(viz_dir: Path, flows_json: Path, *, max_paths: int = 20) -> None:
    from plotly import graph_objects as go

    paths = json.loads(flows_json.read_text(encoding="utf-8"))[:max_paths]
    if not paths:
        return

    # Aggregate layer transitions (compressed consecutive layers)
    trans_count: dict[tuple[str, str], int] = defaultdict(int)
    layer_set: set[str] = set()
    for p in paths:
        seq = _compress_layers(p["nodes"])
        if len(seq) < 2:
            continue
        for a, b in zip(seq[:-1], seq[1:]):
            trans_count[(a, b)] += 1
            layer_set.add(a)
            layer_set.add(b)

    if not trans_count:
        return

    labels = sorted(layer_set)
    label_idx = {lb: i for i, lb in enumerate(labels)}
    sources, targets, values = [], [], []
    for (a, b), v in trans_count.items():
        sources.append(label_idx[a])
        targets.append(label_idx[b])
        values.append(v)

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(label=labels, pad=15, thickness=18),
                link=dict(source=sources, target=targets, value=values),
            )
        ]
    )
    fig.update_layout(title_text="Cross-layer flow transitions (aggregated)", height=500)
    fig.write_html(str(viz_dir / "layer_sankey.html"), include_plotlyjs="cdn")


def render_flow_overlap_heatmap(viz_dir: Path, flow_overlap_json: Path) -> None:
    import json

    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    data = json.loads(flow_overlap_json.read_text(encoding="utf-8"))
    mat = np.asarray(data["pairwise_jaccard"], dtype=np.float64)
    if mat.size == 0:
        return
    labels = [f"P{i}" for i in range(mat.shape[0])]
    fig, ax = plt.subplots(figsize=(max(5, mat.shape[0] * 0.35), max(4, mat.shape[0] * 0.32)))
    sns.heatmap(
        mat,
        xticklabels=labels,
        yticklabels=labels,
        cmap="YlOrRd",
        vmin=0.0,
        vmax=1.0,
        square=True,
        ax=ax,
        cbar_kws={"label": "Jaccard"},
    )
    ax.set_title("Path overlap (Jaccard)")
    fig.tight_layout()
    fig.savefig(viz_dir / "flow_overlap_heatmap.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

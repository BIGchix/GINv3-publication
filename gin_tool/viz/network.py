"""Subgraph network: NetworkX layout + Plotly interactive HTML (+ optional PNG)."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse

from gin_tool.mapping import gin_entity_layer
from gin_tool.viz._select import select_viz_node_indices_named
from gin_tool.viz._style import LAYER_COLORS, MODULE_COLORS, SIGN_EDGE_COLORS
from gin_tool.viz.export import export_cytoscape_sif


def _node_sizes(scores: np.ndarray, indices: list[int]) -> np.ndarray:
    s = scores[indices]
    smax = float(s.max()) if s.size and float(s.max()) > 0 else 1.0
    return 12.0 + 28.0 * np.log1p(s / (smax * 1e-6 + 1e-12))


def render_subgraph_network(
    viz_dir: Path,
    sub_sym: sparse.csr_matrix,
    sub_dir: sparse.csr_matrix,
    names: list[str],
    scores: np.ndarray,
    seed_indices: list[int],
    *,
    max_nodes: int = 200,
    edge_mode: str = "none",
    write_static: bool = True,
    module_by_node: dict[int, int] | None = None,
) -> None:
    from plotly import graph_objects as go

    seed_set = set(seed_indices)
    node_set = select_viz_node_indices_named(
        sub_sym, scores, seed_indices, names, max_nodes=max_nodes
    )
    idx_list = sorted(node_set)
    n_sub = sub_sym.shape[0]
    if len(node_set) < n_sub:
        folded = n_sub - len(node_set)
        print(
            f"viz: showing {len(node_set)}/{n_sub} nodes "
            f"(--viz-max-nodes {max_nodes}; folded {folded})",
            flush=True,
        )

    G = nx.Graph()
    for i in idx_list:
        G.add_node(i)
    sym = sub_sym.tocsr()
    for i in idx_list:
        lo, hi = sym.indptr[i], sym.indptr[i + 1]
        for j in sym.indices[lo:hi]:
            j = int(j)
            if j in node_set and j != i:
                G.add_edge(i, j)

    if G.number_of_nodes() == 0:
        return

    iterations = 50 if G.number_of_nodes() <= 300 else 30
    pos = nx.spring_layout(G, seed=42, iterations=iterations)

    edge_traces = []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        color = "#cccccc"
        if edge_mode == "sign":
            from gin_tool.viz.export import _dominant_sign

            sign = _dominant_sign(sub_dir, u, v)
            color = SIGN_EDGE_COLORS.get(sign, "#cccccc")
        edge_traces.append(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode="lines",
                line=dict(width=0.6, color=color),
                hoverinfo="none",
                showlegend=False,
            )
        )

    layers = [gin_entity_layer(names[i]) for i in idx_list]
    layer_unique = sorted(set(layers))
    node_traces = []
    for layer in layer_unique:
        in_layer = [i for i in idx_list if gin_entity_layer(names[i]) == layer]
        xs = [pos[i][0] for i in in_layer]
        ys = [pos[i][1] for i in in_layer]
        sizes = _node_sizes(scores, in_layer)
        texts = [names[i] for i in in_layer]
        custom = [
            f"score={scores[i]:g}<br>layer={layer}<br>seed={i in seed_set}"
            for i in in_layer
        ]
        marker_kw: dict = dict(
            size=sizes,
            color=LAYER_COLORS.get(layer, "#7f8c8d"),
            line=dict(width=1, color="#333"),
        )
        if module_by_node:
            mcols = [
                MODULE_COLORS[module_by_node.get(i, 0) % len(MODULE_COLORS)]
                for i in in_layer
            ]
            marker_kw["line"] = dict(width=2, color=mcols)
        for i in in_layer:
            if i in seed_set:
                marker_kw.setdefault("symbol", [])
        node_traces.append(
            go.Scatter(
                x=xs,
                y=ys,
                mode="markers+text",
                text=[t if i in seed_set or scores[i] > np.percentile(scores, 95) else "" for i, t in zip(in_layer, texts)],
                textposition="top center",
                textfont=dict(size=9),
                marker=marker_kw,
                name=layer,
                hovertext=custom,
                hoverinfo="text",
            )
        )

    fig = go.Figure(data=edge_traces + node_traces)
    fig.update_layout(
        title="GIN k-hop subgraph (node size ~ RWR score)",
        showlegend=True,
        hovermode="closest",
        margin=dict(l=20, r=20, t=50, b=20),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        width=1000,
        height=700,
    )
    html_path = viz_dir / "subgraph_network.html"
    fig.write_html(str(html_path), include_plotlyjs="cdn")

    if write_static:
        try:
            fig.write_image(str(viz_dir / "subgraph_network.png"), scale=2)
        except Exception as exc:
            print(f"viz: could not write PNG (install kaleido): {exc}", flush=True)

    export_cytoscape_sif(
        viz_dir / "subgraph_for_cytoscape.sif",
        sub_sym,
        sub_dir,
        names,
        scores,
        seed_set,
        node_set,
        module_by_node,
        edge_mode=edge_mode,
    )

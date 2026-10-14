"""Visualization pipeline for GIN CLI outputs (optional ``[viz]`` extras)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import sparse

from gin_tool.mapping import MappingQCReport


def _require_viz_extras() -> None:
    try:
        import matplotlib  # noqa: F401
        import plotly  # noqa: F401
        import seaborn  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Visualization requires optional dependencies. Install with: "
            "pip install -e '.[viz]'"
        ) from exc


def run_visualizations(
    viz_dir: Path,
    *,
    sub_sym: sparse.csr_matrix,
    sub_dir: sparse.csr_matrix,
    scores: np.ndarray,
    names: list[str],
    seed_indices: list[int],
    report: MappingQCReport,
    seed_gin_names: list[str],
    main_out_dir: Path,
    analyze: bool = False,
    max_nodes: int = 200,
    edge_mode: str = "none",
    write_static: bool = True,
    null_precompute_dir: Path | None = None,
    k: int | None = None,
    n_seeds: int | None = None,
    rwr_alpha: float | None = None,
) -> None:
    """
    Generate plots under ``viz_dir`` from an in-memory subgraph run.

    Reads ``modules.json`` / ``flows.json`` from ``main_out_dir`` when ``analyze`` is True.
    """
    _require_viz_extras()
    import matplotlib

    matplotlib.use("Agg")
    viz_dir.mkdir(parents=True, exist_ok=True)

    from gin_tool.viz.network import render_subgraph_network
    from gin_tool.viz.qc import render_qc_dashboard
    from gin_tool.viz.scores import render_ecdf_null_comparison, render_score_distribution

    seed_set = set(seed_indices)
    seed_name_set = set(seed_gin_names)

    module_by_node: dict[int, int] | None = None
    modules_json = main_out_dir / "modules.json"
    if analyze and modules_json.is_file():
        from gin_tool.viz.modules import load_module_partition

        module_by_node = load_module_partition(modules_json, names)

    render_qc_dashboard(viz_dir, report)
    render_score_distribution(viz_dir, names, scores, seed_set)

    null_json = main_out_dir / "null_debiased_context.json"
    render_ecdf_null_comparison(
        viz_dir,
        scores,
        null_precompute_dir=null_precompute_dir,
        k=k,
        n_seeds=n_seeds or len(seed_gin_names),
        rwr_alpha=rwr_alpha,
    )
    if null_json.is_file():
        from gin_tool.viz.null_charts import render_null_summary_bullet

        render_null_summary_bullet(viz_dir, null_json)

    render_subgraph_network(
        viz_dir,
        sub_sym,
        sub_dir,
        names,
        scores,
        seed_indices,
        max_nodes=max_nodes,
        edge_mode=edge_mode,
        write_static=write_static,
        module_by_node=module_by_node,
    )

    if analyze:
        dist_tsv = main_out_dir / "module_pairwise_distances.tsv"
        if modules_json.is_file() and dist_tsv.is_file():
            from gin_tool.viz.modules import render_module_bubble, render_module_heatmap

            render_module_heatmap(viz_dir, dist_tsv)
            render_module_bubble(viz_dir, modules_json, dist_tsv, seed_name_set)

        flows_json = main_out_dir / "flows.json"
        if flows_json.is_file():
            from gin_tool.viz.flows import render_layer_sankey

            render_layer_sankey(viz_dir, flows_json)

        overlap_json = main_out_dir / "flow_overlap.json"
        if overlap_json.is_file():
            from gin_tool.viz.flows import render_flow_overlap_heatmap

            render_flow_overlap_heatmap(viz_dir, overlap_json)

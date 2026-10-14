"""Orchestrate module detection, flow ranking, module–flow coupling, and SVG schematic."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from gin_tool.flows import (
    _path_layer_score,
    enumerate_simple_paths_from_starts,
    rank_cross_layer_paths,
    write_flow_overlap_json,
    write_flows_json,
)
from gin_tool.module_flow import (
    interface_real_nodes,
    path_module_incidence_matrix,
    write_module_flow_json,
)
from gin_tool.modules import (
    louvain_communities_subgraph,
    module_pairwise_min_geodesic,
    modules_to_json,
    write_module_distance_tsv,
)
from gin_tool.schematic import render_schematic_svg
from gin_tool.subgraph import (
    extract_khop_region,
    extract_khop_region_mixed,
    remap_seed_indices_to_subgraph,
)


def _safe_louvain(adj_sym: sparse.csr_matrix) -> list[frozenset[int]]:
    if adj_sym.shape[0] == 0:
        return []
    try:
        return louvain_communities_subgraph(adj_sym)
    except Exception:
        n_comp, lab = connected_components(adj_sym, directed=False)
        return [frozenset(np.flatnonzero(lab == c).astype(int)) for c in range(n_comp)]


def _observed_louvain_modularity(
    adj_sym: sparse.csr_matrix, modules: list[frozenset[int]]
) -> tuple[int, float | None, float | None]:
    """Number of communities, largest fraction of vertices, Newman modularity."""
    n = int(adj_sym.shape[0])
    if not modules or n == 0:
        return 0, None, None
    sizes = [len(c) for c in modules]
    largest_frac = max(sizes) / float(n)
    try:
        G = nx.from_scipy_sparse_array(adj_sym, create_using=nx.Graph)
        comms = [set(c) for c in modules]
        q = float(nx.community.modularity(G, comms, weight="weight", resolution=1))
    except Exception:
        q = None
    return len(modules), largest_frac, q


def run_module_flow_analysis(
    out_dir: Path,
    sifi_path: str,
    resolved_seed_names: list[str],
    k: int | None = None,
    *,
    layer_map: dict[str, str] | None = None,
    user_k: int | None = None,
    max_lines: int | None = None,
    max_path_depth: int = 22,
    top_flows: int = 35,
    max_total_paths: int = 8000,
) -> dict[str, Any]:
    """
    Write analysis artifacts under ``out_dir``:

    - ``modules.json`` — Louvain (or fallback connected) communities on the k-hop **sym** graph.
    - ``module_pairwise_distances.tsv`` — min geodesic between module pairs.
    - ``flows.json`` — top cross-layer directed paths on the induced **dir** subgraph.
    - ``flow_overlap.json`` — Jaccard matrix + parallel / intersect labels.
    - ``module_flow_relation.json`` — path×module incidence + interface real nodes.
    - ``schematic.svg`` — rounded module boxes, flow lines, interface labels.
    Returns observed module/path statistics for optional null debiasing (see
    ``gin_tool.debias.run_debias_pipeline``).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if layer_map is not None:
        sub_sym, sub_dir, old_to_new, new_to_old, ed = extract_khop_region_mixed(
            sifi_path,
            resolved_seed_names,
            layer_map,
            max_lines=max_lines,
            user_k=user_k,
        )
    else:
        if k is None:
            raise ValueError("k is required when layer_map is not provided")
        sub_sym, sub_dir, old_to_new, new_to_old, ed = extract_khop_region(
            sifi_path, resolved_seed_names, k, max_lines=max_lines
        )
    n = sub_sym.shape[0]
    sub_names = [ed.index_to_name[new_to_old[i]] for i in range(n)]
    name_fn = lambda i: sub_names[int(i)]

    gidx = [ed.name_to_index[name] for name in resolved_seed_names if name in ed.name_to_index]
    sub_seeds = remap_seed_indices_to_subgraph(old_to_new, gidx)
    if not sub_seeds:
        raise ValueError("No seed indices in subgraph after mapping.")

    modules = _safe_louvain(sub_sym)
    modules_to_json(str(out_dir / "modules.json"), modules, sub_names)
    n_mod, largest_frac, mod_q = _observed_louvain_modularity(sub_sym, modules)

    if len(modules) >= 2:
        dist, unr = module_pairwise_min_geodesic(sub_sym, modules)
        write_module_distance_tsv(
            str(out_dir / "module_pairwise_distances.tsv"), modules, dist, unr
        )

    raw = enumerate_simple_paths_from_starts(
        sub_dir,
        sub_seeds,
        max_depth=max_path_depth,
        max_paths_per_start=500,
        max_total_paths=max_total_paths,
    )
    ranked = rank_cross_layer_paths(raw, name_fn, top_k=top_flows)
    write_flows_json(str(out_dir / "flows.json"), ranked, name_fn)
    write_flow_overlap_json(str(out_dir / "flow_overlap.json"), ranked)

    inc = path_module_incidence_matrix(ranked, modules)
    ifaces = interface_real_nodes(ranked, modules, sub_sym, name_fn)
    write_module_flow_json(
        str(out_dir / "module_flow_relation.json"),
        inc,
        ifaces,
        sub_names,
    )

    svg_w = max(960.0, 140.0 * max(len(modules), 1) + 220.0)
    render_schematic_svg(
        out_dir / "schematic.svg",
        modules,
        ranked,
        sub_names,
        ifaces,
        width=svg_w,
    )

    best_score: tuple[int, int] | None = None
    if ranked:
        best_score = _path_layer_score(ranked[0], name_fn)
    return {
        "obs_louvain_n": n_mod,
        "obs_louvain_largest_frac": largest_frac,
        "obs_louvain_modularity": mod_q,
        "obs_raw_paths": len(raw),
        "obs_ranked_paths": len(ranked),
        "obs_best_path_score": best_score,
    }

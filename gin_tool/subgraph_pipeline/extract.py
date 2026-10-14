"""Orchestrate scBionet-style subgraph extraction with bipartite (virtual-preserving) output."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from gin_tool.active_tf import ActiveTFRequest, infer_active_tfs
from gin_tool.gin_defaults import DEFAULT_NODE_LIST_DIR, DEFAULT_SIFI
from gin_tool.io import build_indexed_edges, collect_unique_node_names, directed_csr_from_indexed
from gin_tool.subgraph_pipeline.connectivity import repair_connectivity
from gin_tool.subgraph_pipeline.edge_topology import apply_edge_topology
from gin_tool.subgraph_pipeline.gin_self_loop import SeedSelfLoopQC, augment_walk_with_seed_self_loops
from gin_tool.subgraph_pipeline.manifest import (
    CentralDogmaMappingQC,
    EdgeTopologyQC,
    SeedManifest,
    SubgraphExtractResult,
)
from gin_tool.subgraph_pipeline.regulatory_train import filter_indexed_edges_regulatory_train
from gin_tool.node_lists import load_node_list_dir
from gin_tool.subgraph_pipeline.seed_gin_mapping import (
    collect_regulatory_seed_nodes,
    collect_walk_input_nodes,
    map_union_seeds,
)
from gin_tool.subgraph_pipeline.signaling_bipartite import signaling_bipartite_subgraph


def _log_pipeline_warning(message: str) -> None:
    print(f"[WARN] {message}", flush=True)


def _warn_if_node_budget_hit(
    *,
    hit_cap: bool,
    max_real_nodes: int,
    real_node_count: int,
    seed_gin_nodes: list[str],
    walk_real_nodes: set[str],
    warnings: list[str],
    stage: str,
) -> None:
    if not hit_cap and real_node_count < max_real_nodes:
        return
    msg = (
        f"[{stage}] hit max_real_nodes={max_real_nodes} "
        f"(real_nodes={real_node_count}); walk truncated"
    )
    missed = sorted(set(seed_gin_nodes) - walk_real_nodes)
    if missed:
        msg += f"; {len(missed)} seed GIN nodes missing (sample: {missed[:10]})"
    warnings.append(msg)
    _log_pipeline_warning(msg)


@dataclass
class SubgraphPipelineConfig:
    max_real_hops: int = 6
    max_real_nodes: int = 10_000
    active_tf_top_n: int = 80
    active_tf_fdr: float = 0.05
    active_tf_min_overlap: int = 2
    active_tf_low_power_threshold: int = 5
    # If set, skip ranking and use these TF symbols / GIN ids.
    active_tf_user_list: list[str] | None = None
    # Pipeline default B (rwr_cohesion); C = rwr_rewire; fisher = legacy.
    active_tf_method: str = "rwr_cohesion"
    active_tf_rwr_n_null: int = 1000
    active_tf_rwr_n_jobs: int = 8
    active_tf_rwr_emp_p: float = 0.05
    active_tf_rwr_rewire_sweeps: float = 0.1
    active_tf_work_dir: str = ""
    prefilter_regulatory: bool = True
    filter_regulatory: bool = False
    regulatory_relaxed: bool = True
    enable_connectivity_l2: bool = True
    bridge_hops: int = 4
    require_seed_self_loops: bool = True
    seed_self_loop_hops: int = 12
    canonical_edges_only: bool = False
    soft_noncan: bool = False
    soft_noncan_max_fraction: float = 0.25
    gin_sifi: str = ""
    node_list_dir: str = ""
    max_lines: int | None = None
    # Online PubChem only for explicit metabolite / unmapped common-name seeds.
    # Default off so KEGG gene batches stay offline and within memory caps.
    allow_online_metabolite: bool = False


def extract_subgraph(
    seed_manifest: SeedManifest,
    config: SubgraphPipelineConfig | None = None,
    *,
    gin_nodes_set: set[str] | None = None,
    indexed_edges=None,
) -> SubgraphExtractResult:
    """
    Extract a bipartite subgraph for ``seed_manifest``.

    Optional ``gin_nodes_set`` / ``indexed_edges`` avoid reloading GINv*.sifi on
    every call (important for batch jobs under a memory cap).
    """
    cfg = config or SubgraphPipelineConfig()
    warnings: list[str] = []

    sifi = Path(cfg.gin_sifi or DEFAULT_SIFI)
    node_list_dir = Path(cfg.node_list_dir or DEFAULT_NODE_LIST_DIR)
    if not sifi.is_file():
        raise FileNotFoundError(f"GIN sifi missing: {sifi}")

    if gin_nodes_set is None:
        gin_nodes_set = collect_unique_node_names(str(sifi))
    gin_is_member = gin_nodes_set.__contains__

    metabolite_pool = set(load_node_list_dir(node_list_dir).get("metabolite", []))
    # GINv2 UC_IDs in the live sifi often exceed nodeList_metabolite (rebuild only
    # auto-classified CAS*). Union real chem nodes from the graph when available.
    try:
        from gin_tool.subgraph_pipeline.metabolite_id_resolve import collect_chem_nodes_from_sifi

        metabolite_pool |= collect_chem_nodes_from_sifi(sifi)
    except OSError:
        pass
    seed_records, cd_qc = map_union_seeds(
        seed_manifest.union_genes,
        gin_is_member,
        seed_kinds=seed_manifest.seed_kinds or None,
        metabolite_pool=metabolite_pool,
        allow_online_metabolite=cfg.allow_online_metabolite,
    )
    seed_manifest.gin_anchor_nodes = sorted(collect_walk_input_nodes(seed_records))
    seed_manifest.seed_kinds = {
        rec.symbol: rec.kind for rec in seed_records if rec.walk_anchor_nodes or rec.kind == "metabolite"
    }
    if getattr(cd_qc, "n_metabolite", 0):
        seed_manifest.stats["n_metabolite_mapped"] = cd_qc.n_metabolite
    if cd_qc.n_unmapped:
        warnings.append(
            f"{cd_qc.n_unmapped} seed genes unmapped (sample: {cd_qc.unmapped_symbols[:10]})"
        )
    if cd_qc.n_triplet_complete < cd_qc.n_input:
        warnings.append(
            f"triplet incomplete: {cd_qc.n_input - cd_qc.n_triplet_complete}/{cd_qc.n_input} seeds"
        )

    walk_input = collect_walk_input_nodes(seed_records)
    regulatory_seeds = collect_regulatory_seed_nodes(seed_records)
    seed_gin_nodes = sorted(walk_input)
    if not seed_gin_nodes:
        raise RuntimeError("No seed genes mapped to GIN nodes")

    work_dir = cfg.active_tf_work_dir or ""
    active_result = infer_active_tfs(
        ActiveTFRequest(
            seed_gin_nodes=seed_gin_nodes,
            sifi_path=str(sifi),
            node_list_dir=node_list_dir,
            user_active_tfs=cfg.active_tf_user_list or None,
            gin_node_is_member=gin_is_member if cfg.active_tf_user_list else None,
            top_n=cfg.active_tf_top_n,
            fdr_alpha=cfg.active_tf_fdr,
            min_overlap=cfg.active_tf_min_overlap,
            low_power_threshold=cfg.active_tf_low_power_threshold,
            max_lines=cfg.max_lines,
            method=cfg.active_tf_method,
            rwr_work_dir=Path(work_dir) if work_dir else Path("gin_out") / "active_tf_rwr_work",
            rwr_n_null=cfg.active_tf_rwr_n_null,
            rwr_n_jobs=cfg.active_tf_rwr_n_jobs,
            rwr_emp_p_cutoff=cfg.active_tf_rwr_emp_p,
            rwr_rewire_sweeps=cfg.active_tf_rwr_rewire_sweeps,
        )
    )
    active_syms = active_result.active_symbols()
    active_gin = active_result.active_gin_nodes()
    if active_result.warnings:
        warnings.extend(active_result.warnings[:10])

    # Full-GIN QC: warn immediately if any seed promoter is not reached by active TFs.
    from gin_tool.active_tf_promoter_qc import (
        emit_active_tf_promoter_qc,
        qc_active_tf_seed_promoter_connectivity,
    )

    seed_genes_for_qc = sorted(
        {rec.symbol for rec in seed_records if rec.kind == "gene" and rec.promoter_node}
    )
    if seed_genes_for_qc and active_syms:
        promoter_qc = qc_active_tf_seed_promoter_connectivity(
            sifi,
            active_syms,
            seed_genes_for_qc,
            max_lines=cfg.max_lines,
        )
        qc_dir = Path(work_dir) if work_dir else Path("gin_out") / "active_tf_rwr_work"
        warnings.extend(
            emit_active_tf_promoter_qc(
                promoter_qc,
                out_json=qc_dir / "active_tf_seed_promoter_connectivity.json",
                print_warning=True,
            )
        )

    from gin_tool.active_tf import load_tf_symbol_pool

    tf_pool = load_tf_symbol_pool(node_list_dir)
    ed = indexed_edges if indexed_edges is not None else build_indexed_edges(
        str(sifi), max_lines=cfg.max_lines
    )
    seed_symbols = set(seed_manifest.union_genes)

    if cfg.prefilter_regulatory:
        work_ed, prefilter_stats = filter_indexed_edges_regulatory_train(
            ed,
            active_tf_names=set(active_syms) | set(active_gin),
            seed_symbols=seed_symbols,
            seed_gin_names=regulatory_seeds,
            tf_pool=tf_pool,
            relaxed=cfg.regulatory_relaxed,
        )
        use_runtime_filter = cfg.filter_regulatory
    else:
        work_ed = ed
        prefilter_stats = None
        use_runtime_filter = cfg.filter_regulatory

    adj_dir = directed_csr_from_indexed(work_ed)
    walk = signaling_bipartite_subgraph(
        set(walk_input),
        adj_dir,
        work_ed.index_to_name,
        work_ed.name_to_index,
        max_real_hops=cfg.max_real_hops,
        max_real_nodes=cfg.max_real_nodes,
        filter_regulatory=use_runtime_filter,
        active_tf_names=set(active_syms) | set(active_gin),
        seed_symbols=seed_symbols,
        seed_gin_names=regulatory_seeds,
    )

    _warn_if_node_budget_hit(
        hit_cap=getattr(walk, "hit_max_real_nodes", False),
        max_real_nodes=cfg.max_real_nodes,
        real_node_count=len(walk.real_nodes),
        seed_gin_nodes=seed_gin_nodes,
        walk_real_nodes=walk.real_nodes,
        warnings=warnings,
        stage="signaling_walk",
    )

    self_loop_qc = SeedSelfLoopQC()
    n_seeds_no_loop = 0
    if cfg.require_seed_self_loops:
        walk, self_loop_qc = augment_walk_with_seed_self_loops(
            walk,
            seed_records,
            adj_dir,
            work_ed.index_to_name,
            work_ed.name_to_index,
            max_real_hops=cfg.max_real_hops,
            max_real_nodes=cfg.max_real_nodes,
            filter_regulatory=use_runtime_filter,
            active_tf_names=set(active_syms) | set(active_gin),
            seed_symbols=seed_symbols,
            seed_gin_names=regulatory_seeds,
            self_loop_hops=cfg.seed_self_loop_hops,
        )
        n_seeds_no_loop = self_loop_qc.n_failed
        if self_loop_qc.n_failed:
            msg = (
                f"{self_loop_qc.n_failed}/{self_loop_qc.n_seeds} seeds lack self-loop "
                f"(sample: {self_loop_qc.failed_symbols_sample[:10]})"
            )
            warnings.append(msg)
            _log_pipeline_warning(msg)
        if getattr(walk, "hit_max_real_nodes", False):
            _warn_if_node_budget_hit(
                hit_cap=True,
                max_real_nodes=cfg.max_real_nodes,
                real_node_count=len(walk.real_nodes),
                seed_gin_nodes=seed_gin_nodes,
                walk_real_nodes=walk.real_nodes,
                warnings=warnings,
                stage="seed_self_loop",
            )

    conn, conn_qc = repair_connectivity(
        walk,
        seed_gin_nodes,
        adj_dir,
        work_ed.index_to_name,
        work_ed.name_to_index,
        enable_l2=cfg.enable_connectivity_l2,
        bridge_hops=cfg.bridge_hops,
        n_seeds_no_loop=n_seeds_no_loop,
    )
    walk = conn.walk

    protect: set[str] = set(seed_gin_nodes) | set(active_gin)
    for rec in seed_records:
        for node in (rec.rna_node, rec.protein_node, rec.promoter_node):
            if node:
                protect.add(node)

    walk, edge_qc, topo_warnings = apply_edge_topology(
        walk,
        protect_nodes=protect,
        tf_pool=tf_pool,
        canonical_edges_only=cfg.canonical_edges_only,
        soft_noncan=cfg.soft_noncan,
        soft_noncan_max_fraction=cfg.soft_noncan_max_fraction,
        n_seeds_no_loop=n_seeds_no_loop,
        adj_dir=adj_dir,
        index_to_name=work_ed.index_to_name,
        name_to_index=work_ed.name_to_index,
    )
    warnings.extend(topo_warnings)

    return SubgraphExtractResult(
        run_id=seed_manifest.run_id,
        seed_manifest=seed_manifest,
        walk=walk,
        active_tf_symbols=active_syms,
        gin_seed_nodes=seed_gin_nodes,
        central_dogma_qc=cd_qc,
        seed_self_loop_qc=self_loop_qc,
        connectivity_qc=conn_qc,
        edge_topology_qc=edge_qc,
        prefilter_stats=prefilter_stats,
        warnings=warnings,
    )

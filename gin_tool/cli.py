"""Command-line entry for US1: seeds -> subgraph TSV + RWR scores + mapping QC JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gin_tool.gin_defaults import (
    DEFAULT_NODE_LIST_DIR,
    DEFAULT_SIFI,
    DEFAULT_SYMBOL_ALIASES,
    DEFAULT_SYMBOL_BLACKLIST,
    REPO_ROOT,
)
from gin_tool.io import collect_unique_node_names
from gin_tool.mapping import (
    load_symbol_aliases,
    load_symbol_blacklist,
)
from gin_tool.analyze import run_module_flow_analysis
from gin_tool.propagate import random_walk_restart, scores_to_named_tsv
from gin_tool.subgraph import (
    extract_khop_region_mixed,
    extract_khop_subgraph_mixed,
    remap_seed_indices_to_subgraph,
    subgraph_edges_as_tuples,
    subgraph_from_named_edges,
    write_subgraph_tsv,
)


def _default_node_list_dir() -> Path | None:
    return DEFAULT_NODE_LIST_DIR if DEFAULT_NODE_LIST_DIR.is_dir() else None


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=(
            "GIN k-hop subgraph + RWR (US1). Per-seed k from nodeList "
            "(protein/tf/hme→1, others→2) unless --k overrides all. "
            "Mixed protein+RNA seeds supported (union subgraph). "
            "Runtime null debias only (no precompute)."
        ),
    )
    try:
        _sifi_help_default = str(DEFAULT_SIFI.relative_to(REPO_ROOT))
    except ValueError:
        _sifi_help_default = str(DEFAULT_SIFI)
    p.add_argument(
        "--sifi",
        type=Path,
        default=DEFAULT_SIFI,
        help=f"Path to GINv*.sifi or .sifi.gz (default: {_sifi_help_default})",
    )
    p.add_argument(
        "--seeds",
        required=True,
        help="Comma-separated gene symbols or GIN node names",
    )
    p.add_argument(
        "--k",
        type=int,
        default=None,
        help="Override k-hop for all seeds (default: per-seed from nodeList layer)",
    )
    p.add_argument(
        "--max-lines",
        type=int,
        default=None,
        help="Only read first N edges (for testing; omit for full file)",
    )
    p.add_argument(
        "--symbol-blacklist",
        type=Path,
        default=DEFAULT_SYMBOL_BLACKLIST,
        help="HGNC symbols to exclude from mapping (one per line; missing file is OK)",
    )
    p.add_argument(
        "--symbol-aliases",
        type=Path,
        default=DEFAULT_SYMBOL_ALIASES,
        help="Optional old_symbol\\tapproved_symbol TSV applied before mapping",
    )
    p.add_argument("--out-dir", type=Path, default=Path("gin_out"), help="Output directory")
    p.add_argument("--alpha", type=float, default=0.3, help="RWR restart probability")
    p.add_argument(
        "--analyze",
        action="store_true",
        help=(
            "Also run module (Louvain) + cross-layer directed flows + "
            "module–flow relation + schematic.svg in out-dir"
        ),
    )
    p.add_argument(
        "--debias",
        action="store_true",
        default=True,
        help="Runtime null debiasing (default: on)",
    )
    p.add_argument(
        "--no-debias",
        action="store_false",
        dest="debias",
        help="Skip null debiasing outputs",
    )
    p.add_argument(
        "--null-precompute-dir",
        type=Path,
        default=None,
        help=argparse.SUPPRESS,
    )
    p.add_argument(
        "--null-tail-alpha",
        type=float,
        default=0.05,
        help="Upper-tail p threshold for null_significant on debiased RWR scores",
    )
    p.add_argument(
        "--null-debias-mode",
        choices=("residual_median", "residual_p90", "divide_median"),
        default="residual_median",
        help="How to compute score_debiased from raw RWR and null pool",
    )
    p.add_argument(
        "--null-runtime-replicates",
        type=int,
        default=25,
        help="Random-seed RWR replicates for runtime null calibration",
    )
    p.add_argument(
        "--node-list-dir",
        type=Path,
        default=None,
        help="nodeList/ directory (default: <repo>/nodeList). Required for auto-k and debias.",
    )
    p.add_argument(
        "--max-path-depth",
        type=int,
        default=22,
        help="Max vertices per enumerated directed path (analyze mode)",
    )
    p.add_argument(
        "--top-flows",
        type=int,
        default=35,
        help="Number of top cross-layer paths to keep (analyze mode)",
    )
    p.add_argument(
        "--viz",
        action="store_true",
        help="Write visualization plots under out-dir/viz/ (requires pip install -e '.[viz]')",
    )
    p.add_argument(
        "--viz-max-nodes",
        type=int,
        default=200,
        help="Max nodes in subgraph network plot (seeds + top RWR; default 200)",
    )
    p.add_argument(
        "--viz-out-dir",
        type=Path,
        default=None,
        help="Visualization output directory (default: <out-dir>/viz/)",
    )
    p.add_argument(
        "--viz-edge-mode",
        choices=("none", "sign"),
        default="none",
        help="Edge coloring in network plot: uniform or signed",
    )
    p.add_argument(
        "--viz-no-static",
        action="store_true",
        help="Skip PNG exports (HTML only)",
    )
    p.add_argument(
        "--active-tfs",
        default=None,
        help="User-supplied active TF symbols or GIN nodes (comma-separated). Skips ranking.",
    )
    p.add_argument(
        "--active-tfs-file",
        type=Path,
        default=None,
        help="File with one active TF symbol/GIN node per line (skips ranking)",
    )
    p.add_argument(
        "--active-tf-method",
        choices=("rwr_cohesion", "rwr_rewire", "fisher"),
        default="rwr_cohesion",
        help=(
            "Active-TF ranking when no --active-tfs/--active-tfs-file: "
            "rwr_cohesion=B (default), rwr_rewire=C (slower/stricter), fisher=legacy"
        ),
    )
    p.add_argument(
        "--active-tf-n-jobs",
        type=int,
        default=8,
        help="CPU workers for RWR null replicates (B/C; default 8)",
    )
    p.add_argument(
        "--active-tf-n-null",
        type=int,
        default=1000,
        help="RWR null replicates for B/C (default 1000)",
    )
    p.add_argument(
        "--skip-active-tf",
        action="store_true",
        help="Do not infer or write active TF outputs",
    )
    p.add_argument("--tf-top", type=int, default=30, help="Max active TFs when enriching (default 30)")
    p.add_argument("--tf-fdr", type=float, default=0.05, help="FDR threshold for regulon enrich")
    p.add_argument(
        "--tf-min-overlap",
        type=int,
        default=2,
        help="Minimum seed overlap with TF regulon (default 2)",
    )
    p.add_argument(
        "--subgraph-mode",
        choices=("context", "signaling", "khop"),
        default="context",
        help=(
            "Subgraph extraction: context (prefilter TF→DNA + seed-only BFS), "
            "signaling (GINtoSPN real-hop + prefilter), or legacy k-hop"
        ),
    )
    p.add_argument(
        "--context-include-active-tf-anchors",
        action="store_true",
        help=(
            "Legacy context: add active TFs as BFS anchors (default: seeds only, "
            "regulatory pre-filter on graph)"
        ),
    )
    p.add_argument(
        "--context-mech-depth",
        type=int,
        default=2,
        help=(
            "Max directed BFS depth on the bipartite GIN (graph hops, default 2). "
            "On strict R↔V GIN: depth=1 reaches only adjacent virtuals; "
            "depth=2 ≡ one real-hop (R→V→R). Prefer quoting both, e.g. "
            "'graph-hop=2 (real-hop=1)'."
        ),
    )
    p.add_argument(
        "--context-max-nodes",
        type=int,
        default=100_000,
        help="Cap nodes during context subgraph expansion (default 100000)",
    )
    p.add_argument(
        "--no-promoter-reverse-walk",
        action="store_true",
        help=(
            "Disable reverse directed BFS from seed promoters (default: on). "
            "When on, seed promoters expand along incoming edges to reach TFs; "
            "other seeds still expand forward."
        ),
    )
    p.add_argument(
        "--no-lcc-bridge",
        action="store_true",
        help=(
            "Disable post-context LCC island bridge (default: on for context mode). "
            "Bridge uses graph-hop on directed GIN (default 2; not real-hop)."
        ),
    )
    p.add_argument(
        "--lcc-bridge-hops",
        type=int,
        default=2,
        help="Graph-hop budget for LCC BFS island bridge (default 2)",
    )
    p.add_argument(
        "--context-pipeline-version",
        type=str,
        default="v4",
        help=(
            "Context build recipe: v4 (default; HME soft-seeds only at 0b, "
            "met hyperedge enzymes post-LCC, then HME/histprune → dogma) "
            "or v3 (fallback: dogma then HME/histprune). "
            "Also accepts full version strings from gin_tool.pipeline_version."
        ),
    )
    p.add_argument(
        "--context-n-jobs",
        type=int,
        default=1,
        help=(
            "CPU workers for context build parallel steps (default 1 = serial). "
            "Parallelizes: per-seed BFS, role-seed met enzymes, met-enzyme "
            "recovery, HME directed∥undirected hop BFS, per-gene dogma "
            "(process pool). 0 or negative → auto (min(CPU,32) / abs)."
        ),
    )
    p.add_argument(
        "--no-dogma-chain-repair",
        action="store_true",
        help=(
            "Disable post-prune repair of missing *.RNA/*.promoter chains for "
            "subgraph proteins (default: on for context mode). Also tries reverse "
            "real-hop1 from promoter.active to TFs already in the subgraph."
        ),
    )
    p.add_argument(
        "--no-hme-chain-repair",
        action="store_true",
        help=(
            "Disable post-prune histone/HME chain repair (default: on for context). "
            "Backfills hist marks for promoter.active, bridges subgraph PPI → outside "
            "HME writers, and restores HME→histone edges."
        ),
    )
    p.add_argument(
        "--hme-max-real-hops",
        type=int,
        default=3,
        help="Max real-hop budget for PPI→HME bridge (default 3)",
    )
    p.add_argument(
        "--hme-upstream-graph-hops",
        type=int,
        default=4,
        help="Reverse graph-hop budget from histone marks to discover HME (default 4)",
    )
    p.add_argument(
        "--no-eraser-chain-repair",
        action="store_true",
        help=(
            "Disable eraser catalytic-chain repair inside HME repair (default: on). "
            "When on, pulls substrate→product demethylation/deacetylation intermediates "
            "and enzyme edges even if products lack promoter;state links in GIN."
        ),
    )
    p.add_argument(
        "--eraser-max-rounds",
        type=int,
        default=3,
        help="Max demethylation/deacetylation ladder rounds for eraser repair (default 3)",
    )
    p.add_argument(
        "--no-hme-active-hist-backfill",
        action="store_true",
        help=(
            "Disable HME-repair pre-step that pulls same-gene histone marks from GIN "
            "for pathway-seed promoter.active (default: seed-restricted backfill on)."
        ),
    )
    p.add_argument(
        "--hme-backfill-all-actives",
        action="store_true",
        help=(
            "Backfill hist marks for every promoter.active (legacy). Default is "
            "seed-gene-only backfill and drop non-seed promoter.active/;state "
            "(RNA/protein kept)."
        ),
    )
    p.add_argument(
        "--no-drop-nonseed-promoter-active",
        action="store_true",
        help=(
            "Keep non-seed promoter.active when doing seed-restricted hist backfill "
            "(default: drop non-seed active/;state)."
        ),
    )
    p.add_argument(
        "--no-hme-sourceless-hist-prune",
        action="store_true",
        help=(
            "Disable post-HME cleanup that drops histone marks with no full-GIN "
            "in-edge (e.g. source-only H3K9ac) and resulting orphan islands "
            "(default: cleanup on)."
        ),
    )
    p.add_argument(
        "--signaling-max-real-hops",
        type=int,
        default=12,
        help="Max real-hop depth for signaling subgraph (default 12)",
    )
    p.add_argument(
        "--signaling-max-real-nodes",
        type=int,
        default=100_000,
        help="Cap real nodes during signaling walk (default 100000)",
    )
    p.add_argument(
        "--signaling-no-prefilter",
        action="store_true",
        help="Disable protein/TF→DNA pre-filter for signaling mode",
    )
    p.add_argument(
        "--export-spn",
        action="store_true",
        help="Export signaling subgraph to mEPN GraphML for BioLayout SPN (signaling mode)",
    )
    p.add_argument(
        "--spn-exclude-chr",
        action="store_true",
        help="Exclude chromatin locus (chr*) places from SPN GraphML export",
    )
    args = p.parse_args(argv)

    if args.null_precompute_dir is not None:
        print(
            "warning: --null-precompute-dir is ignored; debias uses runtime null only.",
            file=__import__("sys").stderr,
        )

    raw_seeds = [s.strip() for s in args.seeds.split(",") if s.strip()]
    blocked = load_symbol_blacklist(args.symbol_blacklist)
    aliases = load_symbol_aliases(args.symbol_aliases)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    node_list_dir = args.node_list_dir or _default_node_list_dir()
    if node_list_dir is None:
        p.error(
            "nodeList/ is required for layer-aware k and runtime null. "
            "Pass --node-list-dir or add nodeList/ at the repo root."
        )

    gin_set = collect_unique_node_names(str(args.sifi), max_lines=args.max_lines)

    # Multi-layer gene anchors: promoter.active + RNA + protein (all present layers).
    # Metabolites stay single-node. No single-layer prefer / first-hit collapse.
    from gin_tool.mapping import MappingQCReport, gin_entity_layer
    from gin_tool.subgraph_pipeline.seed_gin_mapping import (
        collect_walk_input_nodes,
        map_union_seeds,
    )

    prepared: list[str] = []
    skipped_blacklisted: list[str] = []
    for s in raw_seeds:
        if s in blocked:
            skipped_blacklisted.append(s)
            continue
        prepared.append(aliases.get(s, s))

    seed_records, dogma_qc = map_union_seeds(
        prepared,
        gin_set.__contains__,
        allow_online_metabolite=False,
    )
    resolved = sorted(collect_walk_input_nodes(seed_records))
    layer_counts_anchors: dict[str, int] = {}
    for gnode in resolved:
        ly = gin_entity_layer(gnode)
        layer_counts_anchors[ly] = layer_counts_anchors.get(ly, 0) + 1

    report = MappingQCReport(
        gin_version=args.sifi.name,
        input_count=len(raw_seeds),
        mapped_count=len(resolved),
        unmapped_ids=list(dogma_qc.unmapped_symbols),
        unmapped_detail=[
            {"id": u, "reason": "no_central_dogma_or_metabolite_node"}
            for u in dogma_qc.unmapped_symbols
        ],
        ambiguous_mappings=[],
        layer_counts=layer_counts_anchors,
        warnings=[
            "seed mapping=central_dogma_multi_layer "
            "(gene→promoter.active+RNA+protein when present; no prefer)",
            *([f"blacklisted skipped: {skipped_blacklisted!r}"] if skipped_blacklisted else []),
        ],
    )
    qc_payload = {
        **json.loads(report.to_json()),
        "strategy": "central_dogma_multi_layer",
        "central_dogma": dogma_qc.to_dict(),
        "n_walk_anchors": len(resolved),
        "seed_records": [r.to_dict() for r in seed_records],
    }
    (args.out_dir / "qc_mapping.json").write_text(
        json.dumps(qc_payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    seed_manifest = {
        "run_id": args.out_dir.name,
        "source": "cli_central_dogma",
        "union_genes": [r.symbol for r in seed_records],
        "seed_kinds": {r.symbol: r.kind for r in seed_records},
        "gin_anchor_nodes": resolved,
        "stats": {
            **dogma_qc.to_dict(),
            "n_input_raw": len(raw_seeds),
            "n_walk_anchors": len(resolved),
            "n_gene": dogma_qc.n_gene,
            "n_metabolite": dogma_qc.n_metabolite,
        },
    }
    (args.out_dir / "seed_manifest.json").write_text(
        json.dumps(seed_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    if dogma_qc.unmapped_symbols:
        print(
            f"warning: skipping {len(dogma_qc.unmapped_symbols)} unmapped seeds "
            f"(see qc_mapping.json): {list(dogma_qc.unmapped_symbols)[:20]!r}",
            file=__import__("sys").stderr,
        )
    if not resolved:
        p.error("No valid seeds after central-dogma multi-layer mapping.")
    print(
        f"[seeds] central_dogma anchors={len(resolved)} "
        f"genes={dogma_qc.n_gene} mets={dogma_qc.n_metabolite} "
        f"rna={dogma_qc.n_with_rna} protein={dogma_qc.n_with_protein} "
        f"promoter={dogma_qc.n_with_promoter} triplet={dogma_qc.n_triplet_complete}",
        flush=True,
    )
    from gin_tool.active_tf import (
        ActiveTFRequest,
        infer_active_tfs,
        load_tf_symbol_pool,
        parse_id_list,
        seed_nodes_to_gene_symbols,
        write_active_tf_qc_json,
        write_active_tfs_tsv,
    )
    from gin_tool.io import build_indexed_edges
    from gin_tool.node_lists import build_pools_and_layer_map, summarize_seed_layers

    active_tf_result = None
    if not args.skip_active_tf:
        user_tfs = parse_id_list(args.active_tfs, path=args.active_tfs_file)
        active_tf_result = infer_active_tfs(
            ActiveTFRequest(
                seed_gin_nodes=resolved,
                sifi_path=str(args.sifi),
                node_list_dir=node_list_dir,
                user_active_tfs=user_tfs or None,
                gin_node_is_member=gin_set.__contains__ if user_tfs else None,
                top_n=args.tf_top,
                fdr_alpha=args.tf_fdr,
                min_overlap=args.tf_min_overlap,
                max_lines=args.max_lines,
                method=args.active_tf_method,
                rwr_work_dir=args.out_dir / "active_tf_rwr",
                rwr_n_null=args.active_tf_n_null,
                rwr_n_jobs=args.active_tf_n_jobs,
            )
        )
        write_active_tfs_tsv(args.out_dir / "active_tfs.tsv", active_tf_result)
        write_active_tf_qc_json(args.out_dir / "active_tf_qc.json", active_tf_result)

        # Full-GIN QC: active TFs must reach seed promoters (warn + suggest expanding list).
        from gin_tool.active_tf_promoter_qc import (
            emit_active_tf_promoter_qc,
            qc_active_tf_seed_promoter_connectivity,
            seed_genes_from_gin_nodes,
        )

        seed_genes_for_qc = seed_genes_from_gin_nodes(resolved)
        if not seed_genes_for_qc:
            seed_genes_for_qc = sorted(seed_nodes_to_gene_symbols(resolved))
        if seed_genes_for_qc and active_tf_result.active_symbols():
            promoter_qc = qc_active_tf_seed_promoter_connectivity(
                args.sifi,
                active_tf_result.active_symbols(),
                seed_genes_for_qc,
                max_lines=args.max_lines,
            )
            emit_active_tf_promoter_qc(
                promoter_qc,
                out_json=args.out_dir / "active_tf_seed_promoter_connectivity.json",
            )
            run_meta_promoter_qc = promoter_qc.to_dict()
        else:
            run_meta_promoter_qc = None
    else:
        run_meta_promoter_qc = None

    ed_nl = build_indexed_edges(str(args.sifi), max_lines=args.max_lines)
    _, layer_map, _, _, _ = build_pools_and_layer_map(
        node_list_dir, ed_nl.index_to_name, ed_nl.name_to_index
    )
    layer_counts, layer_by_name, k_by_name = summarize_seed_layers(
        resolved, layer_map, user_k=args.k
    )
    unique_ks = set(k_by_name.values())
    mixed_layers = len({ly for ly in layer_counts if layer_counts[ly] > 0}) > 1
    mixed_k = len(unique_ks) > 1
    k_notes: list[str] = []
    if args.k is not None:
        k_notes.append(f"User override: all seeds use k={args.k}.")
    elif mixed_k:
        k_notes.append(
            "Per-seed k-hop union subgraph: "
            + ", ".join(f"{n}→k={k_by_name[n]}" for n in resolved)
        )
    else:
        k_notes.append(f"Uniform k={next(iter(unique_ks))} for all seeds.")

    k_repr = max(k_by_name.values())

    from gin_tool.pipeline_version import (
        CONTEXT_PIPELINE_ALIASES,
        CONTEXT_PIPELINE_VERSION,
        CONTEXT_PIPELINE_VERSION_V3,
        CONTEXT_PIPELINE_VERSION_V4,
        SPN_PIPELINE_VERSION,
    )

    ctx_ver_raw = str(getattr(args, "context_pipeline_version", "v4") or "v4").strip()
    if ctx_ver_raw not in CONTEXT_PIPELINE_ALIASES:
        p.error(
            f"Unknown --context-pipeline-version {ctx_ver_raw!r}; "
            f"expected one of {sorted(CONTEXT_PIPELINE_ALIASES)}"
        )
    context_pipeline_version = CONTEXT_PIPELINE_ALIASES[ctx_ver_raw]
    is_context_v4 = context_pipeline_version == CONTEXT_PIPELINE_VERSION_V4
    from gin_tool.subgraph_pipeline.parallel_util import resolve_n_jobs

    context_n_jobs = resolve_n_jobs(getattr(args, "context_n_jobs", 1))

    run_meta = {
        "context_pipeline_version": context_pipeline_version,
        "context_pipeline_version_default": CONTEXT_PIPELINE_VERSION,
        "context_pipeline_version_v3_fallback": CONTEXT_PIPELINE_VERSION_V3,
        "context_n_jobs": context_n_jobs,
        "spn_pipeline_version_default": SPN_PIPELINE_VERSION,
        "k_policy": "per_seed" if mixed_k else "uniform",
        "k": k_repr,
        "k_by_seed": k_by_name,
        "seed_layer_counts": layer_counts,
        "layer_by_seed": {n: layer_by_name[n] for n in resolved},
        "mixed_layers": mixed_layers,
        "mixed_k": mixed_k,
        "k_selection_notes": k_notes,
        "node_list_dir": str(node_list_dir.resolve()),
        "null_policy": "runtime_only",
        "debias": args.debias,
    }
    if active_tf_result is not None:
        run_meta["active_tf"] = {
            "method": active_tf_result.method,
            "skipped": active_tf_result.skipped,
            "skip_reason": active_tf_result.skip_reason,
            "active_tf_count": len(active_tf_result.records),
            "active_tf_symbols": active_tf_result.active_symbols(),
            "warnings": active_tf_result.warnings,
        }
    if run_meta_promoter_qc is not None:
        run_meta["active_tf_seed_promoter_qc"] = {
            "n_active_tfs": run_meta_promoter_qc["n_active_tfs"],
            "n_seed_genes": run_meta_promoter_qc["n_seed_genes"],
            "gh1_covered": (
                f"{run_meta_promoter_qc['graphhop1']['n_genes_covered']}/"
                f"{run_meta_promoter_qc['graphhop1']['n_seed_genes']}"
            ),
            "gh2_covered": (
                f"{run_meta_promoter_qc['graphhop2']['n_genes_covered']}/"
                f"{run_meta_promoter_qc['graphhop2']['n_seed_genes']}"
            ),
            "has_uncovered": run_meta_promoter_qc["has_uncovered"],
            "uncovered_genes": run_meta_promoter_qc["graphhop2"]["uncovered_genes"],
            "suggestion": run_meta_promoter_qc.get("suggestion"),
        }
    run_meta["subgraph_mode"] = args.subgraph_mode

    active_tf_gin = (
        active_tf_result.active_gin_nodes()
        if active_tf_result is not None
        else []
    )
    active_tf_symbols = (
        set(active_tf_result.active_symbols())
        if active_tf_result is not None
        else set()
    )
    seed_symbols = (
        active_tf_result.seed_symbols
        if active_tf_result is not None
        else seed_nodes_to_gene_symbols(resolved)
    )
    tf_pool = load_tf_symbol_pool(node_list_dir)

    primary_resolved = list(resolved)
    if args.subgraph_mode == "context" and is_context_v4:
        from gin_tool.subgraph_pipeline.role_seed_enrichment import enrich_role_soft_seeds

        print(f"[context] n_jobs={context_n_jobs} version={context_pipeline_version}", flush=True)
        resolved, role_qc = enrich_role_soft_seeds(
            ed_nl,
            primary_resolved,
            node_list_dir=node_list_dir,
            include_hme=True,
            include_met_enzymes=False,  # met hyperedge enzymes only at post-LCC 3a
            n_jobs=context_n_jobs,
        )
        run_meta["role_seed_enrichment"] = role_qc.to_dict()
        (args.out_dir / "role_seed_enrichment_qc.json").write_text(
            json.dumps(role_qc.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"[role-seeds] soft={len(role_qc.soft_seeds)} "
            f"hme={role_qc.n_hme_soft} met_enz={role_qc.n_met_enzyme_soft} "
            f"hub_skip={role_qc.n_hub_mets_skipped} "
            f"anchors={len(resolved)}",
            flush=True,
        )

    if args.subgraph_mode == "context":
        from gin_tool.context_subgraph import extract_context_region

        if not is_context_v4:
            print(
                f"[context] n_jobs={context_n_jobs} version={context_pipeline_version}",
                flush=True,
            )
        sub_sym, sub_dir, old_to_new, new_to_old, ctx_stats = extract_context_region(
            ed_nl,
            resolved,
            active_tf_gin,
            seed_symbols,
            mech_max_depth=args.context_mech_depth,
            max_nodes=args.context_max_nodes,
            active_tf_symbols=active_tf_symbols,
            tf_pool=tf_pool,
            include_active_tf_anchors=args.context_include_active_tf_anchors,
            promoter_reverse_walk=not args.no_promoter_reverse_walk,
            n_jobs=context_n_jobs,
        )
        ed = ed_nl
        # Export / SPN / RWR use directed edges. Do NOT export sub_sym: it is an
        # undirected support (A+A.T) and previously invented reverse arcs
        # (e.g. transcription→promoter) in subgraph_edges.tsv.
        sub = sub_dir
        run_meta["context_subgraph"] = ctx_stats.to_dict()
        run_meta["context_subgraph"]["edge_export"] = "directed"
    elif args.subgraph_mode == "signaling":
        from gin_tool.signaling_subgraph import (
            build_sym_subgraph_from_walk,
            signaling_real_subgraph_from_edges,
        )

        walk_result, _adj_filt, pf_stats = signaling_real_subgraph_from_edges(
            ed_nl,
            set(resolved),
            max_real_hops=args.signaling_max_real_hops,
            max_real_nodes=args.signaling_max_real_nodes,
            prefilter_regulatory=not args.signaling_no_prefilter,
            filter_regulatory=False,
            active_tf_names=active_tf_symbols | set(active_tf_gin),
            seed_symbols=seed_symbols,
            seed_gin_names=set(resolved),
            tf_pool=tf_pool,
        )
        sub_sym, name_to_local, ordered_names = build_sym_subgraph_from_walk(walk_result)
        old_to_new = {
            ed_nl.name_to_index[n]: name_to_local[n]
            for n in ordered_names
            if n in ed_nl.name_to_index
        }
        new_to_old = {j: ed_nl.name_to_index[n] for n, j in name_to_local.items()}
        sub_dir = sub_sym
        sub = sub_sym
        ed = ed_nl
        run_meta["signaling_subgraph"] = {
            "path_pairs": len(walk_result.path_pairs),
            "real_nodes": len(walk_result.real_nodes),
            "compressed_edges": len(walk_result.compressed_edges),
            "virtual_nodes_visited": walk_result.virtual_nodes_visited,
            "status": walk_result.status,
            "seed_only_input": True,
            "prefilter": pf_stats.to_dict() if pf_stats is not None else None,
        }
        signaling_edges_path = args.out_dir / "signaling_compressed_edges.tsv"
        with signaling_edges_path.open("w", encoding="utf-8") as sf:
            for a, b in sorted(walk_result.compressed_edges):
                sf.write(f"{a}\t{b}\n")

        if args.export_spn:
            from gin_tool.regulatory_filter import filter_indexed_edges_regulatory
            from gin_tool.spn_export import export_signaling_to_graphml

            export_ed = ed_nl
            if not args.signaling_no_prefilter:
                export_ed, _ = filter_indexed_edges_regulatory(
                    ed_nl,
                    active_tf_names=active_tf_symbols | set(active_tf_gin),
                    seed_symbols=seed_symbols,
                    seed_gin_names=set(resolved),
                    tf_pool=tf_pool,
                )
            real_nodes = set(walk_result.real_nodes)
            compressed = set(walk_result.compressed_edges)
            if args.spn_exclude_chr:
                real_nodes = {n for n in real_nodes if not n.startswith("chr")}
                compressed = {
                    (a, b)
                    for a, b in compressed
                    if a in real_nodes and b in real_nodes
                }
            spn_result = export_signaling_to_graphml(
                compressed,
                real_nodes,
                export_ed,
                args.out_dir / "spn.graphml",
                layer_map=layer_map,
                seed_gin_names=set(resolved),
                qc_json_path=args.out_dir / "spn_export_qc.json",
            )
            run_meta["spn_export"] = spn_result.qc_dict()
    elif args.viz:
        sub_sym, sub_dir, old_to_new, new_to_old, ed = extract_khop_region_mixed(
            str(args.sifi),
            resolved,
            layer_map,
            max_lines=args.max_lines,
            user_k=args.k,
        )
        sub = sub_sym
    else:
        sub, old_to_new, new_to_old, ed = extract_khop_subgraph_mixed(
            str(args.sifi),
            resolved,
            layer_map,
            max_lines=args.max_lines,
            user_k=args.k,
        )
        sub_sym, sub_dir = sub, sub

    edges = subgraph_edges_as_tuples(sub, new_to_old, ed.index_to_name)
    if args.subgraph_mode == "context" and not args.no_lcc_bridge:
        from gin_tool.subgraph_pipeline.lcc_multisource_bridge import apply_lcc_bridge_to_edges

        bridge_res = apply_lcc_bridge_to_edges(edges, ed, max_hops=args.lcc_bridge_hops)
        edges = bridge_res.edges
        sub, old_to_new, new_to_old = subgraph_from_named_edges(edges, ed)
        bridge_qc = bridge_res.to_dict()
        run_meta["lcc_bridge"] = bridge_qc
        (args.out_dir / "lcc_bridge_qc.json").write_text(
            json.dumps(bridge_qc, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    def _run_dogma_repair() -> None:
        nonlocal edges, sub, old_to_new, new_to_old
        from gin_tool.dogma_chain_repair import repair_dogma_chains_on_edges

        edges, dogma_qc = repair_dogma_chains_on_edges(
            edges, ed, tf_pool=tf_pool, n_jobs=context_n_jobs
        )
        sub, old_to_new, new_to_old = subgraph_from_named_edges(edges, ed)
        run_meta["dogma_chain_repair"] = dogma_qc.to_dict()
        (args.out_dir / "dogma_chain_repair_qc.json").write_text(
            json.dumps(dogma_qc.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"[dogma-repair] incomplete={dogma_qc.n_incomplete_proteins} "
            f"missing_tx={dogma_qc.n_missing_transcription} "
            f"prom_only={dogma_qc.n_promoter_only} "
            f"repaired={dogma_qc.n_repaired_genes} "
            f"(tx={dogma_qc.n_repaired_transcription} prom_only={dogma_qc.n_repaired_promoter_only}) "
            f"tf_linked={dogma_qc.n_with_tf_link} "
            f"+nodes={dogma_qc.n_nodes_added} +edges={dogma_qc.n_edges_added}",
            flush=True,
        )

    def _run_hme_repair() -> None:
        nonlocal edges, sub, old_to_new, new_to_old
        from gin_tool.hme_chain_repair import repair_hme_chains_with_node_lists

        # Primary seeds only: soft-role anchors must still reach primary seeds
        # under histprune / reverse ancestor prune.
        anchors = set(primary_resolved) if primary_resolved else set(resolved)
        from gin_tool.hme_chain_repair import gene_symbols_from_seed_inputs

        seed_kinds = {r.symbol: r.kind for r in seed_records}
        backfill_seeds = None
        if not args.hme_backfill_all_actives:
            backfill_seeds = gene_symbols_from_seed_inputs(
                [r.symbol for r in seed_records],
                seed_kinds=seed_kinds,
            )
        edges, hme_qc = repair_hme_chains_with_node_lists(
            edges,
            ed,
            node_list_dir,
            max_real_hops=args.hme_max_real_hops,
            upstream_graph_hops=args.hme_upstream_graph_hops,
            repair_eraser_chains=not args.no_eraser_chain_repair,
            eraser_max_rounds=args.eraser_max_rounds,
            prune_sourceless_histones=not args.no_hme_sourceless_hist_prune,
            backfill_active_hist_marks=not args.no_hme_active_hist_backfill,
            backfill_seed_genes=backfill_seeds,
            do_drop_nonseed_promoter_active=not args.no_drop_nonseed_promoter_active,
            anchor_nodes=anchors,
            n_jobs=context_n_jobs,
        )
        sub, old_to_new, new_to_old = subgraph_from_named_edges(edges, ed)
        run_meta["hme_chain_repair"] = hme_qc.to_dict()
        (args.out_dir / "hme_chain_repair_qc.json").write_text(
            json.dumps(hme_qc.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"[hme-repair] hist={hme_qc.n_hist_marks} "
            f"(before_bf={hme_qc.n_hist_marks_before_backfill} "
            f"+bf={hme_qc.n_hist_marks_backfilled} "
            f"active_genes={hme_qc.n_active_genes_with_backfill}/"
            f"{hme_qc.n_active_genes_scanned}) "
            f"seed_only={hme_qc.backfill_seed_genes_only} "
            f"drop_nonseed_active={hme_qc.n_nonseed_active_dropped} "
            f"outside_hme={hme_qc.n_hme_outside} "
            f"bridged={hme_qc.n_hme_bridged} hme→hist={hme_qc.n_hme_to_hist_paths} "
            f"eraser_steps={hme_qc.n_eraser_steps} eraser_enz={hme_qc.n_eraser_enzymes} "
            f"drop_src_hist={hme_qc.n_sourceless_hist_removed} "
            f"drop_dead_hist={hme_qc.n_dead_end_hist_removed} "
            f"drop_island={hme_qc.n_island_nodes_removed} "
            f"drop_stubs={hme_qc.n_catalysis_stubs_removed} "
            f"reprune={hme_qc.n_reprune_nodes_removed} "
            f"+nodes={hme_qc.n_nodes_added} +edges={hme_qc.n_edges_added}",
            flush=True,
        )

    def _run_met_enzyme_recovery() -> None:
        nonlocal edges, sub, old_to_new, new_to_old
        from gin_tool.subgraph_pipeline.met_enzyme_recovery import recover_met_enzymes_on_edges

        edges, mez_qc = recover_met_enzymes_on_edges(
            edges, ed, n_jobs=context_n_jobs
        )
        sub, old_to_new, new_to_old = subgraph_from_named_edges(edges, ed)
        run_meta["met_enzyme_recovery"] = mez_qc.to_dict()
        (args.out_dir / "met_enzyme_recovery_qc.json").write_text(
            json.dumps(mez_qc.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"[met-enzyme] hyper={mez_qc.n_hyperedges_scanned} "
            f"ref={mez_qc.n_enzymes_referenced} "
            f"+enz={mez_qc.n_enzymes_added} +edges={mez_qc.n_edges_added}",
            flush=True,
        )

    if args.subgraph_mode == "context":
        if is_context_v4:
            # v4: soft-seeded roles → metenz → dogma (promoter.active population)
            # → HME/active-hist backfill/eraser/histprune. Dogma before HME so
            # hist backfill sees the full promoter.active set, not only seeds.
            _run_met_enzyme_recovery()
            if not args.no_dogma_chain_repair:
                _run_dogma_repair()
            if not args.no_hme_chain_repair:
                _run_hme_repair()
        else:
            # v3 fallback: dogma then HME/histprune
            if not args.no_dogma_chain_repair:
                _run_dogma_repair()
            if not args.no_hme_chain_repair:
                _run_hme_repair()
    write_subgraph_tsv(str(args.out_dir / "subgraph_edges.tsv"), edges)

    (args.out_dir / "run_meta.json").write_text(
        json.dumps(run_meta, indent=2), encoding="utf-8"
    )

    rwr_seeds = primary_resolved if args.subgraph_mode == "context" else resolved
    gidx = [ed.name_to_index[name] for name in rwr_seeds if name in ed.name_to_index]
    sub_seeds = remap_seed_indices_to_subgraph(old_to_new, gidx)
    scores = random_walk_restart(sub, sub_seeds, alpha=args.alpha)
    scores_to_named_tsv(
        str(args.out_dir / "propagation_scores.tsv"),
        scores,
        ed.index_to_name,
        new_to_old=new_to_old,
    )
    analyze_stats = None
    if args.analyze:
        analyze_stats = run_module_flow_analysis(
            args.out_dir,
            str(args.sifi),
            rwr_seeds,
            layer_map=layer_map,
            user_k=args.k,
            max_lines=args.max_lines,
            max_path_depth=args.max_path_depth,
            top_flows=args.top_flows,
        )
    if args.debias:
        from gin_tool.debias import run_debias_pipeline

        sub_names = [ed.index_to_name[new_to_old[i]] for i in range(sub.shape[0])]
        kw: dict = {}
        if analyze_stats:
            kw = {
                "obs_louvain_n": analyze_stats.get("obs_louvain_n"),
                "obs_louvain_largest_frac": analyze_stats.get("obs_louvain_largest_frac"),
                "obs_louvain_modularity": analyze_stats.get("obs_louvain_modularity"),
                "obs_raw_paths": analyze_stats.get("obs_raw_paths"),
                "obs_ranked_paths": analyze_stats.get("obs_ranked_paths"),
                "obs_best_path_score": analyze_stats.get("obs_best_path_score"),
            }
        from gin_tool.node_lists import UNKNOWN_SEED_LAYER

        dominant_layer = None
        if layer_counts:
            dominant_layer = max(layer_counts.items(), key=lambda x: (x[1], x[0]))[0]
            if dominant_layer == UNKNOWN_SEED_LAYER:
                dominant_layer = None
        run_debias_pipeline(
            args.out_dir,
            sifi_path=str(args.sifi),
            k=k_repr,
            n_seeds=len(resolved),
            resolved_seed_names=resolved,
            scores=scores,
            sub_names=sub_names,
            sub_sym_n=int(sub.shape[0]),
            sub_sym_nnz=int(sub.nnz),
            null_precompute_dir=None,
            rwr_alpha=args.alpha,
            max_lines=args.max_lines,
            debias_mode=args.null_debias_mode,
            tail_alpha=args.null_tail_alpha,
            runtime_replicates=args.null_runtime_replicates,
            node_list_dir=node_list_dir,
            seed_layer=dominant_layer,
            seed_layer_counts=layer_counts,
            k_by_name=k_by_name,
            runtime_only=True,
            analyze=args.analyze,
            **kw,
        )
    if args.viz:
        viz_dir = args.viz_out_dir or (args.out_dir / "viz")
        sub_names = [ed.index_to_name[new_to_old[i]] for i in range(sub.shape[0])]
        try:
            from gin_tool.viz import run_visualizations

            run_visualizations(
                viz_dir,
                sub_sym=sub_sym,
                sub_dir=sub_dir,
                scores=scores,
                names=sub_names,
                seed_indices=sub_seeds,
                report=report,
                seed_gin_names=resolved,
                main_out_dir=args.out_dir,
                analyze=args.analyze,
                max_nodes=args.viz_max_nodes,
                edge_mode=args.viz_edge_mode,
                write_static=not args.viz_no_static,
                null_precompute_dir=None,
                k=k_repr,
                n_seeds=len(resolved),
                rwr_alpha=args.alpha,
            )
        except ImportError as exc:
            p.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

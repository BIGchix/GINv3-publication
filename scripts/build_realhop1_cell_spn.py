#!/usr/bin/env python3
"""Fixed pipeline: context → per-cell GRN allowlist prune → SigPN export.

SigPN = signaling Petri net (not stochastic PN). Historical filenames may use
``spn_*`` (e.g. ``spn_export_qc.json``).

Stages (under ``--out-root/{cell}/``):
  1. GRN allowlist prune from shared ``subgraph_edges.tsv``
  2. SigPN GraphML + ``spn_export_qc.json`` (``spn_std_h1_b0`` export knobs)

Does not run long MC unless ``--runs`` > 0.
See docs/CHIP_ALLOWLIST.md and docs/SIGPN_CUSTOM_SIMULATIONS.md.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from gin_tool.active_tf import load_tf_symbol_pool
from gin_tool.gin_defaults import DEFAULT_NODE_LIST_DIR
from gin_tool.subgraph_pipeline.grn_allowlist_prune import (
    chip_align_cell_dir,
    default_chip_align_root,
    load_measured_tfs,
    load_seed_genes,
    load_subgraph_edges_tsv,
    prune_grn_allowlist,
    write_pruned_cell_dir,
)
from gin_tool.subgraph_pipeline.spn_standard_policy import standard_h1_b0_simulate_kwargs
from gin_tool.subgraph_pipeline.spn_simulate import export_and_simulate_subgraph_dir


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--context-edges",
        type=Path,
        default=root / "gin_out/p53_fig3e_v4_realhop1_gh2/subgraph_edges.tsv",
        help="Shared context subgraph_edges.tsv (real-hop 1 ≡ graph-hop 2)",
    )
    ap.add_argument(
        "--context-dir",
        type=Path,
        default=root / "gin_out/p53_fig3e_v4_realhop1_gh2",
        help="Optional: copy seed_manifest.json into cell dir",
    )
    ap.add_argument(
        "--seed-file",
        type=Path,
        default=root / "seedLists/kegg_p53_signaling.txt",
    )
    ap.add_argument("--cell", default="h1_wa01")
    ap.add_argument(
        "--out-root",
        type=Path,
        default=root / "gin_out/p53_rh1_cell_spn_pipeline",
    )
    ap.add_argument("--chip-root", type=Path, default=default_chip_align_root())
    ap.add_argument("--node-list-dir", type=Path, default=DEFAULT_NODE_LIST_DIR)
    ap.add_argument(
        "--chromhmm-emission-tsv",
        type=Path,
        default=root / "runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv",
    )
    ap.add_argument("--runs", type=int, default=0, help="MC runs (0 = export only)")
    ap.add_argument("--blocks", type=int, default=10)
    ap.add_argument("--reference-qc", type=Path, default=None, help="Compare petri counts to this QC JSON")
    args = ap.parse_args()

    context_edges = args.context_edges.resolve()
    if not context_edges.is_file():
        raise SystemExit(f"missing context edges: {context_edges}")

    edges = load_subgraph_edges_tsv(context_edges)
    tf_symbols = load_tf_symbol_pool(args.node_list_dir)
    seed_genes = load_seed_genes(args.seed_file)

    cell_dir = chip_align_cell_dir(args.chip_root, args.cell)
    allowlist = cell_dir / "allowlist_union.tsv"
    tf_gate = cell_dir / "tf_gate_tfs.txt"
    peaks = cell_dir / "peaks"
    for p in (allowlist, tf_gate):
        if not p.is_file():
            raise SystemExit(f"missing {p}")

    measured = load_measured_tfs(tf_gate)
    prune_result = prune_grn_allowlist(
        edges,
        cell=args.cell,
        source_edges=str(context_edges),
        allowlist_path=allowlist,
        peaks_dir=str(peaks),
        measured_tfs=measured,
        tf_symbols=tf_symbols,
        seed_genes=seed_genes,
        keep_seed_tf_to_seed_gene=False,
    )
    out_cell = args.out_root / args.cell
    write_pruned_cell_dir(out_cell, prune_result)

    manifest_src = args.context_dir / "seed_manifest.json"
    if manifest_src.is_file():
        shutil.copy2(manifest_src, out_cell / "seed_manifest.json")

    pipeline_meta = {
        "pipeline": "realhop1_grn_prune_spn_export_v1",
        "context_edges": str(context_edges),
        "context_dir": str(args.context_dir),
        "seed_file": str(args.seed_file),
        "cell": args.cell,
        "prune_qc": prune_result.to_qc_dict(),
    }
    args.out_root.mkdir(parents=True, exist_ok=True)
    (out_cell / "pipeline_meta.json").write_text(
        json.dumps(pipeline_meta, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(
        f"[prune] {args.cell} edges {prune_result.n_edges_in} -> {prune_result.n_edges_out} "
        f"nodes {prune_result.n_nodes_in} -> {prune_result.n_nodes_out}",
        flush=True,
    )

    if not args.chromhmm_emission_tsv.is_file():
        raise SystemExit(f"missing chromHMM TSV: {args.chromhmm_emission_tsv}")

    export_runs = max(1, args.runs) if args.runs <= 0 else args.runs
    sim_kw = standard_h1_b0_simulate_kwargs(
        chromhmm_emission_tsv=str(args.chromhmm_emission_tsv),
        blocks=args.blocks if args.runs <= 0 else args.blocks,
        runs=export_runs,
        n_jobs=1,
        sim_subdir="spn_sim_repro",
        stem=f"p53_rh1_{args.cell}_export",
    )
    export_and_simulate_subgraph_dir(
        out_cell,
        keep_trajectories=args.runs > 0,
        **sim_kw,
    )

    qc_path = out_cell / "spn_export_qc.json"
    if not qc_path.is_file():
        raise SystemExit("SPN export failed: no spn_export_qc.json")

    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    pn = qc.get("petri_net", {})
    print(
        f"[spn] places={pn.get('n_places')} transitions={pn.get('n_transitions')} "
        f"arcs={pn.get('n_arcs')} graphml={out_cell / 'spn.graphml'}",
        flush=True,
    )

    if args.reference_qc:
        ref = json.loads(Path(args.reference_qc).read_text(encoding="utf-8"))
        rpn = ref.get("petri_net", {})
        diff = {
            k: (pn.get(k), rpn.get(k))
            for k in ("n_places", "n_transitions", "n_arcs")
            if pn.get(k) != rpn.get(k)
        }
        pipeline_meta["reference_qc"] = str(args.reference_qc)
        pipeline_meta["spn_compare"] = {"match": not diff, "diff": diff}
        (out_cell / "pipeline_meta.json").write_text(
            json.dumps(pipeline_meta, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"[compare] match={not diff} diff={diff}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

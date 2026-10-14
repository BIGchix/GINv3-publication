#!/usr/bin/env python3
"""Laptop-scale SigPN demo: subgraph_edges.tsv → export + short Monte-Carlo.

SigPN = signaling Petri net (not stochastic PN).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--subgraph-dir",
        type=Path,
        default=root / "gin_out/example_sigpn_demo",
        help="Directory containing subgraph_edges.tsv",
    )
    ap.add_argument(
        "--place-levels",
        type=Path,
        default=None,
        help="Optional place\\tlevel TSV (default: <subgraph-dir>/place_levels.tsv)",
    )
    ap.add_argument(
        "--chromhmm-emission-tsv",
        type=Path,
        default=root
        / "runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv",
    )
    ap.add_argument("--runs", type=int, default=2)
    ap.add_argument("--blocks", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-jobs", type=int, default=1)
    args = ap.parse_args()

    from gin_tool.subgraph_pipeline.spn_simulate import (
        export_and_simulate_subgraph_dir,
        load_place_levels_tsv,
    )
    from gin_tool.subgraph_pipeline.spn_standard_policy import (
        standard_h1_b0_simulate_kwargs,
    )

    sub = args.subgraph_dir.resolve()
    edges = sub / "subgraph_edges.tsv"
    if not edges.is_file():
        raise SystemExit(f"missing {edges}")
    if not args.chromhmm_emission_tsv.is_file():
        raise SystemExit(f"missing {args.chromhmm_emission_tsv}")

    pl_path = args.place_levels
    if pl_path is None:
        cand = sub / "place_levels.tsv"
        pl_path = cand if cand.is_file() else None
    place_levels = load_place_levels_tsv(pl_path) if pl_path else None

    kw = standard_h1_b0_simulate_kwargs(
        chromhmm_emission_tsv=str(args.chromhmm_emission_tsv),
        blocks=args.blocks,
        runs=args.runs,
        n_jobs=args.n_jobs,
        seed=args.seed,
        sim_subdir="spn_sim_demo",
        stem="sigpn_demo",
        place_levels=place_levels,
    )
    summary = export_and_simulate_subgraph_dir(sub, keep_trajectories=True, **kw)

    qc_path = sub / "spn_export_qc.json"
    if qc_path.is_file():
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        pn = qc.get("petri_net", {})
        print(
            f"[demo] places={pn.get('n_places')} transitions={pn.get('n_transitions')} "
            f"arcs={pn.get('n_arcs')}",
            flush=True,
        )
    print(f"[demo] summary keys: {sorted(summary)[:12]}…", flush=True)
    print(f"[demo] subgraph_dir={sub}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

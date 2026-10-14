#!/usr/bin/env python3
"""Summarize Fig3-style context subgraph hop stats (total/real nodes + edges)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def is_virtual(name: str) -> bool:
    return ";" in name


def summarize_edges(path: Path) -> dict:
    nodes: set[str] = set()
    n_edges = 0
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            nodes.add(a)
            nodes.add(b)
            n_edges += 1
    real = {n for n in nodes if not is_virtual(n)}
    virt = nodes - real
    return {
        "n_total_nodes": len(nodes),
        "n_real_nodes": len(real),
        "n_virtual_nodes": len(virt),
        "n_edges": n_edges,
    }


def summarize_root(root: Path) -> dict:
    edges = root / "subgraph_edges.tsv"
    if not edges.is_file():
        raise FileNotFoundError(edges)
    stats = summarize_edges(edges)
    meta_path = root / "run_meta.json"
    meta: dict = {}
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    ctx = meta.get("context_subgraph") or {}
    name = root.name
    # parse realhop / gh from dirname if present
    real_hop = None
    graph_hop = ctx.get("mech_max_depth")
    if "realhop" in name:
        try:
            real_hop = int(name.split("realhop")[1].split("_")[0])
        except (IndexError, ValueError):
            real_hop = None
    if "gh" in name and graph_hop is None:
        try:
            graph_hop = int(name.split("gh")[1].split("_")[0])
        except (IndexError, ValueError):
            pass
    return {
        "root": str(root),
        "label": name,
        "real_hop": real_hop,
        "graph_hop": graph_hop,
        "context_pipeline_version": meta.get("context_pipeline_version"),
        "subgraph_mode": meta.get("subgraph_mode"),
        **stats,
        "notes": (
            "Shared pathway context only: no cell-specific GRN allowlist prune; "
            "no chromHMM (simulation-time). Virtual = node name contains ';'."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--roots", nargs="+", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    rows = [summarize_root(Path(r)) for r in args.roots]
    rows.sort(key=lambda r: (r.get("real_hop") is None, r.get("real_hop") or 0))
    payload = {
        "figure": "fig3_panel_E_style",
        "definition": {
            "real_hop": "R→V→R on bipartite GIN",
            "graph_hop": "directed BFS depth (--context-mech-depth); graph-hop=2 ≡ real-hop=1",
            "excludes": ["GRN allowlist prune", "chromHMM emission weights"],
            "includes_v4_topology": [
                "dogma chain repair",
                "HME soft-seed + histprune",
                "drop non-seed promoter.active/;state",
                "LCC bridge",
                "met-enzyme recovery",
            ],
        },
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print(f"{'label':40} {'rh':>3} {'gh':>3} {'total':>8} {'real':>8} {'virt':>8} {'edges':>8}")
    for r in rows:
        print(
            f"{r['label'][:40]:40} "
            f"{str(r.get('real_hop') or '-'):>3} "
            f"{str(r.get('graph_hop') or '-'):>3} "
            f"{r['n_total_nodes']:8d} "
            f"{r['n_real_nodes']:8d} "
            f"{r['n_virtual_nodes']:8d} "
            f"{r['n_edges']:8d}"
        )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Validate that this publication bundle has required paths and lock fingerprints.

Does not run Stage A or Monte-Carlo. Safe for CI / reviewer smoke checks.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REQUIRED_PATHS = [
    "data/GINv3.11a.sifi.gz",
    "data/GINv3.11a.sifi.gz.sha256",
    "gin_tool/cli.py",
    "scbionet/spn",
    "nodeList/nodeList_protein.txt",
    "seedLists/p53_pi3k_curated/genes.txt",
    "seedLists/p53_pi3k_curated/metabolites.txt",
    "seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt",
    "scripts/reproduce_publication.py",
    "scripts/build_p53_pi3k_comma_seeds.sh",
    "docs/locks/p53_pub_gin11a_hs0_20260929.json",
    "expected/gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv",
    "expected/gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/subgraph_edges.tsv",
    "runtime/chip_gin_align/h1_wa01/allowlist_union.tsv",
    "runtime/omics_seed/cellline_levels/h1_wa01/h1_wa01_rna_plus_met_sources_place_levels.tsv",
]

EDGE_LOCKS = {
    "expected/gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv": (
        "350b7cb9acc8eaa7b1454f3c912f5aa1d94742a5f05b7492ed01b447039e6a43",
        11368,
    ),
    "expected/gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/subgraph_edges.tsv": (
        "562d79bc65749d06ea1291a33f12c6f774d6135e027d6e717438a3e406b0f15a",
        11363,
    ),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    errors: list[str] = []
    for rel in REQUIRED_PATHS:
        p = ROOT / rel
        if not p.exists():
            errors.append(f"missing {rel}")

    for rel, (exp_sha, exp_n) in EDGE_LOCKS.items():
        p = ROOT / rel
        if not p.is_file():
            continue
        got = sha256_file(p)
        n = sum(1 for _ in p.open())
        if got != exp_sha:
            errors.append(f"sha mismatch {rel}: got {got[:16]}…")
        if n != exp_n:
            errors.append(f"line count {rel}: got {n} expected {exp_n}")

    lock = ROOT / "docs/locks/p53_pub_gin11a_hs0_20260929.json"
    if lock.is_file():
        meta = json.loads(lock.read_text())
        sa = meta.get("stage_a") or {}
        if sa.get("sha256_edges") != EDGE_LOCKS[
            "expected/gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv"
        ][0]:
            errors.append("lock JSON stage_a sha does not match EDGE_LOCKS table")

    if errors:
        print("FAIL validate_bundle:")
        for e in errors:
            print(" -", e)
        return 1
    print("OK validate_bundle: required paths present; Stage A / H1 fingerprints match lock")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

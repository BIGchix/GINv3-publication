#!/usr/bin/env python3
"""Pipeline active-TF selection (step 1 before context subgraph).

Default null = **B / cohesion** (fast enough for routine runs).
Optional **C / rewire** via ``--method c`` (slower, stricter / more accurate).

Writes:
  {out-dir}/active_tfs_emp05.txt
  {out-dir}/active_tfs_emp05.meta.json
  plus underlying ``active_tf_rwr_*_{summary,full,emp05}.*`` under out-dir.

Example (default B):
  python scripts/select_active_tfs.py \\
    --gene-file seedLists/p53_pi3k_curated/genes.txt \\
    --out-dir gin_out/p53_active_tf_B

Optional C:
  python scripts/select_active_tfs.py \\
    --gene-file seedLists/p53_pi3k_curated/genes.txt \\
    --out-dir gin_out/p53_active_tf_C --method c --n-jobs 16
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


_METHOD_TO_STRENGTH = {
    "b": "cohesion",
    "cohesion": "cohesion",
    "rwr_cohesion": "cohesion",
    "c": "rewire",
    "rewire": "rewire",
    "rwr_rewire": "rewire",
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-file", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument(
        "--method",
        choices=sorted(_METHOD_TO_STRENGTH),
        default="b",
        help="Active-TF null: b/cohesion (default) or c/rewire (slower, stricter)",
    )
    ap.add_argument("--n-rwr-null", type=int, default=1000)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--cpus", dest="n_jobs", type=int, help="Alias for --n-jobs")
    ap.add_argument("--emp-p-cutoff", type=float, default=0.05)
    ap.add_argument("--rewire-sweeps", type=float, default=0.1)
    ap.add_argument("--real-hops", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    strength = _METHOD_TO_STRENGTH[args.method]
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tag = "B_cohesion" if strength == "cohesion" else "C_rewire"
    prefix = out / f"active_tf_rwr_{tag}"

    script = Path(__file__).resolve().parent / "active_tf_rwr_null_strength.py"
    cmd = [
        sys.executable,
        "-u",
        str(script),
        "--gene-file",
        str(args.gene_file),
        "--null-strength",
        strength,
        "--n-rwr-null",
        str(args.n_rwr_null),
        "--n-jobs",
        str(args.n_jobs),
        "--emp-p-cutoff",
        str(args.emp_p_cutoff),
        "--real-hops",
        str(args.real_hops),
        "--seed",
        str(args.seed),
        "--out-prefix",
        str(prefix),
    ]
    if strength == "rewire":
        cmd.extend(["--rewire-sweeps", str(args.rewire_sweeps)])

    print(f"[select_active_tfs] method={args.method} null_strength={strength}", flush=True)
    print(" ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        return rc

    emp_src = Path(str(prefix) + "_emp05.txt")
    meta_src = Path(str(prefix) + "_emp05.meta.json")
    emp_dst = out / "active_tfs_emp05.txt"
    meta_dst = out / "active_tfs_emp05.meta.json"
    shutil.copyfile(emp_src, emp_dst)
    meta = json.loads(meta_src.read_text(encoding="utf-8"))
    meta["pipeline_entry"] = "scripts/select_active_tfs.py"
    meta["method_arg"] = args.method
    meta["active_tfs_file"] = str(emp_dst)
    meta_dst.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    n = sum(1 for ln in emp_dst.read_text(encoding="utf-8").splitlines() if ln.strip())
    print(f"[select_active_tfs] wrote {emp_dst} n={n}", flush=True)
    print(
        f"[select_active_tfs] next: gin_tool.cli --active-tfs-file {emp_dst} ...",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

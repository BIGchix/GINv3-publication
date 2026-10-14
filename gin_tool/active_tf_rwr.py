"""Pipeline-integrated directed-RWR active-TF nulls (B/cohesion, C/rewire)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence

_STRENGTH = {
    "rwr_cohesion": "cohesion",
    "cohesion": "cohesion",
    "b": "cohesion",
    "rwr_rewire": "rewire",
    "rewire": "rewire",
    "c": "rewire",
}


def resolve_null_strength(method: str) -> str:
    key = method.strip().lower()
    if key not in _STRENGTH:
        raise ValueError(
            f"unknown RWR active-TF method {method!r}; "
            "expected rwr_cohesion/b or rwr_rewire/c"
        )
    return _STRENGTH[key]


def run_rwr_null_active_tfs(
    gene_symbols: Sequence[str],
    *,
    method: str,
    work_dir: Path,
    n_rwr_null: int = 1000,
    n_jobs: int = 8,
    emp_p_cutoff: float = 0.05,
    rewire_sweeps: float = 0.1,
    real_hops: int = 1,
    seed: int = 42,
) -> tuple[list[dict], dict]:
    """
    Run directed RWR null (B/cohesion or C/rewire) and return emp_p-passing records.

    Returns (records_sorted, summary) where each record has symbol / rwr_score /
    rwr_emp_p / rwr_emp_q / rwr_ge_count.
    """
    strength = resolve_null_strength(method)
    genes = sorted({g.strip() for g in gene_symbols if g and g.strip()})
    if not genes:
        raise ValueError("RWR active-TF selection needs non-empty gene symbols")

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    gene_file = work_dir / "rwr_active_tf_genes.txt"
    gene_file.write_text("\n".join(genes) + "\n", encoding="utf-8")

    tag = "B_cohesion" if strength == "cohesion" else "C_rewire"
    prefix = work_dir / f"active_tf_rwr_{tag}"
    script = Path(__file__).resolve().parents[1] / "scripts" / "active_tf_rwr_null_strength.py"
    if not script.is_file():
        raise FileNotFoundError(f"missing {script}")

    cmd = [
        sys.executable,
        "-u",
        str(script),
        "--gene-file",
        str(gene_file),
        "--null-strength",
        strength,
        "--n-rwr-null",
        str(n_rwr_null),
        "--n-jobs",
        str(n_jobs),
        "--emp-p-cutoff",
        str(emp_p_cutoff),
        "--real-hops",
        str(real_hops),
        "--seed",
        str(seed),
        "--out-prefix",
        str(prefix),
    ]
    if strength == "rewire":
        cmd.extend(["--rewire-sweeps", str(rewire_sweeps)])

    print(f"[active_tf_rwr] method={method} strength={strength} n_genes={len(genes)}", flush=True)
    print(" ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        raise RuntimeError(f"active_tf_rwr_null_strength failed rc={rc}")

    full_path = Path(str(prefix) + "_full.json")
    payload = json.loads(full_path.read_text(encoding="utf-8"))
    summary = payload.get("summary") or {}
    records = payload.get("records") or []
    ordered = sorted(records, key=lambda r: (r["rwr_emp_p"], -r["rwr_score"]))
    selected = [r for r in ordered if r["rwr_emp_p"] <= emp_p_cutoff]

    emp_dst = work_dir / "active_tfs_emp05.txt"
    emp_dst.write_text(
        "\n".join(r["symbol"] for r in selected) + ("\n" if selected else ""),
        encoding="utf-8",
    )
    meta = {
        "null_strength": strength,
        "method": method,
        "emp_p_cutoff": emp_p_cutoff,
        "n_active": len(selected),
        "source_full": str(full_path),
        "active_tfs_file": str(emp_dst),
    }
    (work_dir / "active_tfs_emp05.meta.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    print(f"[active_tf_rwr] selected n={len(selected)} -> {emp_dst}", flush=True)
    return selected, summary

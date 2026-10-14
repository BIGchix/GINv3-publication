#!/usr/bin/env python3
"""Assert on-disk Fig4 C/D noprune artifacts match the publication freeze lock."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = ROOT / "docs/locks/p53_pub_fig4cd_noprune_gin11a_20261001.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    args = ap.parse_args()
    lock_path = args.lock if args.lock.is_absolute() else ROOT / args.lock
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    errors: list[str] = []

    checks: list[tuple[str, str, str]] = [
        ("stage_a_edges", lock["stage_a"]["edges_tsv"], lock["stage_a"]["sha256_edges"]),
        ("work_edges", lock["work_dir"]["edges_tsv"], lock["work_dir"]["sha256_edges"]),
        ("corr_json", lock["work_dir"]["corr_json"], lock["work_dir"]["sha256_corr"]),
        ("fig_C_pdf", lock["panels"]["C"]["pdf"], lock["panels"]["C"]["sha256_pdf"]),
        ("fig_C_png", lock["panels"]["C"]["png"], lock["panels"]["C"]["sha256_png"]),
        ("fig_D_pdf", lock["panels"]["D"]["pdf"], lock["panels"]["D"]["sha256_pdf"]),
        ("fig_D_png", lock["panels"]["D"]["png"], lock["panels"]["D"]["sha256_png"]),
    ]
    for lab, info in lock["place_levels"].items():
        checks.append((f"levels_{lab}", info["path"], info["sha256"]))
    for lab, info in lock["trajectories_metsrc"].items():
        checks.append((f"traj_{lab}", info["path"], info["sha256"]))

    for key, rel, exp in checks:
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"missing {rel}")
            continue
        got = sha256(path)
        if got != exp:
            errors.append(f"{key} sha mismatch: got={got} expected={exp}")

    if not lock["work_dir"].get("edges_match_stage_a", True):
        errors.append("lock claims edges_match_stage_a=false")
    else:
        a = ROOT / lock["stage_a"]["edges_tsv"]
        w = ROOT / lock["work_dir"]["edges_tsv"]
        if a.is_file() and w.is_file() and sha256(a) != sha256(w):
            errors.append("work_dir edges no longer match Stage A")

    qc_path = ROOT / "gin_out/p53_fig4_gin11a_noprune_cd/spn_export_qc.json"
    if qc_path.is_file():
        qc = json.loads(qc_path.read_text())
        exp_places = lock["work_dir"]["spn"]["n_places"]
        got_places = qc["petri_net"]["n_places"]
        if got_places != exp_places:
            errors.append(f"n_places={got_places} expected={exp_places}")
        if qc.get("spn_pipeline_version") != lock["pipeline"]["spn_pipeline_version"]:
            errors.append("spn_pipeline_version mismatch vs lock")
        if qc.get("maintain_seeds") is not True:
            errors.append("maintain_seeds expected True for noprune Fig4 C/D")
    else:
        errors.append(f"missing {qc_path.relative_to(ROOT)}")

    backup = ROOT / lock["backup"]["snapshot_dir"]
    if not backup.is_dir():
        errors.append(f"missing backup snapshot {lock['backup']['snapshot_dir']}")
    bundle = ROOT / lock["backup"]["bundle_tar_gz"]
    if not bundle.is_file():
        errors.append(f"missing bundle {lock['backup']['bundle_tar_gz']}")
    elif sha256(bundle) != lock["backup"]["bundle_sha256"]:
        errors.append("backup bundle sha mismatch")

    if errors:
        print("FAIL fig4cd noprune lock:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1

    print(
        json.dumps(
            {
                "ok": True,
                "version_id": lock["version_id"],
                "status": lock["status"],
                "role": lock["role"],
                "lock": str(lock_path.relative_to(ROOT)),
                "panels": {
                    "C": {
                        "vmin": lock["panels"]["C"]["vmin"],
                        "offdiag": lock["panels"]["C"]["offdiag_range"],
                    },
                    "D": {
                        "vmin": lock["panels"]["D"]["vmin"],
                        "offdiag": lock["panels"]["D"]["offdiag_range"],
                    },
                },
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

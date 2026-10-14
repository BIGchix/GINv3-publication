#!/usr/bin/env python3
"""Assert on-disk artifacts match the frozen p53 GINv3.11a publication lock."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Prefer hash-seed refreeze lock; fall back to first freeze if absent.
_HS0 = ROOT / "docs/locks/p53_pub_gin11a_hs0_20260929.json"
_V1 = ROOT / "docs/locks/p53_pub_gin11a_20260929.json"
DEFAULT_LOCK = _HS0 if _HS0.is_file() else _V1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    ap.add_argument(
        "--skip-sifi",
        action="store_true",
        help="Skip hashing GINv3.11a.sifi (large); still check edges",
    )
    args = ap.parse_args()
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    errors: list[str] = []

    checks: list[tuple[str, str, bool]] = [
        ("stage_a", lock["stage_a"]["edges_tsv"], False),
        ("stage_b_h1", lock["stage_b_h1"]["edges_tsv"], False),
    ]
    if not args.skip_sifi:
        checks.insert(0, ("gin", lock["gin"]["sifi"], True))

    expected_sha = {
        "gin": lock["gin"]["sha256"],
        "stage_a": lock["stage_a"]["sha256_edges"],
        "stage_b_h1": lock["stage_b_h1"]["sha256_edges"],
    }
    expected_n = {
        "stage_a": lock["stage_a"]["n_edges"],
        "stage_b_h1": lock["stage_b_h1"]["n_edges"],
    }

    for key, rel, _ in checks:
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"missing {rel}")
            continue
        got = sha256(path)
        exp = expected_sha[key]
        if got != exp:
            errors.append(f"{key} sha mismatch: got={got} expected={exp}")
        if key in expected_n:
            n = sum(1 for line in path.open() if line.strip())
            if n != expected_n[key]:
                errors.append(f"{key} n_edges={n} expected={expected_n[key]}")

    places_exp = lock["stage_b_h1"]["spn"]["n_places"]
    qc_path = ROOT / "gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/spn_export_qc.json"
    if qc_path.is_file():
        places = json.loads(qc_path.read_text())["petri_net"]["n_places"]
        if places != places_exp:
            errors.append(f"n_places={places} expected={places_exp}")
    else:
        errors.append(f"missing {qc_path.relative_to(ROOT)}")

    if errors:
        print("FAIL p53_pub_gin11a lock:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1

    print(
        json.dumps(
            {
                "ok": True,
                "version_id": lock["version_id"],
                "publication_rng_seed": lock["stage_c_pulse"]["publication_rng_seed"],
                "lock": str(args.lock.relative_to(ROOT)),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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

# Uncompressed GINv3.11a.sifi (lock JSON) and bundled gzip (sidecar).
LOCK_SIFI_SHA = "192b2cb5a8147324fc41064a3f20c943a3d82c03972cff544e3d26219e7ace83"
LOCK_SIFI_GZ_SHA = "f9bb23c44c30d0dd6eab73225600085b378d1912e59f1ff0d8a3a0b37cfd0a58"
SIFI_PLAIN = ROOT / "GINv3.11a.sifi"
SIFI_GZ = ROOT / "data" / "GINv3.11a.sifi.gz"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_gin_prior(errors: list[str], lock_sha: str) -> None:
    """Default: verify bundled ``.sifi.gz``; else plain ``.sifi``."""
    if SIFI_GZ.is_file():
        got = sha256(SIFI_GZ)
        if got != LOCK_SIFI_GZ_SHA:
            errors.append(
                f"gin.gz sha mismatch: got={got} expected={LOCK_SIFI_GZ_SHA}"
            )
        return
    if SIFI_PLAIN.is_file():
        got = sha256(SIFI_PLAIN)
        exp = lock_sha or LOCK_SIFI_SHA
        if got != exp:
            errors.append(f"gin sha mismatch: got={got} expected={exp}")
        return
    errors.append(f"missing {SIFI_GZ.relative_to(ROOT)} or {SIFI_PLAIN.name}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    ap.add_argument(
        "--skip-sifi",
        action="store_true",
        help="Skip hashing GIN prior (plain or .gz); still check edges",
    )
    args = ap.parse_args()
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    errors: list[str] = []

    if not args.skip_sifi:
        check_gin_prior(errors, lock["gin"]["sha256"])

    checks: list[tuple[str, str]] = [
        ("stage_a", lock["stage_a"]["edges_tsv"]),
        ("stage_b_h1", lock["stage_b_h1"]["edges_tsv"]),
    ]
    expected_sha = {
        "stage_a": lock["stage_a"]["sha256_edges"],
        "stage_b_h1": lock["stage_b_h1"]["sha256_edges"],
    }
    expected_n = {
        "stage_a": lock["stage_a"]["n_edges"],
        "stage_b_h1": lock["stage_b_h1"]["n_edges"],
    }

    for key, rel in checks:
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"missing {rel}")
            continue
        got = sha256(path)
        exp = expected_sha[key]
        if got != exp:
            errors.append(f"{key} sha mismatch: got={got} expected={exp}")
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
                "gin_prior": (
                    str(SIFI_GZ.relative_to(ROOT))
                    if SIFI_GZ.is_file()
                    else str(SIFI_PLAIN.relative_to(ROOT))
                ),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

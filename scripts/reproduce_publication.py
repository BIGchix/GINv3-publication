#!/usr/bin/env python3
"""Publication recipe: RECOMPUTE manuscript artifacts with gin_tool / SigPN.

This driver calls the same CLI path as the freeze (it does not copy results
into place as the primary workflow). For custom biology, use
``python -m gin_tool.cli`` — see docs/USER_GUIDE.md and examples/.

Steps (lock p53_pub_gin11a_hs0_20260929 + fig4cd 20261001):

  prepare   unpack data/GINv3.11a.sifi.gz and check SHA-256
  stage_a   recompute context v4 · real-hop 1 · no GRN prune
  fig3e     recompute real-hop 1/2/3 series
  fig4cd    recompute six-line noprune SigPN (200×200)
  prune_h1  H1 ChIP allowlist prune (−5 edges)
  pulse     H1 B0 pulse seed=43 · 350×350
  assert    compare on-disk artifacts to publication locks

``--from-expected`` only installs reference fingerprints for CI / packaging
checks — it is NOT the reviewer reproduction path.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIFI = ROOT / "GINv3.11a.sifi"
SIFI_GZ = ROOT / "data" / "GINv3.11a.sifi.gz"
LOCK_SIFI_SHA = "192b2cb5a8147324fc41064a3f20c943a3d82c03972cff544e3d26219e7ace83"
LOCK_SIFI_BYTES = 824892563

ALL_STEPS = ("prepare", "stage_a", "fig3e", "fig4cd", "prune_h1", "pulse", "assert")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str], *, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    merged = os.environ.copy()
    merged["PYTHONHASHSEED"] = merged.get("PYTHONHASHSEED", "0")
    merged["PYTHONPATH"] = (
        str(ROOT) + os.pathsep + merged["PYTHONPATH"]
        if merged.get("PYTHONPATH")
        else str(ROOT)
    )
    if env:
        merged.update(env)
    subprocess.check_call(cmd, cwd=str(ROOT), env=merged)


def seeds_csv() -> str:
    out = subprocess.check_output(
        ["bash", str(ROOT / "scripts" / "build_p53_pi3k_comma_seeds.sh")],
        cwd=str(ROOT),
        text=True,
    )
    return out.strip()


def step_prepare() -> None:
    if SIFI.is_file() and SIFI.stat().st_size == LOCK_SIFI_BYTES:
        got = sha256_file(SIFI)
        if got != LOCK_SIFI_SHA:
            raise SystemExit(f"GINv3.11a.sifi sha mismatch: {got}")
        print(f"OK existing {SIFI} ({LOCK_SIFI_BYTES} bytes)")
        return
    if not SIFI_GZ.is_file():
        raise SystemExit(f"missing {SIFI_GZ}")
    print(f"unpack {SIFI_GZ} → {SIFI}")
    with gzip.open(SIFI_GZ, "rb") as src, SIFI.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    if SIFI.stat().st_size != LOCK_SIFI_BYTES:
        raise SystemExit(f"unpacked size {SIFI.stat().st_size} != {LOCK_SIFI_BYTES}")
    got = sha256_file(SIFI)
    if got != LOCK_SIFI_SHA:
        raise SystemExit(f"unpacked sha mismatch: {got}")
    print(f"OK {SIFI} sha256={got}")


def _context_cmd(out: Path, *, mech_depth: int, n_jobs: int) -> list[str]:
    atf = ROOT / "seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt"
    return [
        sys.executable,
        "-u",
        "-m",
        "gin_tool.cli",
        "--sifi",
        str(SIFI),
        "--seeds",
        seeds_csv(),
        "--out-dir",
        str(out),
        "--subgraph-mode",
        "context",
        "--context-pipeline-version",
        "v4",
        "--context-mech-depth",
        str(mech_depth),
        "--context-n-jobs",
        str(n_jobs),
        "--active-tfs-file",
        str(atf),
        "--no-debias",
    ]


def step_stage_a(*, n_jobs: int) -> None:
    out = ROOT / "gin_out/p53_shared_context_v4_rh1_gin11a"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    run(_context_cmd(out, mech_depth=2, n_jobs=n_jobs))


def step_fig3e(*, n_jobs: int) -> None:
    mapping = ((1, 2), (2, 4), (3, 6))
    roots = []
    for rh, gh in mapping:
        out = ROOT / f"gin_out/p53_fig3e_v4_gin11a_realhop{rh}_gh{gh}"
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        run(_context_cmd(out, mech_depth=gh, n_jobs=n_jobs))
        roots.append(str(out.relative_to(ROOT)))
    stats = ROOT / "benchmarks/subgraph-speed/fig3e_v4_gin11a_context_hop_stats.json"
    stats.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            sys.executable,
            "-u",
            "scripts/summarize_fig3_context_hops.py",
            "--roots",
            *roots,
            "--out",
            str(stats),
        ]
    )
    plot = ROOT / "scripts/plot_fig3e_v4_hops.py"
    if plot.is_file():
        stem = ROOT / "manuscript/figures/fig3/fig3_E_v4_gin11a_hops"
        stem.parent.mkdir(parents=True, exist_ok=True)
        try:
            run([sys.executable, "-u", str(plot), "--stats", str(stats), "--out-stem", str(stem)])
        except subprocess.CalledProcessError as exc:
            print(f"warning: fig3e plot skipped ({exc})", flush=True)


def step_fig4cd(*, n_jobs: int, runs: int, blocks: int) -> None:
    cmd = [
        sys.executable,
        "-u",
        "scripts/run_fig4cd_noprune_gin11a_stageA.py",
        "--n-jobs",
        str(n_jobs),
    ]
    # publication runner may not expose --runs; pass through if supported
    help_txt = subprocess.check_output(
        [sys.executable, "scripts/run_fig4cd_noprune_gin11a_stageA.py", "-h"],
        cwd=str(ROOT),
        text=True,
    )
    if "--runs" in help_txt:
        cmd.extend(["--runs", str(runs), "--blocks", str(blocks)])
    run(cmd)


def step_prune_h1() -> None:
    shared = ROOT / "gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv"
    if not shared.is_file():
        raise SystemExit("run stage_a first")
    cell_root = ROOT / "gin_out/p53_stageA_gin11a_h1_grn_prune"
    cell = cell_root / "h1_wa01"
    if cell.exists():
        shutil.rmtree(cell)
    cmd = [
        sys.executable,
        "-u",
        "scripts/build_realhop1_cell_spn.py",
        "--context-edges",
        str(shared),
        "--context-dir",
        str(shared.parent),
        "--out-root",
        str(cell_root),
        "--cell",
        "h1_wa01",
        "--runs",
        "0",
        "--node-list-dir",
        str(ROOT / "nodeList"),
        "--chip-root",
        str(ROOT / "runtime/chip_gin_align"),
        "--chromhmm-emission-tsv",
        str(ROOT / "runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv"),
    ]
    ref = ROOT / "expected/gin_out/p53_cell_grn_prune_v4_seedhist/h1_wa01/spn_export_qc.json"
    if not ref.is_file():
        ref = ROOT / "expected/gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/spn_export_qc.json"
    if ref.is_file():
        cmd.extend(["--reference-qc", str(ref)])
    run(cmd)


def step_pulse(*, n_jobs: int, runs: int, blocks: int, seed: int) -> None:
    cell = ROOT / "gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01"
    if not (cell / "subgraph_edges.tsv").is_file():
        raise SystemExit("run prune_h1 first")
    out = ROOT / "runtime/spn_phys5"
    out.mkdir(parents=True, exist_ok=True)
    tag = f"std_stageA_gin11a_s{seed}"
    run(
        [
            sys.executable,
            "-u",
            "scripts/run_p53_b0_std_pulse.py",
            "--force",
            "--subgraph-dir",
            str(cell),
            "--tag",
            tag,
            "--out-dir",
            str(out),
            "--seed",
            str(seed),
            "--runs",
            str(runs),
            "--n-jobs",
            str(n_jobs),
            "--hold",
            "40",
            "--onset",
            "150",
            "--blocks",
            str(blocks),
        ]
    )


def step_assert(*, skip_sifi: bool, topology_only: bool = False) -> None:
    cmd = [sys.executable, "-u", "scripts/assert_p53_pub_gin11a_lock.py"]
    if skip_sifi or not SIFI.is_file():
        cmd.append("--skip-sifi")
    run(cmd)
    if topology_only:
        print("skip fig4cd lock (topology-only / --from-expected)")
        return
    fig4_work = ROOT / "gin_out/p53_fig4_gin11a_noprune_cd"
    fig_c = ROOT / "manuscript/figures/fig4/fig4_C_gin11a_noprune_histone.pdf"
    traj = fig4_work / "spn_sim"
    if fig4_work.is_dir() and fig_c.is_file() and traj.is_dir():
        run([sys.executable, "-u", "scripts/assert_p53_pub_fig4cd_noprune_lock.py"])
    else:
        print("skip fig4cd lock (panels/trajectories not rebuilt)")


def _install_tree(src: Path, dst_root: Path, *, label: str) -> None:
    if not src.is_dir():
        return
    for path in src.rglob("*"):
        if path.is_file():
            target = dst_root / path.relative_to(src)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            print(f"{label} → {target.relative_to(ROOT)}")


def install_expected_fingerprints() -> None:
    gin_src = ROOT / "expected/gin_out"
    if not gin_src.is_dir():
        raise SystemExit("missing expected/gin_out")
    _install_tree(gin_src, ROOT / "gin_out", label="expected")
    # Pulse JSON reports (seed 43) for offline inspection / assert helpers
    _install_tree(
        ROOT / "expected/runtime/spn_phys5",
        ROOT / "runtime/spn_phys5",
        label="expected-pulse",
    )
    _install_tree(
        ROOT / "expected/benchmarks/subgraph-speed",
        ROOT / "benchmarks/subgraph-speed",
        label="expected-bench",
    )
    # Reference PDFs into manuscript/figures when present
    for panel in ("fig3", "fig4"):
        _install_tree(
            ROOT / "expected/figures" / panel,
            ROOT / "manuscript/figures" / panel,
            label="expected-fig",
        )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--steps",
        default="all",
        help="comma list or 'all': " + ",".join(ALL_STEPS),
    )
    ap.add_argument("--n-jobs", type=int, default=int(os.environ.get("N_JOBS", "4")))
    ap.add_argument("--seed", type=int, default=43, help="publication pulse RNG (lock=43)")
    ap.add_argument("--runs", type=int, default=350)
    ap.add_argument("--blocks", type=int, default=350)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="tiny MC sizes (does not match publication lock numbers)",
    )
    ap.add_argument(
        "--from-expected",
        action="store_true",
        help=(
            "CI/packaging only: copy expected/ fingerprints into gin_out and "
            "assert topology. Not a substitute for recomputing Stage A."
        ),
    )
    args = ap.parse_args(argv)

    runs, blocks = args.runs, args.blocks
    if args.smoke:
        runs, blocks = 2, 4
        print("SMOKE: runs=2 blocks=4 (not publication-scale)", flush=True)

    if args.from_expected:
        install_expected_fingerprints()
        # Fig. 4C/D full lock needs trajectory dumps (large); not bundled.
        # Topology / H1 / pulse-JSON fingerprint check is enough offline.
        step_assert(skip_sifi=not SIFI.is_file(), topology_only=True)
        return 0

    wanted = list(ALL_STEPS) if args.steps == "all" else [s.strip() for s in args.steps.split(",") if s.strip()]
    unknown = [s for s in wanted if s not in ALL_STEPS]
    if unknown:
        raise SystemExit(f"unknown steps: {unknown}")

    print(json.dumps({"root": str(ROOT), "steps": wanted, "n_jobs": args.n_jobs}, indent=2))
    for step in wanted:
        print(f"\n===== {step} =====", flush=True)
        if step == "prepare":
            step_prepare()
        elif step == "stage_a":
            step_stage_a(n_jobs=args.n_jobs)
        elif step == "fig3e":
            step_fig3e(n_jobs=args.n_jobs)
        elif step == "fig4cd":
            step_fig4cd(n_jobs=args.n_jobs, runs=min(runs, 200) if not args.smoke else runs, blocks=min(blocks, 200) if not args.smoke else blocks)
        elif step == "prune_h1":
            step_prune_h1()
        elif step == "pulse":
            step_pulse(n_jobs=args.n_jobs, runs=runs, blocks=blocks, seed=args.seed)
        elif step == "assert":
            step_assert(skip_sifi=not SIFI.is_file(), topology_only=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

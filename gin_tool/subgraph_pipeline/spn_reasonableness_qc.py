"""Post-SPN reasonableness QC (trajectory + structural gates)."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from gin_tool.spn_export import is_markable_seed_place, is_virtual_gin_node


def _place_layer(name: str) -> str:
    low = name.lower()
    if ";" in name:
        return "Complex"
    if ".promoter" in low:
        return "Promoter"
    if low.endswith(".rna"):
        return "RNA"
    if low.endswith(".protein"):
        return "Protein"
    if name.startswith(("SID", "CID", "CHEBI", "CAS", "UC_")):
        return "Metabolite"
    if name.startswith("chr"):
        return "Chromatin"
    return "Other"


def _pctiles(vals: np.ndarray, ps: Sequence[float]) -> dict[str, float]:
    if vals.size == 0:
        return {f"p{p}": float("nan") for p in ps} | {
            "min": float("nan"),
            "max": float("nan"),
            "mean": float("nan"),
            "median": float("nan"),
        }
    out = {
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
        "mean": float(np.mean(vals)),
        "median": float(np.median(vals)),
    }
    for p in ps:
        out[f"p{p}"] = float(np.percentile(vals, p))
    return out


def _safe_cv(std: np.ndarray, mean: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    m = np.asarray(mean, dtype=np.float64)
    s = np.asarray(std, dtype=np.float64)
    out = np.full_like(m, np.nan, dtype=np.float64)
    mask = np.abs(m) > eps
    out[mask] = s[mask] / m[mask]
    return out


@dataclass
class SpnReasonablenessQC:
    stem: str
    ok: bool
    warnings: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)
    report: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "stem": self.stem,
            "ok": self.ok,
            "warnings": list(self.warnings),
            "failures": list(self.failures),
            **self.report,
        }


def count_transcription_to_promoter_edges(
    edges: Iterable[tuple[str, str]],
) -> int:
    return sum(
        1
        for a, b in edges
        if is_virtual_gin_node(a)
        and "transcription" in a.lower()
        and ".promoter" in b.lower()
    )


def _load_edges_tsv(path: Path) -> list[tuple[str, str]]:
    edges: list[tuple[str, str]] = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            edges.append((parts[0], parts[1]))
    return edges


def qc_spn_reasonableness_from_arrays(
    *,
    stem: str,
    place_names: Sequence[str],
    mean_markings: np.ndarray,
    std_markings: np.ndarray,
    n_runs: int,
    n_blocks: int,
    n_marked_seeds: int,
    n_interactions: int,
    maintain_seeds: bool,
    marked_seed_names: Sequence[str] | None = None,
    mark_mode: str | None = None,
    warm_start_block: int = 10,
    runaway_growth_warn: float = 5.0,
    runaway_growth_fail: float = 50.0,
    subgraph_edges: list[tuple[str, str]] | None = None,
) -> SpnReasonablenessQC:
    """
    Reasonableness QC on mean/std trajectories (across MC runs).

    Structural gates (product decisions):
      - marked seeds are protein-only when mark_mode is protein-scoped
      - no transcription→promoter arcs in subgraph edges (if provided)
    Trajectory gates:
      - no NaN/Inf
      - no runaway growth (final/mean@warm > runaway_growth among high finals)
      - warm temporal CV not exploding
    """
    places = list(place_names)
    mean = np.asarray(mean_markings, dtype=np.float64)
    std = np.asarray(std_markings, dtype=np.float64)
    if mean.ndim != 2 or mean.shape[1] != len(places):
        raise ValueError(f"mean_markings shape {mean.shape} incompatible with n_places={len(places)}")
    if std.shape != mean.shape:
        raise ValueError("std_markings shape must match mean_markings")

    warnings: list[str] = []
    failures: list[str] = []

    if not np.isfinite(mean).all() or not np.isfinite(std).all():
        failures.append("non-finite values in mean/std markings")

    final = mean[-1]
    final_std = std[-1]
    active = final > 0
    always_zero = np.all(mean <= 0, axis=0)
    ever_active = np.any(mean > 0, axis=0)

    layers = np.array([_place_layer(p) for p in places])
    by_layer: dict[str, dict[str, float | int]] = {}
    for layer in sorted(set(layers.tolist())):
        idx = np.where(layers == layer)[0]
        fin = final[idx]
        act = fin > 0
        by_layer[layer] = {
            "n": int(idx.size),
            "active": int(np.sum(act)),
            "zero": int(np.sum(~act)),
            "median_active": float(np.median(fin[act])) if np.any(act) else float("nan"),
            "max": float(np.max(fin)) if idx.size else float("nan"),
        }

    active_final = final[active]
    dist = _pctiles(active_final, (90, 95, 99, 99.9))

    cv_runs = _safe_cv(final_std, final)
    cv_active = cv_runs[active & np.isfinite(cv_runs)]
    run_cv = {
        "median": float(np.median(cv_active)) if cv_active.size else float("nan"),
        "p90": float(np.percentile(cv_active, 90)) if cv_active.size else float("nan"),
        "p95": float(np.percentile(cv_active, 95)) if cv_active.size else float("nan"),
        "frac_cv_gt_0.5": float(np.mean(cv_active > 0.5)) if cv_active.size else float("nan"),
        "frac_cv_gt_1.0": float(np.mean(cv_active > 1.0)) if cv_active.size else float("nan"),
    }

    warm0 = min(max(warm_start_block, 0), max(mean.shape[0] - 1, 0))
    warm = mean[warm0:]
    warm_std = np.std(warm, axis=0, ddof=1) if warm.shape[0] > 1 else np.zeros(len(places))
    warm_mean = np.mean(warm, axis=0)
    tcv = _safe_cv(warm_std, warm_mean)
    tcv_active = tcv[active & np.isfinite(tcv)]
    temporal = {
        "warm_start_block": warm0,
        "median": float(np.median(tcv_active)) if tcv_active.size else float("nan"),
        "p90": float(np.percentile(tcv_active, 90)) if tcv_active.size else float("nan"),
        "frac_tcv_gt_0.5": float(np.mean(tcv_active > 0.5)) if tcv_active.size else float("nan"),
    }
    if tcv_active.size and float(np.median(tcv_active)) > 0.5:
        warnings.append(
            f"warm temporal CV median {float(np.median(tcv_active)):.3f} high — trajectory not plateaued"
        )

    # pileup / runaway (maintain_seeds can accumulate; warn early, fail only extreme)
    p99 = dist.get("p99", float("nan"))
    hi = final >= p99 if math.isfinite(p99) else np.zeros_like(final, dtype=bool)
    base = np.maximum(mean[warm0], 1e-12)
    growth = final / base
    runaway_warn = hi & (growth > runaway_growth_warn)
    runaway_fail = hi & (growth > runaway_growth_fail)
    n_runaway_warn = int(np.sum(runaway_warn))
    n_runaway_fail = int(np.sum(runaway_fail))
    if n_runaway_fail:
        failures.append(
            f"{n_runaway_fail} high-token places with growth>{runaway_growth_fail} from block{warm0}"
        )
    elif n_runaway_warn:
        warnings.append(
            f"{n_runaway_warn} high-token places with growth>{runaway_growth_warn} from block{warm0} "
            f"(accumulating under maintain_seeds; check hub proteins)"
        )

    order = np.argsort(-final)
    top_pileup = []
    for i in order[:20]:
        if final[i] <= 0:
            continue
        top_pileup.append(
            {
                "place": places[i],
                "layer": layers[i],
                "final_mean": float(final[i]),
                "final_std": float(final_std[i]),
                "cv_runs": float(cv_runs[i]) if math.isfinite(cv_runs[i]) else None,
                "tcv_warm": float(tcv[i]) if math.isfinite(tcv[i]) else None,
                "growth_b10": float(growth[i]),
            }
        )

    # stratified run-CV
    def _med_cv(mask: np.ndarray) -> float:
        v = cv_runs[mask & np.isfinite(cv_runs)]
        return float(np.median(v)) if v.size else float("nan")

    run_cv_strat = {
        "note": "run-CV inflated for low-count places; interpret by stratum",
        "mean_gt_1_median_cv": _med_cv(final > 1),
        "mean_gt_5_median_cv": _med_cv(final > 5),
        "mean_gt_10_median_cv": _med_cv(final > 10),
        "top20_median_cv": _med_cv(np.isin(np.arange(len(places)), order[:20])),
    }

    # structural: protein-only mark required only for protein-scoped mark modes
    marked = list(marked_seed_names or [])
    mode = (mark_mode or "seed_protein").strip().lower()
    protein_scoped = mode in ("seed_protein", "source_protein", "ancestor_protein")
    non_protein_marked = [m for m in marked if not is_markable_seed_place(m)]
    if protein_scoped and non_protein_marked:
        failures.append(
            f"{len(non_protein_marked)} marked seeds are not .protein "
            f"(sample: {non_protein_marked[:8]})"
        )
    if maintain_seeds and n_marked_seeds <= 0:
        failures.append("maintain_seeds=True but n_marked_seeds=0")

    n_trans_prom = None
    if subgraph_edges is not None:
        n_trans_prom = 0
        for a, b in subgraph_edges:
            if is_virtual_gin_node(a) and "transcription" in a.lower() and ".promoter" in b.lower():
                n_trans_prom += 1
        if n_trans_prom:
            failures.append(f"transcription→promoter edges in subgraph: {n_trans_prom}")

    if run_cv["median"] == run_cv["median"] and run_cv["median"] > 1.5:
        warnings.append(
            f"median run-CV {run_cv['median']:.3f} high — unstable MC (stochastic SPN may be OK)"
        )

    hints: list[str] = []
    frac_tok = float(np.mean(active)) if len(places) else 0.0
    if frac_tok > 0.95:
        hints.append("NOTE: nearly all places have tokens — expected if connected + maintain seeds")
    if temporal["median"] == temporal["median"] and temporal["median"] < 0.5:
        hints.append(f"OK: warm temporal CV median {temporal['median']:.3f} — roughly plateaued")
    if n_runaway_fail == 0 and n_runaway_warn == 0:
        hints.append(
            f"OK: no runaway pileup (growth>{runaway_growth_warn} among ≥p99); max final={dist.get('max')}"
        )
    dead_frac = float(np.mean(always_zero)) if len(places) else 0.0
    if dead_frac > 0.5:
        warnings.append(
            f"{dead_frac:.0%} places always-zero ({int(np.sum(always_zero))}/{len(places)}) — "
            "many histone/promoter marks may stay inert"
        )
    if not non_protein_marked:
        hints.append("OK: marked seeds are protein-only")
    if n_trans_prom == 0:
        hints.append("OK: transcription→promoter edges = 0")

    ok = not failures
    overall = "ok" if ok and not warnings else ("acceptable_with_warnings" if ok else "failed")
    verdict = {
        "overall": overall,
        "dead_places": f"{int(np.sum(always_zero))}/{len(places)} always-zero",
        "pileup": (
            "no runaway"
            if n_runaway_fail == 0 and n_runaway_warn == 0
            else (
                f"{n_runaway_fail} extreme / {n_runaway_warn} elevated growth≥p99"
            )
        ),
        "temporal_stability": (
            f"warm temporal CV median {temporal['median']:.3f}"
            if temporal["median"] == temporal["median"]
            else "n/a"
        ),
        "mc_noise": (
            f"run-to-run CV median {run_cv['median']:.3f}"
            if run_cv["median"] == run_cv["median"]
            else "n/a"
        ),
        "structure": {
            "mark_mode": mode,
            "marked_protein_only": (not bool(non_protein_marked)) if protein_scoped else None,
            "n_marked_seeds": n_marked_seeds,
            "n_non_protein_marked": len(non_protein_marked),
            "transcription_to_promoter_edges": n_trans_prom,
            "maintain_seeds": maintain_seeds,
        },
    }

    report = {
        "n_places": len(places),
        "n_blocks_axis": int(mean.shape[0]),
        "n_blocks": n_blocks,
        "n_runs": n_runs,
        "n_marked_seeds": n_marked_seeds,
        "n_interactions": n_interactions,
        "maintain_seeds": maintain_seeds,
        "token_presence_final": {
            "with_token_gt0": int(np.sum(active)),
            "zero_token": int(np.sum(~active)),
            "frac_with_token": frac_tok,
            "ever_active_any_block": int(np.sum(ever_active)),
            "always_zero": int(np.sum(always_zero)),
        },
        "by_layer_final_active": by_layer,
        "final_token_distribution_active": dist,
        "run_to_run_cv_final_active": run_cv,
        "temporal_cv_warm_active": temporal,
        "pileup_flags": {
            "p99_threshold": p99,
            "n_ge_p99": int(np.sum(hi)),
            "n_runaway_growth_gt5_and_hi": n_runaway_warn,
            "n_runaway_growth_gt50_and_hi": n_runaway_fail,
            "top_pileup": top_pileup,
        },
        "run_to_run_cv_stratified": run_cv_strat,
        "verdict_hints": hints,
        "verdict": verdict,
    }
    return SpnReasonablenessQC(
        stem=stem,
        ok=ok,
        warnings=warnings,
        failures=failures,
        report=report,
    )


def qc_spn_reasonableness_from_sim_dir(
    sim_dir: str | Path,
    stem: str,
    *,
    subgraph_dir: str | Path | None = None,
    marked_seed_names: Sequence[str] | None = None,
    warm_start_block: int = 10,
) -> SpnReasonablenessQC:
    """Load ``{stem}_trajectory.json`` + summary and run QC."""
    sim_dir = Path(sim_dir)
    traj_path = sim_dir / f"{stem}_trajectory.json"
    summary_path = sim_dir / f"{stem}_sim_summary.json"
    if not traj_path.is_file():
        raise FileNotFoundError(traj_path)
    traj = json.loads(traj_path.read_text(encoding="utf-8"))
    summary = (
        json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.is_file() else {}
    )
    places = traj["places"]
    mean = np.asarray(traj["mean_markings"], dtype=np.float64)
    std = np.asarray(traj["std_markings"], dtype=np.float64)

    edges = None
    if subgraph_dir is not None:
        sd = Path(subgraph_dir)
        for name in ("subgraph_edges.tsv", "gin_edges.tsv"):
            ep = sd / name
            if ep.is_file():
                edges = _load_edges_tsv(ep)
                break
        if marked_seed_names is None:
            # Prefer places recorded as fixed_high in sim summary when present;
            # else fall back to seed_manifest protein anchors (seed_protein mode).
            fh = summary.get("marked_seed_names") or summary.get("fixed_high")
            if fh:
                marked_seed_names = list(fh)
            else:
                man = sd / "seed_manifest.json"
                if man.is_file():
                    data = json.loads(man.read_text(encoding="utf-8"))
                    anchors = set(data.get("gin_anchor_nodes") or [])
                    place_set = set(places)
                    marked_seed_names = sorted(
                        n for n in anchors if n in place_set and is_markable_seed_place(n)
                    )

    return qc_spn_reasonableness_from_arrays(
        stem=stem,
        place_names=places,
        mean_markings=mean,
        std_markings=std,
        n_runs=int(traj.get("n_runs") or summary.get("runs") or 0),
        n_blocks=int(traj.get("n_blocks") or summary.get("blocks") or max(mean.shape[0] - 1, 0)),
        n_marked_seeds=int(summary.get("n_marked_seeds") or 0),
        n_interactions=int(summary.get("n_interactions") or 0),
        maintain_seeds=bool(summary.get("maintain_seeds", True)),
        marked_seed_names=marked_seed_names,
        mark_mode=str(summary.get("mark_mode") or "seed_protein"),
        warm_start_block=warm_start_block,
        subgraph_edges=edges,
    )


def emit_spn_reasonableness_qc(
    qc: SpnReasonablenessQC,
    *,
    out_json: Path | None = None,
    print_summary: bool = True,
) -> SpnReasonablenessQC:
    if out_json is not None:
        out_json = Path(out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(
            json.dumps(qc.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    if print_summary:
        status = "OK" if qc.ok else "FAIL"
        v = qc.report.get("verdict", {})
        print(
            f"[SPN-QC] {status} overall={v.get('overall')} "
            f"marked={qc.report.get('n_marked_seeds')} "
            f"dead={v.get('dead_places')} pileup={v.get('pileup')}",
            flush=True,
        )
        for w in qc.warnings[:8]:
            print(f"[SPN-QC][WARN] {w}", flush=True)
        for f in qc.failures[:8]:
            print(f"[SPN-QC][FAIL] {f}", flush=True)
        if out_json is not None:
            print(f"[SPN-QC] wrote {out_json}", flush=True)
    return qc

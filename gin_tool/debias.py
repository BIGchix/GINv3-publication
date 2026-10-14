"""Null-calibrated (debiased) scores and filtered analyze artifacts."""

from __future__ import annotations

import json
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from gin_tool.null_context import build_null_debiased_summary
from gin_tool.null_rwr import collect_null_rwr_ecdf_samples
from gin_tool.null_context import read_manifest, resolve_runs_path
from gin_tool.propagate import random_walk_restart

@dataclass
class RwrNullCalibration:
    """Pooled null node-score distribution for per-node calibration."""

    source: str  # precomputed | runtime | none
    k: int
    n_seeds_requested: int
    n_seeds_null_used: int
    rwr_alpha: float
    n_replicates: int
    pooled_scores: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64), repr=False)
    null_median: float = 0.0
    null_p90: float = 0.0
    null_max: float = 0.0
    notes: list[str] = field(default_factory=list)

    def has_pool(self) -> bool:
        return self.pooled_scores.size > 0


def _pool_from_precomputed(
    null_dir: Path,
    k: int,
    n_seeds: int,
    rwr_alpha: float,
) -> RwrNullCalibration:
    nearest_n, samples, _maxs = collect_null_rwr_ecdf_samples(
        null_dir, k, n_seeds, rwr_alpha=rwr_alpha
    )
    if samples:
        pooled = np.concatenate(samples)
        return RwrNullCalibration(
            source="precomputed",
            k=k,
            n_seeds_requested=n_seeds,
            n_seeds_null_used=nearest_n,
            rwr_alpha=rwr_alpha,
            n_replicates=len(samples),
            pooled_scores=pooled,
            null_median=float(np.median(pooled)),
            null_p90=float(np.quantile(pooled, 0.9)),
            null_max=float(np.max(pooled)),
        )
    return RwrNullCalibration(
        source="none",
        k=k,
        n_seeds_requested=n_seeds,
        n_seeds_null_used=nearest_n,
        rwr_alpha=rwr_alpha,
        n_replicates=0,
        notes=[
            f"No rwr_ecdf_sample in precompute for k={k}, n_seeds={nearest_n}."
        ],
    )


def _runtime_null_pool(
    sifi_path: str,
    k: int,
    n_seeds: int,
    *,
    rwr_alpha: float = 0.3,
    n_replicates: int = 25,
    max_sub_nodes: int = 25_000,
    max_lines: int | None = None,
    base_seed: int = 42,
    seed_layer: str | None = None,
    seed_layer_counts: dict[str, int] | None = None,
    node_list_dir: Path | None = None,
    k_by_name: dict[str, int] | None = None,
) -> RwrNullCalibration:
    from gin_tool.io import (
        build_indexed_edges,
        directed_csr_from_indexed,
        symmetrize_undirected_unweighted,
    )
    from gin_tool.mapping import real_entity_indices
    from gin_tool.node_lists import (
        RUNTIME_ONLY_LAYERS,
        UNKNOWN_SEED_LAYER,
        build_pools_and_layer_map,
        default_k_for_layer,
        max_k_for_layer,
        validate_k_for_seed_layer,
    )
    from gin_tool.runtime_null import (
        run_runtime_null_replicate,
        runtime_null_pool_stratified,
    )

    notes: list[str] = []
    k_eff = int(max(k_by_name.values())) if k_by_name else int(k)

    ed = build_indexed_edges(sifi_path, max_lines=max_lines)
    sym = symmetrize_undirected_unweighted(directed_csr_from_indexed(ed))
    global_pool = real_entity_indices(ed.index_to_name)

    buckets: dict[str, list[int]] = {}
    if node_list_dir is not None and node_list_dir.is_dir():
        buckets, _, _, _, _ = build_pools_and_layer_map(
            node_list_dir, ed.index_to_name, ed.name_to_index
        )

    active_layers = {
        ly: c for ly, c in (seed_layer_counts or {}).items() if int(c) > 0
    }
    use_stratified = len(active_layers) > 1

    if use_stratified:
        chunks, strat_notes, used, skipped = runtime_null_pool_stratified(
            sym,
            ed.index_to_name,
            buckets,
            active_layers,
            n_seeds,
            n_replicates=n_replicates,
            rwr_alpha=rwr_alpha,
            max_sub_nodes=max_sub_nodes,
            base_seed=base_seed,
            global_pool=global_pool,
        )
        notes.extend(strat_notes)
        if skipped:
            notes.append(
                f"runtime null: {skipped} stratified replicates skipped "
                f"(empty/large subgraph)."
            )
    else:
        if active_layers:
            seed_layer = max(active_layers.items(), key=lambda x: (x[1], x[0]))[0]
            if seed_layer == UNKNOWN_SEED_LAYER:
                seed_layer = None

        if seed_layer is not None:
            err = validate_k_for_seed_layer(seed_layer, k_eff)
            if err:
                notes.append(err)
                return RwrNullCalibration(
                    source="none",
                    k=k_eff,
                    n_seeds_requested=n_seeds,
                    n_seeds_null_used=n_seeds,
                    rwr_alpha=rwr_alpha,
                    n_replicates=0,
                    notes=notes,
                )
            cap = max_k_for_layer(seed_layer)
            if k_eff > cap:
                k_eff = cap
                notes.append(
                    f"runtime null: clamped k to {k_eff} for layer {seed_layer!r}."
                )

        pool: list[int]
        if seed_layer is not None and buckets:
            pool = list(buckets.get(seed_layer, []))
            k_rep = default_k_for_layer(seed_layer)
            if not pool:
                notes.append(f"runtime null: empty nodeList pool for layer {seed_layer!r}.")
            else:
                notes.append(
                    f"runtime null: sampling from nodeList layer {seed_layer!r} "
                    f"(pool={len(pool)}, k={k_rep})."
                )
        else:
            if seed_layer in RUNTIME_ONLY_LAYERS:
                notes.append(
                    f"runtime null: seed_layer={seed_layer!r} but no node_list_dir; "
                    "using all real entities (not recommended)."
                )
            pool = global_pool
            k_rep = k_eff

        if not pool:
            notes.append("No real (non-virtual) nodes for runtime null.")
            return RwrNullCalibration(
                source="none",
                k=k_eff,
                n_seeds_requested=n_seeds,
                n_seeds_null_used=n_seeds,
                rwr_alpha=rwr_alpha,
                n_replicates=0,
                notes=notes,
            )

        chunks = []
        used = 0
        skipped = 0
        for rep in range(n_replicates):
            rng = random.Random(base_seed + rep * 1_000_003 + k_rep * 9_999 + n_seeds)
            rec = run_runtime_null_replicate(
                sym,
                ed.index_to_name,
                pool,
                k_rep,
                n_seeds,
                rng,
                rwr_alpha=rwr_alpha,
                max_sub_nodes=max_sub_nodes,
            )
            if not rec["ok"]:
                skipped += 1
                continue
            chunks.append(rec["rwr_scores_sample"])
            used += 1
        if skipped:
            notes.append(
                f"runtime null: {skipped}/{n_replicates} replicates skipped "
                f"(empty/large subgraph or no seeds)."
            )

    if not chunks:
        notes.append("runtime null produced no score pool.")
        return RwrNullCalibration(
            source="none",
            k=k_eff,
            n_seeds_requested=n_seeds,
            n_seeds_null_used=n_seeds,
            rwr_alpha=rwr_alpha,
            n_replicates=0,
            notes=notes,
        )

    pooled = np.concatenate(chunks)
    return RwrNullCalibration(
        source="runtime",
        k=k_eff,
        n_seeds_requested=n_seeds,
        n_seeds_null_used=n_seeds,
        rwr_alpha=rwr_alpha,
        n_replicates=used,
        pooled_scores=pooled,
        null_median=float(np.median(pooled)),
        null_p90=float(np.quantile(pooled, 0.9)),
        null_max=float(np.max(pooled)),
        notes=notes,
    )


def build_rwr_null_calibration(
    null_dir: Path | None,
    k: int,
    n_seeds: int,
    *,
    rwr_alpha: float = 0.3,
    sifi_path: str | None = None,
    max_lines: int | None = None,
    runtime_replicates: int = 25,
    runtime_max_sub_nodes: int = 25_000,
    obs_subgraph_nodes: int | None = None,
    seed_layer: str | None = None,
    seed_layer_counts: dict[str, int] | None = None,
    node_list_dir: Path | None = None,
    k_by_name: dict[str, int] | None = None,
) -> RwrNullCalibration:
    cal = RwrNullCalibration(
        source="none",
        k=k,
        n_seeds_requested=n_seeds,
        n_seeds_null_used=n_seeds,
        rwr_alpha=rwr_alpha,
        n_replicates=0,
        notes=[],
    )
    if null_dir is not None and null_dir.is_dir():
        cal = _pool_from_precomputed(null_dir, k, n_seeds, rwr_alpha)

    if cal.has_pool():
        return cal

    if sifi_path is None:
        cal.notes.append("No sifi_path for runtime null fallback.")
        return cal

    if obs_subgraph_nodes is not None and obs_subgraph_nodes > runtime_max_sub_nodes:
        cal.notes.append(
            f"Observed subgraph n={obs_subgraph_nodes} > {runtime_max_sub_nodes}; "
            "skipped runtime RWR null."
        )
        return cal

    runtime = _runtime_null_pool(
        sifi_path,
        k,
        n_seeds,
        rwr_alpha=rwr_alpha,
        n_replicates=runtime_replicates,
        max_sub_nodes=runtime_max_sub_nodes,
        max_lines=max_lines,
        seed_layer=seed_layer,
        seed_layer_counts=seed_layer_counts,
        node_list_dir=node_list_dir,
        k_by_name=k_by_name,
    )
    runtime.notes = cal.notes + runtime.notes
    return runtime


def empirical_upper_tail_p(values: np.ndarray, null_pool: np.ndarray) -> np.ndarray:
    if null_pool.size == 0:
        return np.full(values.shape, np.nan)
    sorted_null = np.sort(null_pool)
    n = len(sorted_null)
    idx = np.searchsorted(sorted_null, values, side="left")
    return (n - idx) / float(n)


def calibrate_rwr_scores(
    scores: np.ndarray,
    cal: RwrNullCalibration,
    *,
    debias_mode: str = "residual_median",
    tail_alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tail_p = empirical_upper_tail_p(scores, cal.pooled_scores)
    if not cal.has_pool():
        return scores.copy(), tail_p, np.zeros(scores.shape, dtype=bool)

    eps = 1e-12
    if debias_mode == "residual_p90":
        debiased = np.maximum(0.0, scores - cal.null_p90)
    elif debias_mode == "divide_median":
        debiased = scores / (cal.null_median + eps)
    else:
        debiased = np.maximum(0.0, scores - cal.null_median)

    sig = np.isfinite(tail_p) & (tail_p <= tail_alpha) & (debiased > 0)
    return debiased, tail_p, sig


def write_propagation_scores_debiased_tsv(
    path: Path,
    names: Sequence[str],
    scores: np.ndarray,
    debiased: np.ndarray,
    tail_p: np.ndarray,
    significant: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("node\tscore\tscore_debiased\ttail_p_null_ge\tnull_significant\n")
        for nm, sc, db, tp, sig in zip(names, scores, debiased, tail_p, significant):
            tp_s = f"{tp:g}" if np.isfinite(tp) else "NA"
            f.write(f"{nm}\t{sc:g}\t{db:g}\t{tp_s}\t{int(bool(sig))}\n")


def _iter_null_rows_all_n(runs_path: Path, k: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    from gin_tool.null_context import _open_runs_text

    for line in _open_runs_text(runs_path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("error"):
            continue
        if int(r["k"]) == int(k):
            rows.append(r)
    return rows


def size_matched_null_pool(
    null_dir: Path,
    k: int,
    obs_nodes: int,
    *,
    band_frac: float = 0.35,
    min_rows: int = 5,
) -> tuple[np.ndarray, list[str]]:
    """
    Pool RWR scores from null replicates with similar subgraph size (any n_seeds).

    Useful when grid n_seeds (e.g. 50) yields huge subgraphs without per-node RWR
    in precompute, but smaller random draws at other grid points have rwr samples.
    """
    manifest = read_manifest(null_dir)
    runs_path = resolve_runs_path(manifest, null_dir)
    lo = max(1, obs_nodes * (1.0 - band_frac))
    hi = obs_nodes * (1.0 + band_frac)
    chunks: list[np.ndarray] = []
    notes: list[str] = []
    for r in _iter_null_rows_all_n(runs_path, k):
        n_sub = r.get("subgraph_nodes")
        if n_sub is None:
            continue
        if not (lo <= float(n_sub) <= hi):
            continue
        samp = r.get("rwr_ecdf_sample")
        if samp:
            chunks.append(np.asarray(samp, dtype=np.float64))
    if len(chunks) < min_rows:
        notes.append(
            f"size_matched pool: only {len(chunks)} replicates in "
            f"[{lo:.0f}, {hi:.0f}] nodes with rwr_ecdf_sample."
        )
        return np.array([], dtype=np.float64), notes
    notes.append(
        f"size_matched pool: {len(chunks)} replicates, "
        f"node band [{lo:.0f}, {hi:.0f}] around obs={obs_nodes}."
    )
    return np.concatenate(chunks), notes


def build_rwr_null_calibration_extended(
    null_dir: Path | None,
    k: int,
    n_seeds: int,
    obs_subgraph_nodes: int,
    *,
    rwr_alpha: float = 0.3,
    sifi_path: str | None = None,
    max_lines: int | None = None,
    runtime_replicates: int = 25,
    runtime_max_sub_nodes: int = 25_000,
    seed_layer: str | None = None,
    seed_layer_counts: dict[str, int] | None = None,
    node_list_dir: Path | None = None,
    k_by_name: dict[str, int] | None = None,
    runtime_only: bool = False,
) -> RwrNullCalibration:
    """Try grid-matched precompute, then size-matched pool, then runtime null."""
    from gin_tool.node_lists import RUNTIME_ONLY_LAYERS

    precompute_dir = None if runtime_only else null_dir
    if seed_layer in RUNTIME_ONLY_LAYERS:
        precompute_dir = None
    if seed_layer_counts and len({ly for ly, c in seed_layer_counts.items() if c > 0}) > 1:
        precompute_dir = None

    cal = build_rwr_null_calibration(
        precompute_dir,
        k,
        n_seeds,
        rwr_alpha=rwr_alpha,
        sifi_path=None,
        obs_subgraph_nodes=obs_subgraph_nodes,
    )
    if cal.has_pool():
        return cal

    if precompute_dir is not None and precompute_dir.is_dir():
        pooled, notes = size_matched_null_pool(precompute_dir, k, obs_subgraph_nodes)
        if pooled.size > 0:
            return RwrNullCalibration(
                source="precomputed_size_matched",
                k=k,
                n_seeds_requested=n_seeds,
                n_seeds_null_used=cal.n_seeds_null_used,
                rwr_alpha=rwr_alpha,
                n_replicates=-1,
                pooled_scores=pooled,
                null_median=float(np.median(pooled)),
                null_p90=float(np.quantile(pooled, 0.9)),
                null_max=float(np.max(pooled)),
                notes=cal.notes + notes,
            )
        cal.notes.extend(notes)

    return build_rwr_null_calibration(
        precompute_dir,
        k,
        n_seeds,
        rwr_alpha=rwr_alpha,
        sifi_path=sifi_path,
        max_lines=max_lines,
        runtime_replicates=runtime_replicates,
        runtime_max_sub_nodes=runtime_max_sub_nodes,
        obs_subgraph_nodes=obs_subgraph_nodes,
        seed_layer=seed_layer,
        seed_layer_counts=seed_layer_counts,
        node_list_dir=node_list_dir,
        k_by_name=k_by_name,
    )


def filter_flows_by_null(
    flows_path: Path,
    out_path: Path,
    *,
    null_summary: dict[str, Any],
    tail_alpha: float = 0.05,
) -> dict[str, Any]:
    """Keep flows only if path enumeration is not explainable by null (when stats exist)."""
    flows = json.loads(flows_path.read_text(encoding="utf-8"))
    block = null_summary.get("flows", {})
    vs = block.get("vs_null_paths")
    meta: dict[str, Any] = {"kept": len(flows), "total": len(flows), "filtered": False}

    if vs and vs.get("raw_path_count"):
        tp = vs["raw_path_count"].get("tail_p_larger_than_null_median")
        if tp is None:
            tp = null_summary.get("subgraph_size_vs_null", {}).get(
                "tail_p_edges_as_large_as_observed"
            )
        if tp is not None and float(tp) > tail_alpha:
            meta["filtered"] = True
            meta["note"] = (
                "Path count is typical under random-seed null; "
                "keeping flows but flagging low structural surprise."
            )
            meta["kept"] = len(flows)
        else:
            meta["note"] = "Path count exceeds random-seed null expectation."
    else:
        meta["note"] = "No path null stats; flows copied unchanged."

    out_path.write_text(json.dumps(flows, indent=2), encoding="utf-8")
    return meta


def run_debias_pipeline(
    out_dir: Path,
    *,
    sifi_path: str,
    k: int,
    n_seeds: int,
    resolved_seed_names: list[str],
    scores: np.ndarray,
    sub_names: list[str],
    sub_sym_n: int,
    sub_sym_nnz: int,
    null_precompute_dir: Path | None = None,
    rwr_alpha: float = 0.3,
    max_lines: int | None = None,
    debias_mode: str = "residual_median",
    tail_alpha: float = 0.05,
    runtime_replicates: int = 25,
    node_list_dir: Path | None = None,
    seed_layer: str | None = None,
    seed_layer_counts: dict[str, int] | None = None,
    k_by_name: dict[str, int] | None = None,
    runtime_only: bool = False,
    analyze: bool = False,
    obs_louvain_n: int | None = None,
    obs_louvain_largest_frac: float | None = None,
    obs_louvain_modularity: float | None = None,
    obs_raw_paths: int | None = None,
    obs_ranked_paths: int | None = None,
    obs_best_path_score: tuple[int, int] | None = None,
) -> dict[str, Any]:
    """
    Write debiased RWR table, null context JSON, and optional filtered analyze copies.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    cal = build_rwr_null_calibration_extended(
        null_precompute_dir,
        k,
        n_seeds,
        sub_sym_n,
        rwr_alpha=rwr_alpha,
        sifi_path=sifi_path,
        max_lines=max_lines,
        runtime_replicates=runtime_replicates,
        seed_layer=seed_layer,
        seed_layer_counts=seed_layer_counts,
        node_list_dir=node_list_dir,
        k_by_name=k_by_name,
        runtime_only=runtime_only,
    )
    debiased, tail_p, sig = calibrate_rwr_scores(
        scores, cal, debias_mode=debias_mode, tail_alpha=tail_alpha
    )

    write_propagation_scores_debiased_tsv(
        out_dir / "propagation_scores_debiased.tsv",
        sub_names,
        scores,
        debiased,
        tail_p,
        sig,
    )

    obs_rwr_max = float(np.max(scores)) if scores.size else None
    obs_rwr_p90 = float(np.quantile(scores, 0.9)) if scores.size else None

    if (
        not runtime_only
        and null_precompute_dir is not None
        and null_precompute_dir.is_dir()
    ):
        summary = build_null_debiased_summary(
            null_precompute_dir,
            k,
            n_seeds,
            sub_sym_n,
            sub_sym_nnz,
            sifi_path=sifi_path,
            obs_louvain_n=obs_louvain_n,
            obs_louvain_largest_frac=obs_louvain_largest_frac,
            obs_louvain_modularity=obs_louvain_modularity,
            obs_raw_paths=obs_raw_paths,
            obs_ranked_paths=obs_ranked_paths,
            obs_best_path_score=obs_best_path_score,
            obs_rwr_max=obs_rwr_max,
            obs_rwr_p90=obs_rwr_p90,
        )
    else:
        summary = {
            "schema_version": "2.0",
            "null_policy": "runtime_only",
            "k": k,
            "n_seeds": n_seeds,
            "seed_layer": seed_layer,
            "seed_layer_counts": seed_layer_counts,
            "k_by_seed": k_by_name,
            "observed": {
                "subgraph_nodes": sub_sym_n,
                "subgraph_edges_sym_nnz": sub_sym_nnz,
                "rwr_max": obs_rwr_max,
                "rwr_p90": obs_rwr_p90,
            },
            "note": "RWR debias uses runtime null from nodeList layer pool; no precompute comparison.",
        }
    summary["schema_version"] = "2.0"
    summary["debiasing"] = {
        "rwr_calibration_source": cal.source,
        "rwr_null_replicates": cal.n_replicates,
        "debias_mode": debias_mode,
        "tail_alpha": tail_alpha,
        "null_median_score": cal.null_median,
        "null_p90_score": cal.null_p90,
        "n_nodes_significant": int(np.sum(sig)),
        "n_nodes_tested": int(scores.size),
        "notes": cal.notes,
        "outputs": {
            "scores_debiased": "propagation_scores_debiased.tsv",
            "ranking_recommendation": (
                "Use score_debiased for ranking; null_significant=1 for strict hits."
            ),
        },
    }

    if analyze:
        flows_src = out_dir / "flows.json"
        if flows_src.is_file():
            fm = filter_flows_by_null(
                flows_src,
                out_dir / "flows_debiased.json",
                null_summary=summary,
                tail_alpha=tail_alpha,
            )
            summary["debiasing"]["flows_filter"] = fm

        mod_src = out_dir / "modules.json"
        if mod_src.is_file() and summary.get("modules", {}).get("vs_null_louvain"):
            summary["debiasing"]["modules_note"] = (
                "modules.json unchanged; interpret with modules.vs_null_louvain in summary."
            )

    (out_dir / "null_debiased_context.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary

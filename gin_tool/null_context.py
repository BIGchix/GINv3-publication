"""Compare observed k-hop subgraph stats to a precomputed random-seed null (JSON Lines)."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator

import numpy as np


def sha256_file(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_manifest(null_dir: Path) -> dict[str, Any]:
    p = null_dir / "manifest.json"
    if not p.is_file():
        raise FileNotFoundError(f"missing manifest: {p}")
    return json.loads(p.read_text(encoding="utf-8"))


def resolve_runs_path(manifest: dict[str, Any], null_dir: Path) -> Path:
    rf = manifest.get("runs_file")
    if not rf:
        raise ValueError("manifest missing runs_file")
    p = Path(rf)
    if p.is_file():
        return p
    name = p.name
    cand = null_dir / name
    if cand.is_file():
        return cand
    raise FileNotFoundError(f"runs file not found: {rf} or {cand}")


def nearest_grid_n(n_user: int, n_min: int, n_max: int, n_step: int) -> tuple[int, str]:
    if n_user <= 0:
        raise ValueError("n_seeds must be positive")
    if n_min <= n_user <= n_max and (n_user - n_min) % n_step == 0:
        return n_user, "exact_grid_match"
    grid: list[int] = []
    x = n_min
    while x <= n_max:
        grid.append(x)
        x += n_step
    if not grid:
        raise ValueError("empty n_seeds grid in manifest")
    nearest = min(grid, key=lambda g: abs(g - n_user))
    return nearest, f"nearest_grid_point(requested={n_user}, using_n_seeds={nearest})"


def _open_runs_text(path: Path) -> Iterator[str]:
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as f:
            yield from f
    else:
        with path.open("r", encoding="utf-8") as f:
            yield from f


def collect_null_metrics(runs_path: Path, k: int, n_seeds: int) -> dict[str, np.ndarray]:
    nodes: list[int] = []
    nnzs: list[int] = []
    lou_n: list[float] = []
    lou_frac: list[float] = []
    lou_q: list[float] = []
    raw_paths: list[float] = []
    ranked_paths: list[float] = []
    best_layers: list[float] = []
    rwr_max: list[float] = []
    rwr_p90: list[float] = []
    for line in _open_runs_text(runs_path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("error"):
            continue
        if int(r["k"]) != int(k) or int(r["n_seeds"]) != int(n_seeds):
            continue
        nodes.append(int(r["subgraph_nodes"]))
        nnzs.append(int(r["subgraph_edges_sym_nnz"]))
        if r.get("louvain_num_communities") is not None:
            lou_n.append(float(r["louvain_num_communities"]))
        if r.get("louvain_largest_comm_frac") is not None:
            lou_frac.append(float(r["louvain_largest_comm_frac"]))
        if r.get("louvain_modularity_estimate") is not None:
            lou_q.append(float(r["louvain_modularity_estimate"]))
        if r.get("directed_raw_path_count") is not None:
            raw_paths.append(float(r["directed_raw_path_count"]))
        if r.get("directed_ranked_path_count") is not None:
            ranked_paths.append(float(r["directed_ranked_path_count"]))
        if r.get("best_cross_layer_layers") is not None:
            best_layers.append(float(r["best_cross_layer_layers"]))
        if r.get("rwr_score_max") is not None:
            rwr_max.append(float(r["rwr_score_max"]))
        if r.get("rwr_score_p90") is not None:
            rwr_p90.append(float(r["rwr_score_p90"]))
    if not nodes:
        raise ValueError(f"No null rows for k={k}, n_seeds={n_seeds} in {runs_path}")
    na = np.array(nodes, dtype=np.float64)
    nz = np.array(nnzs, dtype=np.float64)
    md = nz / np.maximum(na, 1.0)
    out: dict[str, np.ndarray] = {
        "subgraph_nodes": na,
        "subgraph_edges_sym_nnz": nz,
        "sym_mean_deg_proxy": md,
    }
    if lou_n:
        out["louvain_num_communities"] = np.array(lou_n, dtype=np.float64)
    if lou_frac:
        out["louvain_largest_comm_frac"] = np.array(lou_frac, dtype=np.float64)
    if lou_q:
        out["louvain_modularity_estimate"] = np.array(lou_q, dtype=np.float64)
    if raw_paths:
        out["directed_raw_path_count"] = np.array(raw_paths, dtype=np.float64)
    if ranked_paths:
        out["directed_ranked_path_count"] = np.array(ranked_paths, dtype=np.float64)
    if best_layers:
        out["best_cross_layer_layers"] = np.array(best_layers, dtype=np.float64)
    if rwr_max:
        out["rwr_score_max"] = np.array(rwr_max, dtype=np.float64)
    if rwr_p90:
        out["rwr_score_p90"] = np.array(rwr_p90, dtype=np.float64)
    return out


def _quantile_summary(value: float, arr: np.ndarray) -> dict[str, Any]:
    a = np.sort(arr)
    q = float(np.mean(a < float(value)))
    return {
        "observed": float(value),
        "null_median": float(np.median(a)),
        "null_p05": float(np.quantile(a, 0.05)),
        "null_p95": float(np.quantile(a, 0.95)),
        "empirical_quantile_lt": q,
        "tail_p_larger_than_null_median": float(np.mean(a >= float(value))),
    }


def _tail_p(value: float, arr: np.ndarray, *, upper: bool) -> float:
    """One-sided: upper=True → P(null >= value); upper=False → P(null <= value)."""
    if upper:
        return float(np.mean(arr >= float(value)))
    return float(np.mean(arr <= float(value)))


def build_null_debiased_summary(
    null_dir: Path,
    k: int,
    n_seeds_user: int,
    obs_nodes: int,
    obs_nnz: int,
    *,
    sifi_path: str | None = None,
    obs_louvain_n: int | None = None,
    obs_louvain_largest_frac: float | None = None,
    obs_louvain_modularity: float | None = None,
    obs_raw_paths: int | None = None,
    obs_ranked_paths: int | None = None,
    obs_best_path_score: tuple[int, int] | None = None,
    obs_rwr_max: float | None = None,
    obs_rwr_p90: float | None = None,
) -> dict[str, Any]:
    """
    Compare observed sym subgraph (nodes, nnz) to precomputed null at same k and
    nearest grid ``n_seeds``. Optional Louvain / path stats are contextual when
    present in the null file or as observed-only.
    """
    manifest = read_manifest(null_dir)
    kv = manifest.get("k_values")
    if isinstance(kv, list) and kv and int(k) not in [int(x) for x in kv]:
        raise ValueError(f"k={k} not in manifest k_values={kv}")
    runs_path = resolve_runs_path(manifest, null_dir)
    n_min = int(manifest["n_min"])
    n_max = int(manifest["n_max"])
    n_step = int(manifest["n_step"])
    nearest_n, match_note = nearest_grid_n(n_seeds_user, n_min, n_max, n_step)
    arrays = collect_null_metrics(runs_path, k, nearest_n)
    n_null = int(arrays["subgraph_nodes"].size)

    obs_deg = float(obs_nnz) / max(float(obs_nodes), 1.0)
    null_deg = arrays["sym_mean_deg_proxy"]

    match_block: dict[str, Any] = {
        "k": k,
        "n_seeds_requested": n_seeds_user,
        "n_seeds_null_used": nearest_n,
        "match_note": match_note,
        "n_null_replicates": n_null,
    }
    if abs(nearest_n - n_seeds_user) > n_step:
        match_block["warning"] = (
            f"requested n_seeds={n_seeds_user} is far from matched grid value {nearest_n} "
            f"(n_step={n_step}); null compares random draws of size {nearest_n}, not {n_seeds_user}."
        )

    out: dict[str, Any] = {
        "schema_version": "1.0",
        "null_dir": str(null_dir.resolve()),
        "manifest_k_values": manifest.get("k_values"),
        "manifest_gin_sha256": manifest.get("gin_sha256"),
        "match": match_block,
        "subgraph_size_vs_null": {
            "subgraph_nodes": _quantile_summary(float(obs_nodes), arrays["subgraph_nodes"]),
            "subgraph_edges_sym_nnz": _quantile_summary(float(obs_nnz), arrays["subgraph_edges_sym_nnz"]),
            "sym_mean_deg_proxy": _quantile_summary(obs_deg, null_deg),
            "tail_p_nodes_as_large_as_observed": _tail_p(float(obs_nodes), arrays["subgraph_nodes"], upper=True),
            "tail_p_edges_as_large_as_observed": _tail_p(float(obs_nnz), arrays["subgraph_edges_sym_nnz"], upper=True),
        },
        "interpretation": {
            "en": (
                "Random-seed null at the same k and (matched) n_seeds counts how large a "
                "k-hop sym subgraph typically is in GIN. Low tail_p for nodes/edges means "
                "your hit set induces a **smaller / sparser** neighborhood than most random "
                "gene-sized seed draws (more **localized**). High tail_p means **broader** "
                "spread than random seeds of the same count."
            ),
            "zh": (
                "在相同 k 与匹配的 n_seeds 下，null 描述「随机种子」在 GIN 上 k-hop 无向支撑子图的"
                "典型规模。节点/边的 upper tail p 很小表示子图比多数随机种子更小/更稀（更局部）；"
                "很大表示比随机种子更发散。"
            ),
        },
    }

    if sifi_path and manifest.get("gin_sha256") not in (None, "skipped"):
        p = Path(sifi_path)
        if p.is_file():
            h = sha256_file(p)
            out["sifi_sha256_check"] = {
                "current": h,
                "manifest": manifest["gin_sha256"],
                "match": h == manifest["gin_sha256"],
            }

    mod_block: dict[str, Any] = {
        "observed": {
            "num_communities": obs_louvain_n,
            "largest_community_frac": obs_louvain_largest_frac,
            "modularity": obs_louvain_modularity,
        },
    }
    if (
        obs_louvain_n is not None
        and "louvain_num_communities" in arrays
        and arrays["louvain_num_communities"].size > 0
    ):
        a = arrays["louvain_num_communities"]
        mod_block["vs_null_louvain"] = {
            "num_communities": _quantile_summary(float(obs_louvain_n), a),
            "tail_p_as_many_communities": _tail_p(float(obs_louvain_n), a, upper=True),
        }
        if obs_louvain_modularity is not None and "louvain_modularity_estimate" in arrays:
            mq = arrays["louvain_modularity_estimate"]
            mod_block["vs_null_louvain"]["modularity"] = _quantile_summary(
                float(obs_louvain_modularity), mq
            )
    else:
        mod_block["vs_null_louvain"] = None
        mod_block["note"] = (
            "Precomputed null has no usable Louvain columns (re-run with default structure stats "
            "and ensure most replicates have n_sub ≤ louvain_max_n, or increase louvain_max_n if RAM allows)."
        )
    out["modules"] = mod_block

    flow_block: dict[str, Any] = {
        "observed": {
            "raw_path_count": obs_raw_paths,
            "ranked_path_count": obs_ranked_paths,
            "best_cross_layer_score": list(obs_best_path_score) if obs_best_path_score else None,
        },
        "vs_null_paths": None,
        "note": None,
    }
    if (
        obs_raw_paths is not None
        and "directed_raw_path_count" in arrays
        and arrays["directed_raw_path_count"].size > 0
    ):
        pr = arrays["directed_raw_path_count"]
        flow_block["vs_null_paths"] = {
            "raw_path_count": _quantile_summary(float(obs_raw_paths), pr),
            "tail_p_as_many_raw_paths": _tail_p(float(obs_raw_paths), pr, upper=True),
        }
        if obs_ranked_paths is not None and "directed_ranked_path_count" in arrays:
            rk = arrays["directed_ranked_path_count"]
            flow_block["vs_null_paths"]["ranked_path_count"] = _quantile_summary(
                float(obs_ranked_paths), rk
            )
        if obs_best_path_score is not None and "best_cross_layer_layers" in arrays:
            bl = arrays["best_cross_layer_layers"]
            flow_block["vs_null_paths"]["best_cross_layer_layers"] = _quantile_summary(
                float(obs_best_path_score[0]), bl
            )
    else:
        flow_block["note"] = (
            "Null has no path-count columns for this (k, n_seeds) slice (paths skipped in null "
            "when n_sub > path_enum_max_n, or --no-structure-stats was used)."
        )
    out["flows"] = flow_block

    rwr_block: dict[str, Any] = {
        "observed": {"score_max": obs_rwr_max, "score_p90": obs_rwr_p90},
        "vs_null_rwr": None,
        "note": None,
    }
    if obs_rwr_max is not None and "rwr_score_max" in arrays and arrays["rwr_score_max"].size > 0:
        rm = arrays["rwr_score_max"]
        rwr_block["vs_null_rwr"] = {
            "score_max": _quantile_summary(float(obs_rwr_max), rm),
            "tail_p_max_as_large": _tail_p(float(obs_rwr_max), rm, upper=True),
        }
        if obs_rwr_p90 is not None and "rwr_score_p90" in arrays:
            rp = arrays["rwr_score_p90"]
            rwr_block["vs_null_rwr"]["score_p90"] = _quantile_summary(
                float(obs_rwr_p90), rp
            )
    else:
        rwr_block["note"] = (
            "Null file has no rwr_score_max rows (re-run precompute_null with schema 1.2, "
            "or ensure n_sub <= rwr_max_n for most replicates)."
        )
    out["rwr_scores"] = rwr_block

    return out

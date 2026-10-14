"""Load RWR null samples from precompute_null JSON Lines (schema 1.2+)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from gin_tool.null_context import nearest_grid_n, read_manifest, resolve_runs_path


def _open_runs_text(path: Path) -> Iterator[str]:
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as f:
            yield from f
    else:
        with path.open("r", encoding="utf-8") as f:
            yield from f


def collect_null_rwr_rows(
    runs_path: Path,
    k: int,
    n_seeds: int,
    *,
    rwr_alpha: float | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in _open_runs_text(runs_path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("error"):
            continue
        if int(r["k"]) != int(k) or int(r["n_seeds"]) != int(n_seeds):
            continue
        if rwr_alpha is not None and r.get("rwr_alpha") is not None:
            if abs(float(r["rwr_alpha"]) - float(rwr_alpha)) > 1e-9:
                continue
        rows.append(r)
    return rows


def collect_null_rwr_ecdf_samples(
    null_dir: Path,
    k: int,
    n_seeds_user: int,
    *,
    rwr_alpha: float | None = None,
) -> tuple[int, list[np.ndarray], list[float]]:
    """
    Return (n_seeds_used, list of per-replicate ECDF sample arrays, null max scores).

    Uses nearest grid ``n_seeds`` from manifest when user count is off-grid.
    """
    manifest = read_manifest(null_dir)
    runs_path = resolve_runs_path(manifest, null_dir)
    n_min = int(manifest["n_min"])
    n_max = int(manifest["n_max"])
    n_step = int(manifest["n_step"])
    nearest_n, _ = nearest_grid_n(n_seeds_user, n_min, n_max, n_step)
    alpha = rwr_alpha
    if alpha is None and manifest.get("rwr_alpha") is not None:
        alpha = float(manifest["rwr_alpha"])

    rows = collect_null_rwr_rows(runs_path, k, nearest_n, rwr_alpha=alpha)
    samples: list[np.ndarray] = []
    max_scores: list[float] = []
    for r in rows:
        samp = r.get("rwr_ecdf_sample")
        if samp:
            samples.append(np.asarray(samp, dtype=np.float64))
        if r.get("rwr_score_max") is not None:
            max_scores.append(float(r["rwr_score_max"]))
    return nearest_n, samples, max_scores


def null_ecdf_median_curve(samples: list[np.ndarray], grid: np.ndarray) -> np.ndarray:
    """Median ECDF across null replicate subsamples at fixed score grid."""
    if not samples:
        return np.full_like(grid, np.nan)
    curves = []
    for s in samples:
        s = np.sort(s)
        y = np.searchsorted(s, grid, side="right") / max(len(s), 1)
        curves.append(y)
    return np.median(np.stack(curves, axis=0), axis=0)

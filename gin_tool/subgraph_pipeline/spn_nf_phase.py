"""TP53–MDM2 (or generic pair) lagged phase / NF pass-rate from MC trajectories."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.stats import pearsonr, t as student_t


def acf1(z: np.ndarray) -> float:
    z = np.asarray(z, float)
    z = z - z.mean()
    if z.size < 3 or float(np.std(z)) == 0:
        return float("nan")
    return float(np.corrcoef(z[:-1], z[1:])[0, 1])


def n_eff_ar1(n: int, a: float, b: float) -> int:
    if not (math.isfinite(a) and math.isfinite(b)):
        return max(3, n)
    denom = 1.0 + a * b
    if abs(denom) < 1e-9:
        return 3
    return max(3, int(round(n * (1.0 - a * b) / denom)))


def p_from_r(r: float, n: int) -> float:
    if n <= 2 or not math.isfinite(r) or abs(r) >= 1.0:
        return float("nan")
    tstat = r * math.sqrt((n - 2) / (1.0 - r * r + 1e-12))
    return float(2.0 * student_t.sf(abs(tstat), n - 2))


def lagged_corr(x: np.ndarray, y: np.ndarray, lag: int) -> tuple[float, int]:
    if lag >= 0:
        a, b = x[: len(x) - lag or None], y[lag:]
    else:
        a, b = x[-lag:], y[: len(y) + lag or None]
    n = min(len(a), len(b))
    if n < 10 or float(np.std(a)) == 0 or float(np.std(b)) == 0:
        return float("nan"), n
    r, _ = pearsonr(a[:n], b[:n])
    return float(r), n


def phase_from_series(
    x: np.ndarray,
    y: np.ndarray,
    *,
    warm: int = 20,
    max_lag: int = 40,
) -> dict[str, Any]:
    """Lagged Pearson phase for one run (x=leader e.g. TP53, y=follower e.g. MDM2)."""
    x = np.asarray(x, float)[warm:]
    y = np.asarray(y, float)[warm:]
    if x.size < max_lag + 10:
        return {"error": "series too short", "induction_strict": False, "lag_gt0_pos": False}

    lags = np.arange(-max_lag, max_lag + 1)
    rs = np.array([lagged_corr(x, y, int(lag))[0] for lag in lags], float)
    if not np.any(np.isfinite(rs)):
        return {"error": "degenerate series", "induction_strict": False, "lag_gt0_pos": False}
    best_i = int(np.nanargmax(rs))
    lag = int(lags[best_i])
    r = float(rs[best_i])
    if lag >= 0:
        a, b = x[: len(x) - lag or None], y[lag:]
    else:
        a, b = x[-lag:], y[: len(y) + lag or None]
    n = min(len(a), len(b))
    ne = n_eff_ar1(n, acf1(a[:n]), acf1(b[:n]))
    p_eff = p_from_r(r, ne)
    same_r = float("nan")
    if float(np.std(x)) and float(np.std(y)):
        same_r = float(pearsonr(x, y)[0])
    return {
        "best_lag": lag,
        "best_r": r,
        "p_eff_AR1": p_eff,
        "same_lag_r": same_r,
        "tp53_plateau": float(np.mean(x)),
        "mdm2_plateau": float(np.mean(y)),
        "lag_gt0_pos": bool(lag > 0 and r > 0),
        "induction_strict": bool(
            lag > 0 and r > 0 and math.isfinite(p_eff) and p_eff < 0.05
        ),
    }


def per_run_phase_from_block_samples(
    block_samples: np.ndarray,
    place_names: list[str],
    pair: tuple[str, str],
    *,
    warm: int = 20,
    max_lag: int = 40,
) -> list[dict[str, Any]]:
    """``block_samples`` shape (n_runs, n_blocks+1, n_places)."""
    a_name, b_name = pair
    try:
        ia = place_names.index(a_name)
        ib = place_names.index(b_name)
    except ValueError as exc:
        raise ValueError(f"NF QC places missing from model: {pair}") from exc

    out: list[dict[str, Any]] = []
    for run in range(block_samples.shape[0]):
        x = block_samples[run, :, ia]
        y = block_samples[run, :, ib]
        ph = phase_from_series(x, y, warm=warm, max_lag=max_lag)
        ph["run"] = int(run)
        out.append(ph)
    return out


def summarize_nf_passrate(per_run: list[dict[str, Any]]) -> dict[str, Any]:
    ok = [r for r in per_run if "error" not in r]
    n = len(per_run)
    n_ok = len(ok)
    if not ok:
        return {"n_runs": n, "n_ok": 0, "error": "no valid runs"}

    strict = [r for r in ok if r.get("induction_strict")]
    lag_pos = [r for r in ok if r.get("lag_gt0_pos")]
    lags = [int(r["best_lag"]) for r in ok if math.isfinite(r.get("best_r", float("nan")))]
    rs = [float(r["best_r"]) for r in ok if math.isfinite(r.get("best_r", float("nan")))]

    return {
        "n_runs": n,
        "n_ok": n_ok,
        "n_induction_strict": len(strict),
        "rate_induction_strict": len(strict) / n_ok,
        "n_lag_gt0_pos": len(lag_pos),
        "rate_lag_gt0_pos": len(lag_pos) / n_ok,
        "best_lag_median": float(np.median(lags)) if lags else float("nan"),
        "best_lag_mean": float(np.mean(lags)) if lags else float("nan"),
        "best_r_median": float(np.median(rs)) if rs else float("nan"),
        "best_r_mean": float(np.mean(rs)) if rs else float("nan"),
    }

"""Sparse diffusion / random-walk-with-restart on a (sub)graph CSR."""

from __future__ import annotations

import numpy as np
from scipy import sparse


def row_normalize_csr(adj: sparse.csr_matrix) -> sparse.csr_matrix:
    """Return row-stochastic transition matrix (zero rows stay zero)."""
    adj = adj.tocsr()
    out = adj.copy().astype(np.float64)
    sums = np.asarray(out.sum(axis=1)).ravel()
    for i in range(out.shape[0]):
        lo, hi = out.indptr[i], out.indptr[i + 1]
        if sums[i] > 0:
            out.data[lo:hi] /= sums[i]
    return out


def random_walk_restart(
    adj: sparse.csr_matrix,
    seed_indices: list[int],
    *,
    alpha: float = 0.3,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> np.ndarray:
    """
    RWR stationary scores: r = (1-alpha)*W.T*r + alpha*s, s uniform on seeds.

    ``adj`` should be nonnegative (e.g. unweighted subgraph). Uses power iteration.
    """
    n = adj.shape[0]
    W = row_normalize_csr(adj)
    s = np.zeros(n, dtype=np.float64)
    for j in seed_indices:
        s[j] += 1.0
    if s.sum() > 0:
        s /= s.sum()
    r = s.copy()
    WT = W.T.tocsr()
    for _ in range(max_iter):
        r_new = (1.0 - alpha) * (WT @ r) + alpha * s
        delta = float(np.linalg.norm(r_new - r, ord=1))
        r = r_new
        if delta < tol:
            break
    return r


def scores_to_named_tsv(
    path: str,
    scores: np.ndarray,
    index_to_name: list[str],
    new_to_old: dict[int, int] | None = None,
) -> None:
    """Write ``node\tscore``; if ``new_to_old`` given, map subgraph indices to global names."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("node\tscore\n")
        if new_to_old is None:
            for i, sc in enumerate(scores):
                f.write(f"{index_to_name[i]}\t{sc:g}\n")
        else:
            for j, sc in enumerate(scores):
                gi = new_to_old[j]
                f.write(f"{index_to_name[gi]}\t{sc:g}\n")

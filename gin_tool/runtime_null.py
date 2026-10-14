"""Runtime random-seed null (k-hop + RWR) for debias and benchmarks."""

from __future__ import annotations

import random
import time
from typing import Any, Sequence

import numpy as np
from scipy import sparse

from gin_tool.io import induced_subgraph_csr, khop_nodes_csr
from gin_tool.node_lists import UNKNOWN_SEED_LAYER, default_k_for_layer
from gin_tool.propagate import random_walk_restart


def run_runtime_null_replicate(
    sym: sparse.csr_matrix,
    index_to_name: Sequence[str],
    pool_indices: Sequence[int],
    k: int,
    n_seeds: int,
    rng: random.Random,
    *,
    rwr_alpha: float = 0.3,
    max_sub_nodes: int = 25_000,
    rwr_sample_n: int = 50,
) -> dict[str, Any]:
    """
    One runtime null draw: sample ``n_seeds`` from ``pool_indices``, k-hop, RWR.

    Returns timing, subgraph size, and optional pooled score sample for ECDF.
    """
    t0 = time.perf_counter()
    pool = list(pool_indices)
    if not pool:
        return {
            "ok": False,
            "reason": "empty_pool",
            "elapsed_sec": round(time.perf_counter() - t0, 6),
            "subgraph_nodes": 0,
            "rwr_scores_sample": None,
        }
    n_pick = min(int(n_seeds), len(pool))
    seeds = rng.sample(pool, n_pick)
    nodes = khop_nodes_csr(sym, seeds, int(k))
    sub_sym, old_to_new, _ = induced_subgraph_csr(sym, nodes)
    n_sub = int(sub_sym.shape[0])
    if n_sub == 0:
        return {
            "ok": False,
            "reason": "empty_subgraph",
            "elapsed_sec": round(time.perf_counter() - t0, 6),
            "subgraph_nodes": 0,
            "rwr_scores_sample": None,
        }
    if n_sub > max_sub_nodes:
        return {
            "ok": False,
            "reason": "n_sub_exceeds_max",
            "elapsed_sec": round(time.perf_counter() - t0, 6),
            "subgraph_nodes": n_sub,
            "rwr_scores_sample": None,
        }
    sub_seeds = [old_to_new[int(s)] for s in seeds if int(s) in old_to_new]
    if not sub_seeds:
        return {
            "ok": False,
            "reason": "no_seeds_in_subgraph",
            "elapsed_sec": round(time.perf_counter() - t0, 6),
            "subgraph_nodes": n_sub,
            "rwr_scores_sample": None,
        }
    scores = random_walk_restart(sub_sym, sub_seeds, alpha=rwr_alpha)
    n_samp = min(int(rwr_sample_n), n_sub)
    if n_sub > n_samp:
        idx = rng.sample(range(n_sub), n_samp)
    else:
        idx = list(range(n_sub))
    sample = scores[np.array(sorted(idx), dtype=int)]
    return {
        "ok": True,
        "reason": None,
        "elapsed_sec": round(time.perf_counter() - t0, 6),
        "subgraph_nodes": n_sub,
        "picked_seeds": n_pick,
        "rwr_scores_sample": sample,
    }


def allocate_replicates_by_layer(
    layer_counts: dict[str, int],
    n_replicates: int,
) -> dict[str, int]:
    """
    Split ``n_replicates`` across layers proportional to user seed counts.

    Each layer with a positive count receives at least one replicate when possible.
    """
    layers = {ly: int(c) for ly, c in layer_counts.items() if int(c) > 0}
    if not layers:
        return {}
    n_layers = len(layers)
    if n_replicates <= 0:
        return {ly: 0 for ly in layers}

    if n_replicates < n_layers:
        top = sorted(layers.keys(), key=lambda ly: (-layers[ly], ly))[:n_replicates]
        return {ly: (1 if ly in top else 0) for ly in layers}

    alloc = {ly: 1 for ly in layers}
    remaining = n_replicates - n_layers
    total = sum(layers.values())
    extras_float = {ly: remaining * layers[ly] / total for ly in layers}
    extras_floor = {ly: int(extras_float[ly]) for ly in layers}
    for ly in layers:
        alloc[ly] += extras_floor[ly]
    rem2 = n_replicates - sum(alloc.values())
    order = sorted(
        layers.keys(),
        key=lambda ly: (-(extras_float[ly] - extras_floor[ly]), ly),
    )
    for ly in order[:rem2]:
        alloc[ly] += 1
    return alloc


def runtime_null_pool_stratified(
    sym: sparse.csr_matrix,
    index_to_name: Sequence[str],
    buckets: dict[str, list[int]],
    layer_counts: dict[str, int],
    n_seeds: int,
    *,
    n_replicates: int = 25,
    rwr_alpha: float = 0.3,
    max_sub_nodes: int = 25_000,
    base_seed: int = 42,
    global_pool: Sequence[int] | None = None,
) -> tuple[list[np.ndarray], list[str], int, int]:
    """
    Runtime null replicates stratified by layers present in the user's seed list.

    Returns ``(score_chunks, notes, n_used, n_skipped)``.
    """
    alloc = allocate_replicates_by_layer(layer_counts, n_replicates)
    notes: list[str] = []
    notes.append(f"stratified runtime null replicate allocation: {alloc}")
    chunks: list[np.ndarray] = []
    used = 0
    skipped = 0
    rep_global = 0
    for layer, n_rep in sorted(alloc.items()):
        if n_rep <= 0:
            continue
        if layer == UNKNOWN_SEED_LAYER:
            pool = list(global_pool or [])
            k_layer = 2
            layer_label = "unlisted/global"
        else:
            pool = list(buckets.get(layer, []))
            k_layer = default_k_for_layer(layer)
            layer_label = layer
        if not pool:
            notes.append(f"runtime null: empty pool for layer {layer_label!r}; skipped {n_rep} reps.")
            skipped += n_rep
            continue
        for _ in range(n_rep):
            rng = random.Random(
                base_seed + rep_global * 1_000_003 + hash(layer) % (2**20) + n_seeds
            )
            rep_global += 1
            rec = run_runtime_null_replicate(
                sym,
                index_to_name,
                pool,
                k_layer,
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
    return chunks, notes, used, skipped

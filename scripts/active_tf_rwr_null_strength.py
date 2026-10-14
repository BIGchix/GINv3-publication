#!/usr/bin/env python3
"""Directed RWR active-TF nulls with strengthened designs A / B / C.

Shared observed statistic: RWR on directed out-envelope (graph_hops=2*real_hops)
from true protein seeds → TF scores.

Null strengths (``--null-strength``):

* ``strata``   — stratified seed swap + per_null envelope (baseline)
* ``env_size`` — (A) strata + reject envelopes outside [env_lo,env_hi]*|E_obs|
* ``cohesion`` — (B, **pipeline default**) (A) + **localized** stratified seeds:
                 pick a random center protein, sample strata-matched seeds from
                 its directed out-envelope (forces co-localization), then keep
                 only if the seed-set envelope also passes the env-size band
* ``rewire``   — (C, optional) keep true seeds + observed envelope;
                 degree-preserving directed edge swaps on the induced subgraph,
                 then RWR. Slower (~6× B) but stricter / more accurate.

Pipeline entry (default B): ``scripts/select_active_tfs.py``.

Example (default B / cohesion):
  python scripts/active_tf_rwr_null_strength.py \\
    --gene-file seedLists/p53_pi3k_curated/genes.txt \\
    --n-rwr-null 1000 --n-jobs 16 \\
    --out-prefix benchmarks/subgraph-speed/p53_rwr_null_B_cohesion

Optional C / rewire (slower, stricter):
  python scripts/active_tf_rwr_null_strength.py \\
    --gene-file seedLists/p53_pi3k_curated/genes.txt \\
    --null-strength rewire --n-rwr-null 1000 --n-jobs 16 \\
    --rewire-sweeps 0.1 \\
    --out-prefix benchmarks/subgraph-speed/p53_rwr_null_C_rewire
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path

import numpy as np
from scipy import sparse

from gin_tool.active_tf import (
    bh_fdr,
    build_tf_regulon_map,
    gin_node_to_gene_symbol,
    load_tf_gin_node_map,
    load_tf_symbol_pool,
    seed_nodes_to_gene_symbols,
)
from gin_tool.gin_defaults import DEFAULT_NODE_LIST_DIR, DEFAULT_SIFI
from gin_tool.io import (
    build_indexed_edges,
    collect_unique_node_names,
    directed_csr_from_indexed,
    induced_subgraph_csr,
    khop_nodes_csr,
)
from gin_tool.mapping import is_virtual_gin_node
from gin_tool.node_lists import load_node_list_dir
from gin_tool.propagate import random_walk_restart
from gin_tool.subgraph_pipeline import build_seed_manifest_from_files
from gin_tool.subgraph_pipeline.metabolite_id_resolve import collect_chem_nodes_from_sifi
from gin_tool.subgraph_pipeline.seed_gin_mapping import (
    collect_walk_input_nodes,
    map_union_seeds,
)

_WORKER: dict = {}


def _gene_symbol_filter(raw) -> set[str]:
    out: set[str] = set()
    for s in raw:
        if s.startswith(("SID", "CID", "CHEBI", "CAS", "KEGG", "HMDB", "UC_")):
            continue
        if s[:1].isdigit():
            continue
        out.add(s)
    return out


def _quantile_bin(values: dict, n_bins: int) -> dict:
    if not values:
        return {}
    arr = np.array(sorted(values.values()), dtype=float)
    edges = np.unique(np.quantile(arr, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 2:
        return {k: 0 for k in values}
    out: dict = {}
    for k, v in values.items():
        b = int(np.searchsorted(edges, v, side="right") - 1)
        out[k] = max(0, min(b, len(edges) - 2))
    return out


def _stratified_sample(rng, key_to_genes: dict, key_counts: dict) -> set:
    picked: list = []
    for key, need in key_counts.items():
        pool = key_to_genes.get(key, [])
        if need <= 0:
            continue
        if need >= len(pool):
            picked.extend(pool)
            continue
        idx = rng.choice(len(pool), size=need, replace=False)
        picked.extend(pool[i] for i in idx)
    need_total = sum(key_counts.values())
    if len(picked) < need_total:
        used = set(picked)
        rest = [g for genes in key_to_genes.values() for g in genes if g not in used]
        extra = need_total - len(picked)
        if rest and extra > 0:
            idx = rng.choice(len(rest), size=min(extra, len(rest)), replace=False)
            picked.extend(rest[i] for i in np.atleast_1d(idx))
    return set(picked)


def _out_neighbors(adj, idx: int) -> set[int]:
    lo = adj.indptr[idx]
    hi = adj.indptr[idx + 1]
    return {int(v) for v in adj.indices[lo:hi]}


def seed_mutual_reach_cohesion(adj, seed_idx, hops: int) -> float:
    """Fraction of other seeds inside each seed's directed out-ball.

    Uses seed-local expansion (only walks from seed out-neighbors), not a full
    multi-source khop over the whole GIN — important for rejection sampling speed.
    """
    seeds = [int(i) for i in seed_idx]
    if len(seeds) < 2:
        return 0.0
    seed_set = set(seeds)
    n1 = [_out_neighbors(adj, s) for s in seeds]
    acc = 0.0
    for i, s in enumerate(seeds):
        reach = set(n1[i])
        if hops >= 2:
            extra: set[int] = set()
            for v in n1[i]:
                extra |= _out_neighbors(adj, v)
            reach |= extra
        if hops > 2:
            frontier = set(reach)
            visited = set(reach) | {s}
            for _ in range(hops - 2):
                nxt: set[int] = set()
                for u in frontier:
                    for v in _out_neighbors(adj, u):
                        if v in visited:
                            continue
                        visited.add(v)
                        nxt.add(v)
                if not nxt:
                    break
                reach |= nxt
                frontier = nxt
        others = seed_set - {s}
        acc += len(others & reach) / len(others)
    return acc / len(seeds)


def seed_out_jaccard_cohesion(
    adj,
    seed_idx,
    *,
    max_pairs: int = 200,
    rng: np.random.Generator | None = None,
) -> float:
    """Mean Jaccard overlap of hop-1 out-neighborhoods over seed pairs."""
    seeds = [int(i) for i in seed_idx]
    if len(seeds) < 2:
        return 0.0
    neigh = [_out_neighbors(adj, i) for i in seeds]
    pairs = [(i, j) for i in range(len(seeds)) for j in range(i + 1, len(seeds))]
    if len(pairs) > max_pairs:
        if rng is None:
            rng = np.random.default_rng(0)
        pick = rng.choice(len(pairs), size=max_pairs, replace=False)
        pairs = [pairs[int(k)] for k in pick]
    acc = 0.0
    for i, j in pairs:
        a, b = neigh[i], neigh[j]
        if not a and not b:
            continue
        inter = len(a & b)
        union = len(a | b)
        acc += (inter / union) if union else 0.0
    return acc / max(1, len(pairs))


def _score_tfs_from_sub(scores, new_to_old, index_to_name, name_to_tf, tf_symbols):
    out: dict[str, float] = {}
    for j, sc in enumerate(scores):
        old = int(new_to_old[j])
        name = index_to_name[old]
        bare = name_to_tf.get(name) or gin_node_to_gene_symbol(name)
        if bare not in tf_symbols:
            continue
        out[bare] = max(out.get(bare, 0.0), float(sc))
    return out


def _rwr_on_seed_envelope(adj_use, seed_idx, graph_hops, max_sub_nodes, alpha):
    if not seed_idx:
        return None
    env = {int(i) for i in khop_nodes_csr(adj_use, seed_idx, graph_hops)}
    if not env or len(env) > max_sub_nodes:
        return None
    sub, old_to_new, new_to_old = induced_subgraph_csr(adj_use, env)
    sub_seeds = [old_to_new[i] for i in seed_idx if i in old_to_new]
    if not sub_seeds:
        return None
    scores = random_walk_restart(sub, sub_seeds, alpha=alpha)
    return scores, new_to_old, int(sub.shape[0]), int(sub.nnz)


def _edge_list_from_csr(adj):
    adj = adj.tocsr()
    rows: list[int] = []
    cols: list[int] = []
    for u in range(adj.shape[0]):
        lo = adj.indptr[u]
        hi = adj.indptr[u + 1]
        for v in adj.indices[lo:hi]:
            if u == int(v):
                continue
            rows.append(u)
            cols.append(int(v))
    return np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32)


def _csr_from_edges(n, rows, cols):
    data = np.ones(len(rows), dtype=np.float64)
    return sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def degree_preserving_directed_rewire(rows, cols, n_nodes, *, n_swaps, rng):
    """Directed configuration-style edge swaps; preserves in/out degrees."""
    rows = rows.copy()
    cols = cols.copy()
    m = len(rows)
    if m < 2 or n_swaps <= 0:
        return rows, cols, 0
    edge_set = {(int(rows[i]), int(cols[i])) for i in range(m)}
    ok = 0
    attempts = 0
    max_attempts = max(n_swaps * 20, n_swaps + 1)
    while ok < n_swaps and attempts < max_attempts:
        attempts += 1
        i, j = rng.integers(0, m, size=2)
        if i == j:
            continue
        a, b = int(rows[i]), int(cols[i])
        c, d = int(rows[j]), int(cols[j])
        if a == c or b == d or a == d or c == b:
            continue
        if (a, d) in edge_set or (c, b) in edge_set:
            continue
        edge_set.discard((a, b))
        edge_set.discard((c, d))
        edge_set.add((a, d))
        edge_set.add((c, b))
        cols[i] = d
        cols[j] = b
        ok += 1
    return rows, cols, ok


def _init_seed_worker(adj_use, index_to_name, name_to_tf, track, graph_hops, max_sub_nodes, alpha):
    _WORKER.clear()
    _WORKER.update(
        mode="seed",
        adj=adj_use,
        index_to_name=index_to_name,
        name_to_tf=name_to_tf,
        track_set=set(track),
        graph_hops=graph_hops,
        max_sub_nodes=max_sub_nodes,
        alpha=alpha,
    )


def _seed_worker_one(null_idx):
    w = _WORKER
    pack = _rwr_on_seed_envelope(
        w["adj"], null_idx, w["graph_hops"], w["max_sub_nodes"], w["alpha"]
    )
    if pack is None:
        return None
    sc, nto, n_env, _ = pack
    null_tf = _score_tfs_from_sub(
        sc, nto, w["index_to_name"], w["name_to_tf"], w["track_set"]
    )
    return null_tf, n_env


def _init_rewire_worker(
    base_rows,
    base_cols,
    n_nodes,
    sub_seeds,
    new_to_old,
    index_to_name,
    name_to_tf,
    track,
    alpha,
    n_swaps,
):
    _WORKER.clear()
    _WORKER.update(
        mode="rewire",
        base_rows=base_rows,
        base_cols=base_cols,
        n_nodes=n_nodes,
        sub_seeds=sub_seeds,
        new_to_old=new_to_old,
        index_to_name=index_to_name,
        name_to_tf=name_to_tf,
        track_set=set(track),
        alpha=alpha,
        n_swaps=n_swaps,
    )


def _rewire_worker_one(rep_seed):
    w = _WORKER
    rng = np.random.default_rng(rep_seed)
    rows, cols, n_ok = degree_preserving_directed_rewire(
        w["base_rows"],
        w["base_cols"],
        w["n_nodes"],
        n_swaps=w["n_swaps"],
        rng=rng,
    )
    if n_ok < max(1, w["n_swaps"] // 100):
        return None
    sub = _csr_from_edges(w["n_nodes"], rows, cols)
    sc = random_walk_restart(sub, w["sub_seeds"], alpha=w["alpha"])
    null_tf = _score_tfs_from_sub(
        sc, w["new_to_old"], w["index_to_name"], w["name_to_tf"], w["track_set"]
    )
    return null_tf, w["n_nodes"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-file", required=True)
    ap.add_argument("--real-hops", type=int, default=1)
    ap.add_argument("--n-rwr-null", type=int, default=1000)
    ap.add_argument("--rwr-alpha", type=float, default=0.3)
    ap.add_argument("--n-bins", type=int, default=5)
    ap.add_argument("--max-sub-nodes", type=int, default=120000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--top-n", type=int, default=80)
    ap.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help=(
            "Parallel worker processes for null replicates (ProcessPool + fork). "
            "Increase to use more CPUs, e.g. --n-jobs 16 or 32. Alias: --cpus."
        ),
    )
    ap.add_argument(
        "--cpus",
        dest="n_jobs",
        type=int,
        help="Alias for --n-jobs (number of CPU worker processes).",
    )
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument(
        "--null-strength",
        choices=("strata", "env_size", "cohesion", "rewire"),
        default="cohesion",
        help=(
            "Null design: cohesion=B (default), rewire=C (slower/stricter). "
            "Aliases via --method: b|cohesion, c|rewire."
        ),
    )
    ap.add_argument(
        "--method",
        choices=("a", "b", "c", "strata", "env_size", "cohesion", "rewire"),
        default=None,
        help="Alias for --null-strength: a=env_size, b=cohesion, c=rewire.",
    )
    ap.add_argument(
        "--emp-p-cutoff",
        type=float,
        default=0.05,
        help="Write {out-prefix}_emp05.txt with symbols at emp_p<=cutoff (default 0.05)",
    )
    ap.add_argument("--env-lo", type=float, default=0.5)
    ap.add_argument("--env-hi", type=float, default=1.5)
    ap.add_argument(
        "--cohesion-lo",
        type=float,
        default=0.5,
        help="B: keep nulls with cohesion >= lo * obs (one-sided by default)",
    )
    ap.add_argument(
        "--cohesion-hi",
        type=float,
        default=0.0,
        help="B: if >0, also require cohesion <= hi * obs; 0 disables upper bound",
    )
    ap.add_argument(
        "--sample-tries-factor",
        type=int,
        default=80,
        help="Max stratified draws = n_rwr_null * factor when rejecting (A/B)",
    )
    ap.add_argument(
        "--sample-progress-every",
        type=int,
        default=500,
        help="Log sampling progress every N tries",
    )
    ap.add_argument(
        "--rewire-sweeps",
        type=float,
        default=0.1,
        help="Target successful swaps ≈ sweeps * n_edges (C); 0.1≈minutes-scale for GIN envelope",
    )
    ap.add_argument(
        "--cohesion-metric",
        choices=("mutual_reach", "out_jaccard"),
        default="mutual_reach",
        help="B: mutual_reach=fraction of other seeds in out-khop; out_jaccard=hop1 overlap",
    )
    args = ap.parse_args()

    alias = {
        "a": "env_size",
        "b": "cohesion",
        "c": "rewire",
        "strata": "strata",
        "env_size": "env_size",
        "cohesion": "cohesion",
        "rewire": "rewire",
    }
    if args.method is not None:
        args.null_strength = alias[args.method]

    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    print("[1] load directed GIN", flush=True)
    ed = build_indexed_edges(str(DEFAULT_SIFI))
    adj = directed_csr_from_indexed(ed)
    index_to_name = ed.index_to_name
    name_to_index = ed.name_to_index

    manifest = build_seed_manifest_from_files(
        gene_file=args.gene_file, run_id="rwr_abc"
    )
    gin_nodes = set(collect_unique_node_names(str(DEFAULT_SIFI)))
    pools = load_node_list_dir(DEFAULT_NODE_LIST_DIR)
    met_pool = set(pools.get("metabolite", [])) | collect_chem_nodes_from_sifi(
        str(DEFAULT_SIFI)
    )
    seed_kinds = {g: "gene" for g in manifest.union_genes}
    recs, _ = map_union_seeds(
        sorted(manifest.union_genes),
        gin_nodes.__contains__,
        seed_kinds=seed_kinds,
        metabolite_pool=met_pool,
        allow_online_metabolite=False,
    )
    walk_anchors = collect_walk_input_nodes(recs)
    protein_anchors = [
        r.protein_node for r in recs if r.protein_node and r.protein_node in name_to_index
    ]
    seed_symbols = _gene_symbol_filter(seed_nodes_to_gene_symbols(walk_anchors))
    print(
        f"  protein_anchors={len(protein_anchors)} seed_symbols={len(seed_symbols)}",
        flush=True,
    )
    seed_idx = [name_to_index[a] for a in protein_anchors]
    if not seed_idx:
        raise SystemExit("no protein seed anchors resolved in GIN")

    graph_hops = max(1, 2 * int(args.real_hops))
    print(f"[2] observed envelope hops={graph_hops}", flush=True)
    obs_env = {int(i) for i in khop_nodes_csr(adj, seed_idx, graph_hops)}
    obs_env_n = len(obs_env)
    print(f"  observed_out_envelope_nodes={obs_env_n}", flush=True)
    if obs_env_n > args.max_sub_nodes:
        raise SystemExit(f"envelope {obs_env_n} > max_sub_nodes {args.max_sub_nodes}")

    if args.cohesion_metric == "mutual_reach":
        obs_cohesion = seed_mutual_reach_cohesion(adj, seed_idx, graph_hops)
    else:
        obs_cohesion = seed_out_jaccard_cohesion(adj, seed_idx, rng=rng)
    print(
        f"  observed_cohesion[{args.cohesion_metric}]={obs_cohesion:.6f}",
        flush=True,
    )

    print("[3] strata pool", flush=True)
    tf_gin_map = load_tf_gin_node_map(DEFAULT_NODE_LIST_DIR)
    tf_pool = load_tf_symbol_pool(DEFAULT_NODE_LIST_DIR)
    tf_symbols = set(tf_pool)
    regulons, universe = build_tf_regulon_map(
        str(DEFAULT_SIFI), tf_symbols, tf_gin_by_bare=tf_gin_map
    )
    gene_tf_deg = Counter()
    for genes in regulons.values():
        for g in genes:
            gene_tf_deg[g] += 1

    gene_to_prot_idx: dict[str, int] = {}
    outdeg: dict[str, int] = {}
    for i, name in enumerate(index_to_name):
        if is_virtual_gin_node(name):
            continue
        if not name.endswith(".protein"):
            continue
        bare = gin_node_to_gene_symbol(name)
        if not bare:
            continue
        prefer = gene_to_prot_idx.get(bare)
        if prefer is None or (adj.indptr[i + 1] - adj.indptr[i]) > (
            adj.indptr[prefer + 1] - adj.indptr[prefer]
        ):
            gene_to_prot_idx[bare] = i
            outdeg[bare] = int(adj.indptr[i + 1] - adj.indptr[i])

    gin_n = len(index_to_name)
    pool_genes = {
        g: gene_to_prot_idx[g]
        for g in gene_to_prot_idx
        if g in universe or g in gene_tf_deg or g in seed_symbols
    }
    tf_bins = _quantile_bin(
        {g: float(gene_tf_deg.get(g, 0)) for g in pool_genes}, args.n_bins
    )
    od_bins = _quantile_bin(
        {g: float(outdeg.get(g, 0)) for g in pool_genes}, args.n_bins
    )
    key_to_genes: dict[tuple, list] = defaultdict(list)
    for g in pool_genes:
        key_to_genes[(tf_bins.get(g, 0), od_bins.get(g, 0))].append(g)

    obs_counts: dict[tuple, int] = Counter()
    for g in seed_symbols:
        if g in pool_genes:
            obs_counts[(tf_bins.get(g, 0), od_bins.get(g, 0))] += 1
    n_pick = sum(obs_counts.values())
    if n_pick <= 0:
        raise SystemExit("no seed symbols in strata pool")

    name_to_tf = {v: k for k, v in tf_gin_map.items()}
    print(f"[4] observed RWR (null_strength={args.null_strength})", flush=True)
    obs_pack = _rwr_on_seed_envelope(
        adj, seed_idx, graph_hops, args.max_sub_nodes, args.rwr_alpha
    )
    if obs_pack is None:
        raise SystemExit("observed RWR failed")
    scores, new_to_old_obs, n_sub, n_edges = obs_pack
    obs_tf = _score_tfs_from_sub(
        scores, new_to_old_obs, index_to_name, name_to_tf, tf_symbols
    )
    track = sorted(obs_tf.keys())
    print(f"  n_tf={len(track)} n_sub={n_sub} n_edges={n_edges}", flush=True)

    ge = {tf: 0 for tf in track}
    env_sizes: list[int] = []
    cohesions: list[float] = []
    ok = 0
    skipped = 0
    n_jobs = max(1, int(args.n_jobs))
    sampling_meta: dict = {}

    def consume(result) -> bool:
        nonlocal ok, skipped
        if result is None:
            skipped += 1
            return False
        null_tf, n_env = result
        env_sizes.append(n_env)
        for tf in ge:
            if null_tf.get(tf, 0.0) >= obs_tf.get(tf, 0.0):
                ge[tf] += 1
        ok += 1
        return True

    if args.null_strength == "rewire":
        sub, old_to_new, new_to_old = induced_subgraph_csr(adj, obs_env)
        sub_seeds = [old_to_new[i] for i in seed_idx if i in old_to_new]
        rows, cols = _edge_list_from_csr(sub)
        n_swaps = max(1, int(args.rewire_sweeps * len(rows)))
        print(
            f"  rewire: n_edges={len(rows)} target_swaps={n_swaps} n_jobs={n_jobs}",
            flush=True,
        )
        rep_seeds = rng.integers(0, 2147483647, size=args.n_rwr_null)
        initargs = (
            rows,
            cols,
            int(sub.shape[0]),
            sub_seeds,
            new_to_old,
            index_to_name,
            name_to_tf,
            track,
            args.rwr_alpha,
            n_swaps,
        )
        t_r = time.time()
        done = 0
        total = args.n_rwr_null
        ctx = get_context("fork")
        with ProcessPoolExecutor(
            max_workers=n_jobs,
            mp_context=ctx,
            initializer=_init_rewire_worker,
            initargs=initargs,
        ) as pool:
            futs = [pool.submit(_rewire_worker_one, int(rs)) for rs in rep_seeds]
            for fut in as_completed(futs):
                consume(fut.result())
                done += 1
                if done % 50 == 0 or done == total:
                    print(
                        f"  rewire {done}/{total} ok={ok} skip={skipped} "
                        f"({time.time() - t_r:.1f}s)",
                        flush=True,
                    )
        sampling_meta = {
            "mode": "rewire",
            "n_swaps": n_swaps,
            "n_edges": len(rows),
        }
    else:
        # strata / env_size / cohesion
        env_lo = args.env_lo * obs_env_n
        env_hi = args.env_hi * obs_env_n
        need_env = args.null_strength in ("env_size", "cohesion")
        localized = args.null_strength == "cohesion"
        max_tries = args.n_rwr_null * args.sample_tries_factor
        print(
            f"  sample nulls strength={args.null_strength} "
            f"env_band=[{env_lo:.0f},{env_hi:.0f}] localized={localized}",
            flush=True,
        )
        null_jobs: list[list[int]] = []
        tries = 0
        reject_env = 0
        reject_local = 0
        reject_short = 0
        prot_idx_pool = list(pool_genes.values())
        while len(null_jobs) < args.n_rwr_null and tries < max_tries:
            tries += 1
            if localized:
                if not prot_idx_pool:
                    reject_local += 1
                    continue
                center = int(rng.choice(prot_idx_pool))
                local_nodes = {int(i) for i in khop_nodes_csr(adj, [center], graph_hops)}
                local_genes = []
                for ni in local_nodes:
                    name = index_to_name[ni]
                    if not name.endswith(".protein"):
                        continue
                    bare = gin_node_to_gene_symbol(name)
                    if bare and bare in pool_genes:
                        local_genes.append(bare)
                if len(local_genes) < n_pick:
                    reject_local += 1
                    continue
                local_key_to_genes: dict[tuple, list] = defaultdict(list)
                for g in local_genes:
                    local_key_to_genes[(tf_bins.get(g, 0), od_bins.get(g, 0))].append(g)
                need = dict(obs_counts)
                # fall back to global keys if local strata too thin
                for key, cnt in list(need.items()):
                    if len(local_key_to_genes.get(key, [])) < cnt:
                        reject_short += 1
                        need = None  # type: ignore[assignment]
                        break
                if need is None:
                    continue
                null_syms = _stratified_sample(rng, local_key_to_genes, need)
            else:
                null_syms = _stratified_sample(rng, key_to_genes, obs_counts)

            if len(null_syms) < n_pick:
                reject_short += 1
                continue
            null_idx = [gene_to_prot_idx[g] for g in null_syms if g in gene_to_prot_idx]
            if len(null_idx) < 2:
                reject_short += 1
                continue

            if need_env:
                env_n = len(khop_nodes_csr(adj, null_idx, graph_hops))
                if env_n < env_lo or env_n > env_hi:
                    reject_env += 1
                    continue
                if localized:
                    if args.cohesion_metric == "mutual_reach":
                        cval = seed_mutual_reach_cohesion(adj, null_idx, graph_hops)
                    else:
                        cval = seed_out_jaccard_cohesion(adj, null_idx, rng=rng)
                    if cval < args.cohesion_lo * obs_cohesion:
                        reject_env += 1
                        continue
                    if args.cohesion_hi > 0 and cval > args.cohesion_hi * obs_cohesion:
                        reject_env += 1
                        continue
                    cohesions.append(cval)

            null_jobs.append(null_idx)
            if (
                args.sample_progress_every
                and tries % args.sample_progress_every == 0
            ):
                print(
                    f"  accepted {len(null_jobs)} tries={tries} "
                    f"rej_env={reject_env} rej_local={reject_local} "
                    f"rej_short={reject_short}",
                    flush=True,
                )

        print(
            f"  sampling tries={tries} accepted={len(null_jobs)} "
            f"rej_env={reject_env} rej_local={reject_local}",
            flush=True,
        )
        if len(null_jobs) < args.n_rwr_null:
            print(
                f"  WARN: only accepted {len(null_jobs)} after {tries} tries",
                flush=True,
            )
        sampling_meta = {
            "tries": tries,
            "reject_env": reject_env,
            "reject_local": reject_local,
            "reject_short": reject_short,
            "accepted": len(null_jobs),
            "localized": localized,
        }

        print(f"  dispatch {len(null_jobs)} RWR jobs n_jobs={n_jobs}", flush=True)
        t_r = time.time()
        initargs = (
            adj,
            index_to_name,
            name_to_tf,
            track,
            graph_hops,
            args.max_sub_nodes,
            args.rwr_alpha,
        )
        ctx = get_context("fork")
        with ProcessPoolExecutor(
            max_workers=n_jobs,
            mp_context=ctx,
            initializer=_init_seed_worker,
            initargs=initargs,
        ) as pool:
            futs = [pool.submit(_seed_worker_one, nj) for nj in null_jobs]
            done = 0
            total = len(futs)
            for fut in as_completed(futs):
                consume(fut.result())
                done += 1
                if done % 50 == 0 or done == total:
                    med = float(np.median(env_sizes)) if env_sizes else 0.0
                    print(
                        f"  null {done}/{total} ok={ok} skip={skipped} "
                        f"env_med={med:.0f} ({time.time() - t_r:.1f}s)",
                        flush=True,
                    )

    n_rep = max(1, ok)
    records = []
    emp_p = []
    for tf in track:
        p = (ge[tf] + 1.0) / (n_rep + 1.0)
        emp_p.append(p)
        records.append(
            {
                "symbol": tf,
                "rwr_score": float(obs_tf.get(tf, 0.0)),
                "rwr_emp_p": float(p),
                "rwr_ge_count": int(ge[tf]),
            }
        )
    qs = bh_fdr(emp_p)
    for r, q in zip(records, qs):
        r["rwr_emp_q"] = float(q)

    ordered = sorted(records, key=lambda r: (r["rwr_emp_p"], -r["rwr_score"]))
    env_arr = np.array(env_sizes, dtype=np.int64) if env_sizes else np.array([], dtype=np.int64)

    focus_syms = (
        "MYC",
        "MDM2",
        "TP53",
        "TP73",
        "PTEN",
        "ATM",
        "CDK2",
        "MAX",
        "CUL4A",
        "CASP3",
        "CCND2",
        "E2F1",
        "SIRT1",
    )
    by = {r["symbol"]: r for r in ordered}
    focus_table = []
    for rank, sym in enumerate(focus_syms, start=1):
        r = by.get(sym)
        if not r:
            continue
        focus_table.append(
            {
                "symbol": sym,
                "emp_p": r["rwr_emp_p"],
                "rank": next(
                    (i + 1 for i, x in enumerate(ordered) if x["symbol"] == sym),
                    None,
                ),
                "score": r["rwr_score"],
                "ge": r["rwr_ge_count"],
            }
        )

    summary = {
        "method": f"directed_RWR_null_strength_{args.null_strength}",
        "null_strength": args.null_strength,
        "real_hops": args.real_hops,
        "graph_hops": graph_hops,
        "observed_envelope_nodes": obs_env_n,
        "observed_seed_out_jaccard": float(obs_cohesion)
        if args.cohesion_metric == "out_jaccard"
        else seed_out_jaccard_cohesion(adj, seed_idx, rng=rng),
        "cohesion_metric": args.cohesion_metric,
        "env_lo": args.env_lo,
        "env_hi": args.env_hi,
        "cohesion_lo": args.cohesion_lo,
        "cohesion_hi": args.cohesion_hi,
        "rewire_sweeps": args.rewire_sweeps if args.null_strength == "rewire" else None,
        "n_rwr_null": args.n_rwr_null,
        "n_jobs": n_jobs,
        "n_replicates_ok": ok,
        "n_replicates_skipped": skipped,
        "null_envelope_nodes_median": float(np.median(env_arr)) if env_arr.size else None,
        "null_envelope_nodes_min": int(env_arr.min()) if env_arr.size else None,
        "null_envelope_nodes_max": int(env_arr.max()) if env_arr.size else None,
        "null_cohesion_median": float(np.median(cohesions)) if cohesions else None,
        "p_le_0_05": sum(1 for r in ordered if r["rwr_emp_p"] <= 0.05),
        "q_le_0_05": sum(1 for r in ordered if r["rwr_emp_q"] <= 0.05),
        "n_tf_scored": len(ordered),
        "top_by_emp_p": [r["symbol"] for r in ordered[: args.top_n]],
        "top_records": ordered[: min(30, len(ordered))],
        "focus": focus_table,
        "sampling": sampling_meta,
        "elapsed_seconds": time.time() - t0,
    }

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(prefix) + "_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    Path(str(prefix) + "_full.json").write_text(
        json.dumps({"summary": summary, "records": ordered}, indent=2),
        encoding="utf-8",
    )
    emp_syms = [r["symbol"] for r in ordered if r["rwr_emp_p"] <= args.emp_p_cutoff]
    emp_path = Path(str(prefix) + "_emp05.txt")
    emp_path.write_text("\n".join(emp_syms) + ("\n" if emp_syms else ""), encoding="utf-8")
    meta = {
        "null_strength": args.null_strength,
        "emp_p_cutoff": args.emp_p_cutoff,
        "n_active": len(emp_syms),
        "n_rwr_null": args.n_rwr_null,
        "n_replicates_ok": ok,
        "source_summary": str(Path(str(prefix) + "_summary.json")),
        "pipeline_default": args.null_strength in ("cohesion", "rewire"),
        "note": "B/cohesion = pipeline default; C/rewire = optional stricter/slower",
    }
    Path(str(prefix) + "_emp05.meta.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    print(
        f"[done] {prefix}_summary.json emp05_n={len(emp_syms)} -> {emp_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Augment bipartite signaling walk with seed central-dogma self-loops."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass, field
from typing import Any

from gin_tool.mapping import is_virtual_gin_node
from gin_tool.signaling_subgraph import _edge_allowed, _out_neighbors_csr
from gin_tool.subgraph_pipeline.manifest import NodeCountStats, SeedSelfLoopQC
from gin_tool.subgraph_pipeline.seed_gin_mapping import CentralDogmaSeed
from gin_tool.subgraph_pipeline.signaling_bipartite import (
    SignalingBipartiteWalkResult,
    add_real_hop_path,
)
from gin_tool.subgraph_pipeline.stats import node_count_stats


@dataclass
class SeedSelfLoopRecord:
    symbol: str
    rna_node: str | None
    rna_self_loop: bool = False
    rna_self_loop_hops: int | None = None
    protein_self_loop: bool = False
    protein_self_loop_hops: int | None = None
    promoter_self_loop: bool = False
    promoter_self_loop_hops: int | None = None
    triplet_complete: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _bfs_self_loop_real_hops(
    s_idx: int,
    adj,
    index_to_name: list[str],
    input_idx: set[int],
    *,
    max_real_hops: int,
    max_real_nodes: int,
    filter_regulatory: bool,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
) -> tuple[int | None, dict[int, set[int]], bool]:
    real_dist: dict[int, int] = {s_idx: 0}
    real_parents: dict[int, set[int]] = {s_idx: set()}
    frontier_real = [s_idx]
    self_loop_hop: int | None = None

    for hop in range(max_real_hops):
        if not frontier_real or self_loop_hop is not None:
            break

        v_from: dict[int, set[int]] = defaultdict(set)
        for r in frontier_real:
            u_name = index_to_name[r]
            for v in _out_neighbors_csr(adj, r):
                v = int(v)
                v_name = index_to_name[v]
                if not is_virtual_gin_node(v_name):
                    continue
                if not _edge_allowed(
                    u_name,
                    v_name,
                    filter_regulatory=filter_regulatory,
                    active_tf_names=active_tf_names,
                    seed_symbols=seed_symbols,
                    seed_gin_names=seed_gin_names,
                ):
                    continue
                v_from[v].add(r)

        if not v_from:
            break

        next_preds: dict[int, set[int]] = defaultdict(set)
        new_hop = hop + 1
        for v, pred_reals in v_from.items():
            v_name = index_to_name[v]
            for r2 in _out_neighbors_csr(adj, v):
                r2 = int(r2)
                r2_name = index_to_name[r2]
                if is_virtual_gin_node(r2_name):
                    continue
                if not _edge_allowed(
                    v_name,
                    r2_name,
                    filter_regulatory=filter_regulatory,
                    active_tf_names=active_tf_names,
                    seed_symbols=seed_symbols,
                    seed_gin_names=seed_gin_names,
                ):
                    continue
                if r2 == s_idx and new_hop > 0:
                    self_loop_hop = new_hop
                    next_preds[r2].update(pred_reals)
                    continue
                old = real_dist.get(r2)
                if old is not None and old <= new_hop:
                    continue
                next_preds[r2].update(pred_reals)

        if self_loop_hop is not None:
            real_parents[s_idx] = set(next_preds.get(s_idx, set()))
            break

        if not next_preds:
            break

        frontier_real = []
        for r2, preds in next_preds.items():
            real_dist[r2] = new_hop
            real_parents[r2] = set(preds)
            frontier_real.append(r2)
            if len(real_dist) >= max_real_nodes:
                break

    hit_cap = len(real_dist) >= max_real_nodes and self_loop_hop is None
    return self_loop_hop, real_parents, hit_cap


def augment_walk_with_seed_self_loops(
    walk: SignalingBipartiteWalkResult,
    seed_records: list[CentralDogmaSeed],
    adj_dir,
    index_to_name: list[str],
    name_to_index: dict[str, int],
    *,
    max_real_hops: int,
    max_real_nodes: int,
    filter_regulatory: bool,
    active_tf_names: set[str],
    seed_symbols: set[str],
    seed_gin_names: set[str],
    self_loop_hops: int | None = None,
) -> tuple[SignalingBipartiteWalkResult, SeedSelfLoopQC]:
    hop_limit = self_loop_hops if self_loop_hops is not None else max(max_real_hops, 12)
    input_idx = {
        name_to_index[n]
        for rec in seed_records
        for n in rec.walk_anchor_nodes
        if n in name_to_index
    }

    loop_records: list[SeedSelfLoopRecord] = []
    hit_cap = False
    for rec in seed_records:
        if rec.kind == "metabolite":
            loop_records.append(SeedSelfLoopRecord(symbol=rec.symbol, rna_node=None, triplet_complete=False))
            continue
        lr = SeedSelfLoopRecord(
            symbol=rec.symbol,
            rna_node=rec.rna_node,
            triplet_complete=rec.triplet_complete,
        )

        for attr, node in (
            ("rna", rec.rna_node),
            ("protein", rec.protein_node),
            ("promoter", rec.promoter_node),
        ):
            if not node or node not in name_to_index:
                continue
            idx = name_to_index[node]
            hop, parents, capped = _bfs_self_loop_real_hops(
                idx,
                adj_dir,
                index_to_name,
                input_idx,
                max_real_hops=hop_limit,
                max_real_nodes=max_real_nodes,
                filter_regulatory=filter_regulatory,
                active_tf_names=active_tf_names,
                seed_symbols=seed_symbols,
                seed_gin_names=seed_gin_names,
            )
            hit_cap = hit_cap or capped
            if hop is not None and hop > 0:
                setattr(lr, f"{attr}_self_loop", True)
                setattr(lr, f"{attr}_self_loop_hops", hop)
                add_real_hop_path(walk, adj_dir, index_to_name, name_to_index, idx, idx, hop, parents)
                walk.pair_min_real_hops[(node, node)] = hop

        loop_records.append(lr)

    walk.hit_max_real_nodes = getattr(walk, "hit_max_real_nodes", False) or hit_cap

    n_rna = sum(1 for r in loop_records if r.rna_self_loop)
    n_prot = sum(1 for r in loop_records if r.protein_self_loop)
    n_prom = sum(1 for r in loop_records if r.promoter_self_loop)
    n_either = sum(
        1 for r in loop_records if r.rna_self_loop or r.protein_self_loop or r.promoter_self_loop
    )
    failed = [
        r.symbol
        for r, rec in zip(loop_records, seed_records)
        if rec.kind == "gene"
        and not (r.rna_self_loop or r.protein_self_loop or r.promoter_self_loop)
    ]

    qc = SeedSelfLoopQC(
        n_seeds=len(seed_records),
        n_triplet_complete=sum(1 for r in seed_records if r.triplet_complete),
        n_rna_self_loop=n_rna,
        n_protein_self_loop=n_prot,
        n_promoter_self_loop=n_prom,
        n_either_self_loop=n_either,
        n_failed=len(failed),
        max_real_hops=hop_limit,
        failed_symbols_sample=failed[:50],
        node_stats=node_count_stats(walk, n_seeds_no_loop=len(failed)),
    )
    return walk, qc

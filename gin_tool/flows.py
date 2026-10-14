"""Long cross-layer directed paths (information / regulatory axes) on a subgraph."""

from __future__ import annotations

import json
from typing import Callable, Sequence

from gin_tool.mapping import gin_entity_layer, is_virtual_gin_node


def _out_neighbors(adj_dir, u: int) -> list[int]:
    lo, hi = adj_dir.indptr[u], adj_dir.indptr[u + 1]
    return [int(v) for v in adj_dir.indices[lo:hi]]


def _path_layer_score(path: Sequence[int], name_fn: Callable[[int], str]) -> tuple[int, int]:
    """
    Score = (number of distinct **real** coarse layers, path length).
    Lexicographic maximize: prefer cross-omics diversity, then longer axes.
    """
    layers: set[str] = set()
    for idx in path:
        nm = name_fn(idx)
        if is_virtual_gin_node(nm):
            continue
        layers.add(gin_entity_layer(nm))
    return (len(layers), len(path))


def enumerate_simple_paths_from_starts(
    adj_dir,
    start_nodes: Sequence[int],
    *,
    max_depth: int = 24,
    max_paths_per_start: int = 400,
    max_total_paths: int = 5000,
) -> list[tuple[int, ...]]:
    """
    Depth-first enumeration of **simple** directed paths from each start.

    Paths are capped by ``max_depth`` (vertex count including the start).
    """
    all_paths: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()

    def record(p: tuple[int, ...]) -> None:
        if p not in seen:
            seen.add(p)
            all_paths.append(p)

    for s in start_nodes:
        if len(all_paths) >= max_total_paths:
            break
        col: list[tuple[int, ...]] = []

        def dfs(path: list[int], visited: set[int]) -> None:
            if len(col) >= max_paths_per_start or len(all_paths) + len(col) >= max_total_paths:
                return
            if len(path) >= max_depth:
                if len(path) >= 2:
                    col.append(tuple(path))
                return
            u = path[-1]
            nbrs = [v for v in _out_neighbors(adj_dir, u) if v not in visited]
            if not nbrs:
                if len(path) >= 2:
                    col.append(tuple(path))
                return
            for v in nbrs:
                if len(col) >= max_paths_per_start:
                    break
                visited.add(v)
                path.append(v)
                dfs(path, visited)
                path.pop()
                visited.remove(v)

        dfs([s], {s})
        for p in col:
            record(p)
            if len(all_paths) >= max_total_paths:
                break

    return all_paths


def rank_cross_layer_paths(
    paths: Sequence[tuple[int, ...]],
    name_fn: Callable[[int], str],
    *,
    top_k: int = 40,
) -> list[tuple[int, ...]]:
    scored = [(_path_layer_score(p, name_fn), p) for p in paths]
    scored.sort(key=lambda x: x[0], reverse=True)
    out: list[tuple[int, ...]] = []
    for _sc, p in scored:
        if len(out) >= top_k:
            break
        out.append(p)
    return out


def jaccard(a: set[int], b: set[int]) -> float:
    if not a and not b:
        return 1.0
    u = len(a & b)
    den = len(a | b)
    return u / den if den else 0.0


def flow_pairwise_overlap(paths: Sequence[Sequence[int]]) -> list[list[float]]:
    """Jaccard overlap on vertex sets; parallel flows ≈ low Jaccard."""
    sets = [set(p) for p in paths]
    m = len(sets)
    mat = [[0.0] * m for _ in range(m)]
    for i in range(m):
        mat[i][i] = 1.0
        for j in range(i + 1, m):
            v = jaccard(sets[i], sets[j])
            mat[i][j] = mat[j][i] = v
    return mat


def classify_parallel_vs_intersect(
    paths: Sequence[Sequence[int]],
    *,
    parallel_threshold: float = 0.12,
    intersect_threshold: float = 0.35,
) -> list[dict[str, object]]:
    """Per unordered pair of paths: label ``parallel`` / ``intersect`` / ``mixed``."""
    mat = flow_pairwise_overlap(paths)
    m = len(paths)
    rows: list[dict[str, object]] = []
    for i in range(m):
        for j in range(i + 1, m):
            v = mat[i][j]
            if v <= parallel_threshold:
                label = "parallel"
            elif v >= intersect_threshold:
                label = "intersect"
            else:
                label = "mixed"
            rows.append({"path_i": i, "path_j": j, "jaccard": v, "label": label})
    return rows


def write_flows_json(
    path: str,
    paths: Sequence[Sequence[int]],
    name_fn: Callable[[int], str],
) -> None:
    out = []
    for pi, p in enumerate(paths):
        layers = sorted(
            {gin_entity_layer(name_fn(i)) for i in p if not is_virtual_gin_node(name_fn(i))}
        )
        out.append(
            {
                "id": pi,
                "length": len(p),
                "distinct_real_layers": layers,
                "nodes": [name_fn(i) for i in p],
            }
        )
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)


def write_flow_overlap_json(path: str, paths: Sequence[Sequence[int]]) -> None:
    detail = classify_parallel_vs_intersect(paths)
    mat = flow_pairwise_overlap(paths)
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"pairwise_jaccard": mat, "pairs": detail}, f, indent=2)

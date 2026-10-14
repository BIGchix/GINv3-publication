"""k-hop subgraph extraction around seeds without materializing full GIN projection."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from scipy import sparse

from gin_tool.io import (
    IndexedEdges,
    build_indexed_edges,
    directed_csr_from_indexed,
    induced_subgraph_csr,
    iter_sifi_edges,
    khop_nodes_csr,
    khop_nodes_per_seed,
    symmetrize_undirected_unweighted,
)
from gin_tool.node_lists import k_for_seed


def seeds_from_names(
    name_to_index: dict[str, int],
    seed_names: Sequence[str],
) -> list[int]:
    out: list[int] = []
    for s in seed_names:
        if s not in name_to_index:
            raise KeyError(f"Unknown GIN node: {s!r}")
        out.append(name_to_index[s])
    return out


def extract_khop_subgraph(
    sifi_path: str,
    seed_gin_names: Sequence[str],
    k: int,
    *,
    max_lines: int | None = None,
) -> tuple[sparse.csr_matrix, dict[int, int], dict[int, int], IndexedEdges]:
    """
    Load edges (optionally first ``max_lines`` for large files), build directed CSR,
    symmetrize for topology, take k-hop closure around seeds, return induced subgraph.

    Does **not** compute the full ``R @ I.T`` entity–entity projection.
    """
    sub_sym, _sub_dir, old_to_new, new_to_old, ed = extract_khop_region(
        sifi_path, seed_gin_names, k, max_lines=max_lines
    )
    return sub_sym, old_to_new, new_to_old, ed


def extract_khop_region(
    sifi_path: str,
    seed_gin_names: Sequence[str],
    k: int,
    *,
    max_lines: int | None = None,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, dict[int, int], dict[int, int], IndexedEdges]:
    """
    Same k-hop node set as ``extract_khop_subgraph``, but return **both**:

    - ``sub_sym`` — undirected unweighted support (neighborhood / modules / distances).
    - ``sub_dir`` — directed signed subgraph on the same vertices (information flows).

    Vertex relabeling ``old_to_new`` / ``new_to_old`` is identical for both matrices.
    """
    ed = build_indexed_edges(sifi_path, max_lines=max_lines)
    seeds = seeds_from_names(ed.name_to_index, seed_gin_names)
    adj_dir = directed_csr_from_indexed(ed)
    sym = symmetrize_undirected_unweighted(adj_dir)
    nodes = khop_nodes_csr(sym, seeds, k)
    sub_sym, old_to_new, new_to_old = induced_subgraph_csr(sym, nodes)
    sub_dir, _, _ = induced_subgraph_csr(adj_dir, nodes)
    return sub_sym, sub_dir, old_to_new, new_to_old, ed


def extract_khop_region_mixed(
    sifi_path: str,
    seed_gin_names: Sequence[str],
    layer_map: dict[str, str],
    *,
    max_lines: int | None = None,
    user_k: int | None = None,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, dict[int, int], dict[int, int], IndexedEdges]:
    """
    k-hop closure with **per-seed** ``k`` (from nodeList layer), then union + induce.

    Protein/TF/HME seeds use k=1; other listed layers use k=2 unless ``user_k`` overrides all.
    """
    ed = build_indexed_edges(sifi_path, max_lines=max_lines)
    seeds = seeds_from_names(ed.name_to_index, seed_gin_names)
    k_list = [k_for_seed(n, layer_map, user_k=user_k) for n in seed_gin_names]
    adj_dir = directed_csr_from_indexed(ed)
    sym = symmetrize_undirected_unweighted(adj_dir)
    nodes = khop_nodes_per_seed(sym, seeds, k_list)
    sub_sym, old_to_new, new_to_old = induced_subgraph_csr(sym, nodes)
    sub_dir, _, _ = induced_subgraph_csr(adj_dir, nodes)
    return sub_sym, sub_dir, old_to_new, new_to_old, ed


def extract_khop_subgraph_mixed(
    sifi_path: str,
    seed_gin_names: Sequence[str],
    layer_map: dict[str, str],
    *,
    max_lines: int | None = None,
    user_k: int | None = None,
) -> tuple[sparse.csr_matrix, dict[int, int], dict[int, int], IndexedEdges]:
    """Sym subgraph only; per-seed k via ``layer_map``."""
    sub_sym, _sub_dir, old_to_new, new_to_old, ed = extract_khop_region_mixed(
        sifi_path, seed_gin_names, layer_map, max_lines=max_lines, user_k=user_k
    )
    return sub_sym, old_to_new, new_to_old, ed


def subgraph_edges_as_tuples(
    sub_adj: sparse.csr_matrix,
    new_to_old: dict[int, int],
    index_to_name: Sequence[str],
) -> list[tuple[str, str, float]]:
    """Enumerate directed edges of induced subgraph with GIN string names."""
    rows, cols = sub_adj.nonzero()
    out: list[tuple[str, str, float]] = []
    for i, j in zip(rows, cols):
        oi, oj = new_to_old[int(i)], new_to_old[int(j)]
        out.append((index_to_name[oi], index_to_name[oj], float(sub_adj[int(i), int(j)])))
    return out


def write_subgraph_tsv(
    path: str,
    edges: Iterable[tuple[str, str, float]],
) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for a, b, w in edges:
            f.write(f"{a}\t{b}\t{w:g}\n")


def subgraph_from_named_edges(
    edges: Sequence[tuple[str, str, float]],
    ed: IndexedEdges,
) -> tuple[sparse.csr_matrix, dict[int, int], dict[int, int]]:
    """Build induced subgraph CSR + index maps from named edge list."""
    node_names = sorted({a for a, _, _ in edges} | {b for _, b, _ in edges})
    name_to_local = {n: i for i, n in enumerate(node_names)}
    old_to_new = {
        ed.name_to_index[n]: name_to_local[n]
        for n in node_names
        if n in ed.name_to_index
    }
    new_to_old = {loc: ed.name_to_index[n] for n, loc in name_to_local.items()}
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for a, b, w in edges:
        if a not in name_to_local or b not in name_to_local:
            continue
        i, j = name_to_local[a], name_to_local[b]
        rows.append(i)
        cols.append(j)
        data.append(float(w))
    n = len(node_names)
    sub = sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    return sub, old_to_new, new_to_old


def remap_seed_indices_to_subgraph(
    old_to_new: dict[int, int],
    global_seed_indices: Sequence[int],
) -> list[int]:
    """Map full-graph node indices to induced subgraph labels (0..n-1)."""
    out: list[int] = []
    for g in global_seed_indices:
        if g in old_to_new:
            out.append(old_to_new[g])
    return out

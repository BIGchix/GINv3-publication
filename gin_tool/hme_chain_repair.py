"""Post-prune histone / HME chain repair for context subgraphs.

Product direction:

  subgraph PPI proteins  →  (real-hop path)  →  outside HME  →  histone marks

Optional pre-step (default on): for **gene pathway seeds** that already have
``GENE.promoter.active`` in the subgraph, pull same-gene histone marks from
the full GIN that feed ``GENE.promoter;state``. Non-seed ``promoter.active``
(and ``promoter;state``) are dropped so peripheral genes do not drive
transcription without chromatin; ``.RNA`` / ``.protein`` are kept (RNA remains
a seed-injection entry). Then histone marks discover which HME enzymes matter
(reverse collect). Connectivity is bridged **outward from subgraph proteins**
to those HME, and HME→histone writer edges are restored.

Eraser catalytic chains (extra, not GIN-topology-limited to in-subgraph hist):

  substrate_mark → {substrate;ENZYME} → product_mark
  ENZYME.protein → {substrate;ENZYME} → product_mark

GIN deliberately omits product→promoter;state for demethylated / deacetylated
intermediates (e.g. bare H3K27). Repair still pulls those product nodes and
eraser edges so HDAC/SIRT/KDM are not left as PPI-only stubs.

Post-repair cleanup (default on):

1. Drop histone marks with **no full-GIN in-edge** (source-only ``H3K9ac`` etc.)
2. Drop remaining histones that do not link to ``promoter;state`` / ``.active``
3. Drop undirected islands with no core keep-seed
4. Drop hist-catalysis **virtual stubs** (``MARK;ENZYME`` with no real histone neighbor)
5. **Reverse (ancestor) prune** from seed anchors so orphan HME enzymes
   (e.g. HDAC with only stub outs) are removed
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np
from scipy import sparse

from gin_tool.active_tf import gin_node_to_gene_symbol, load_tf_symbol_pool
from gin_tool.io import IndexedEdges, directed_csr_from_indexed
from gin_tool.mapping import is_virtual_gin_node
from gin_tool.node_lists import load_node_list_dir

# Heuristic eraser families in nodeList_HME (writers use other symbols).
_ERASER_TOKENS = ("HDAC", "SIRT", "KDM", "JMJD", "PHF2", "PHF8")


def _outs(adj: sparse.csr_matrix, i: int) -> np.ndarray:
    return adj.indices[adj.indptr[i] : adj.indptr[i + 1]]


def _out_edges(adj: sparse.csr_matrix, i: int) -> Iterator[tuple[int, float]]:
    lo, hi = adj.indptr[i], adj.indptr[i + 1]
    for j, w in zip(adj.indices[lo:hi], adj.data[lo:hi]):
        yield int(j), float(w)


def _is_histone_mark(name: str) -> bool:
    u = name.upper()
    return any(h in u for h in ("H3K", "H4K", "H2AK", "H2BK", "H3AC", "H4AC"))


def _is_protein(name: str) -> bool:
    return name.lower().endswith(".protein") and ";" not in name


def _gene_from_promoter_active(name: str) -> str | None:
    if not name.endswith(".promoter.active"):
        return None
    return name[: -len(".promoter.active")]


def _gene_from_promoter_state(name: str) -> str | None:
    low = name.lower()
    if not low.endswith(".promoter;state") and ".promoter;state" not in low:
        return None
    # Canonical: GENE.promoter;state
    if name.endswith(".promoter;state"):
        return name[: -len(".promoter;state")]
    return None


def is_gene_seed_symbol(symbol: str) -> bool:
    """True for gene-like seed symbols (exclude SID/CID/CHEBI metabolite tokens)."""
    s = (symbol or "").strip()
    if not s or is_virtual_gin_node(s):
        return False
    u = s.upper()
    if u.startswith(("SID", "CID", "CHEBI", "PUBCHEM")):
        return False
    if s.endswith(".protein") or s.endswith(".RNA") or ".promoter" in s.lower():
        # Allow bare gene only; callers should pass symbols not places.
        return False
    return True


def gene_symbols_from_seed_inputs(
    symbols: Iterable[str],
    *,
    seed_kinds: dict[str, str] | None = None,
) -> set[str]:
    """Select gene seed symbols from union seed list / kinds map."""
    out: set[str] = set()
    kinds = seed_kinds or {}
    for raw in symbols:
        sym = (raw or "").strip()
        if not sym:
            continue
        kind = kinds.get(sym) or kinds.get(sym.upper())
        if kind == "metabolite":
            continue
        if kind == "gene" or (kind is None and is_gene_seed_symbol(sym)):
            # Normalize away accidental .protein suffix
            if sym.endswith(".protein"):
                sym = sym[: -len(".protein")]
            out.add(sym)
    return out


def drop_nonseed_promoter_active(
    edge_set: dict[tuple[str, str], float],
    nodes: set[str],
    *,
    seed_genes: set[str],
) -> dict[str, object]:
    """Remove non-seed ``promoter.active`` / ``promoter;state``; keep RNA/protein.

    Does not remove ``.RNA`` or ``.protein``. Hist marks left dangling are
    cleaned by subsequent histprune / reverse-ancestor prune.
    """
    drop_nodes: set[str] = set()
    for n in list(nodes):
        if is_virtual_gin_node(n):
            continue
        g_act = _gene_from_promoter_active(n)
        if g_act is not None and g_act not in seed_genes:
            drop_nodes.add(n)
            continue
        g_st = _gene_from_promoter_state(n)
        if g_st is not None and g_st not in seed_genes:
            drop_nodes.add(n)
    if not drop_nodes:
        return {
            "n_nonseed_active_dropped": 0,
            "n_nonseed_state_dropped": 0,
            "n_edges_removed": 0,
            "dropped_active_sample": [],
        }
    n_active = sum(1 for n in drop_nodes if n.endswith(".promoter.active"))
    n_state = sum(1 for n in drop_nodes if "promoter;state" in n.lower())
    rm_edges = [e for e in edge_set if e[0] in drop_nodes or e[1] in drop_nodes]
    for e in rm_edges:
        del edge_set[e]
    nodes -= drop_nodes
    sample = sorted(n for n in drop_nodes if n.endswith(".promoter.active"))[:40]
    return {
        "n_nonseed_active_dropped": n_active,
        "n_nonseed_state_dropped": n_state,
        "n_edges_removed": len(rm_edges),
        "dropped_active_sample": sample,
    }


def backfill_active_promoter_hist_marks(
    edge_set: dict[tuple[str, str], float],
    nodes: set[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
    adj_rev: sparse.csr_matrix,
    seed_genes: set[str] | None = None,
) -> dict[str, object]:
    """
    For each ``GENE.promoter.active`` in the subgraph (optionally restricted to
    ``seed_genes``), materialize same-gene histone marks from the full GIN that
    have a directed edge into ``GENE.promoter;state`` (or into ``.promoter.active``),
    plus those edges.

    Also restores missing ``promoter;state → promoter.active`` when the state
    node exists in GIN. Does **not** invent hist→state links absent from GIN.
    """
    stats: dict[str, object] = {
        "n_active_genes_scanned": 0,
        "n_active_genes_with_backfill": 0,
        "n_hist_marks_backfilled": 0,
        "n_hist_state_edges_backfilled": 0,
        "n_state_nodes_added": 0,
        "hist_backfill_sample": [],
        "backfill_seed_genes_only": seed_genes is not None,
        "n_seed_genes_filter": len(seed_genes) if seed_genes is not None else 0,
    }
    actives = sorted(
        n
        for n in nodes
        if n.endswith(".promoter.active") and not is_virtual_gin_node(n)
    )
    if seed_genes is not None:
        actives = [
            n
            for n in actives
            if (_gene_from_promoter_active(n) or "") in seed_genes
        ]
    stats["n_active_genes_scanned"] = len(actives)
    before_hist = {
        n
        for n in nodes
        if _is_histone_mark(n) and not is_virtual_gin_node(n)
    }
    genes_with_new: set[str] = set()
    new_marks: list[str] = []

    for pa in actives:
        gene = _gene_from_promoter_active(pa)
        if not gene:
            continue
        state = f"{gene}.promoter;state"
        prefix = f"{gene}.promoter."

        if state in name_to_index and state not in nodes:
            nodes.add(state)
            stats["n_state_nodes_added"] = int(stats["n_state_nodes_added"]) + 1
            _set_edge_from_gin(
                edge_set, state, pa, name_to_index=name_to_index, adj=adj
            )
        elif state in nodes and (state, pa) not in edge_set:
            _set_edge_from_gin(
                edge_set, state, pa, name_to_index=name_to_index, adj=adj
            )

        anchors = [a for a in (state, pa) if a in name_to_index]
        if not anchors:
            continue

        for anch in anchors:
            i = name_to_index[anch]
            for pred_i in _outs(adj_rev, i):
                pred = index_to_name[int(pred_i)]
                if is_virtual_gin_node(pred):
                    continue
                if not _is_histone_mark(pred):
                    continue
                if not pred.startswith(prefix):
                    continue
                if pred not in nodes:
                    nodes.add(pred)
                    if pred not in before_hist:
                        new_marks.append(pred)
                        genes_with_new.add(gene)
                if _set_edge_from_gin(
                    edge_set, pred, anch, name_to_index=name_to_index, adj=adj
                ):
                    stats["n_hist_state_edges_backfilled"] = (
                        int(stats["n_hist_state_edges_backfilled"]) + 1
                    )
                    genes_with_new.add(gene)

    # Deduplicate while preserving order for sample.
    seen_m: set[str] = set()
    uniq_new: list[str] = []
    for m in new_marks:
        if m not in seen_m:
            seen_m.add(m)
            uniq_new.append(m)
    stats["n_hist_marks_backfilled"] = len(uniq_new)
    stats["n_active_genes_with_backfill"] = len(genes_with_new)
    stats["hist_backfill_sample"] = uniq_new[:40]
    return stats


def is_eraser_hme_symbol(symbol: str) -> bool:
    """True for deacetylase / demethylase-like HME gene symbols."""
    u = (symbol or "").upper()
    if u.endswith(".PROTEIN"):
        u = u[: -len(".PROTEIN")]
    return any(tok in u for tok in _ERASER_TOKENS)


def eraser_symbols_from_hme(hme_symbols: Iterable[str]) -> set[str]:
    return {s for s in hme_symbols if is_eraser_hme_symbol(s)}


def _enzyme_from_catalysis_virtual(virtual: str, eraser_symbols: set[str]) -> str | None:
    """Parse ``MARK;ENZYME.protein`` (or swapped) → ``ENZYME.protein`` if eraser."""
    if ";" not in virtual:
        return None
    for part in virtual.split(";"):
        sym = gin_node_to_gene_symbol(part) or part
        if sym.lower().endswith(".protein"):
            sym = sym[: -len(".protein")]
        if sym in eraser_symbols:
            return f"{sym}.protein"
    return None


@dataclass(frozen=True)
class EraserCatalyticStep:
    substrate: str
    enzyme: str
    virtual: str
    product: str
    w_sub_v: float
    w_v_prod: float
    w_enz_v: float


def collect_eraser_catalytic_steps(
    substrate_marks: Sequence[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
    eraser_symbols: set[str],
    max_rounds: int = 3,
) -> list[EraserCatalyticStep]:
    """
    Walk GIN eraser reactions starting from subgraph substrate marks.

    Each step is ``substrate → virtual → product`` with ``enzyme → virtual``.
    Newly produced marks become substrates for further rounds (demethylation
    ladder), without inventing product→promoter;state edges.
    """
    if not eraser_symbols:
        return []
    frontier = {m for m in substrate_marks if m in name_to_index}
    seen_keys: set[tuple[str, str, str, str]] = set()
    steps: list[EraserCatalyticStep] = []

    for _ in range(max(1, int(max_rounds))):
        if not frontier:
            break
        next_frontier: set[str] = set()
        for mark in sorted(frontier):
            i = name_to_index[mark]
            for v, w_sv in _out_edges(adj, i):
                vn = index_to_name[v]
                if not is_virtual_gin_node(vn):
                    continue
                enz = _enzyme_from_catalysis_virtual(vn, eraser_symbols)
                if enz is None or enz not in name_to_index:
                    continue
                # enzyme must feed this virtual in GIN
                enz_i = name_to_index[enz]
                enz_outs = {j: w for j, w in _out_edges(adj, enz_i)}
                if v not in enz_outs:
                    continue
                w_ev = float(enz_outs[v])
                for r, w_vp in _out_edges(adj, v):
                    rn = index_to_name[r]
                    if not _is_histone_mark(rn) or is_virtual_gin_node(rn):
                        continue
                    if rn == mark:
                        continue
                    key = (mark, enz, vn, rn)
                    if key in seen_keys:
                        continue
                    seen_keys.add(key)
                    steps.append(
                        EraserCatalyticStep(
                            substrate=mark,
                            enzyme=enz,
                            virtual=vn,
                            product=rn,
                            w_sub_v=float(w_sv),
                            w_v_prod=float(w_vp),
                            w_enz_v=w_ev,
                        )
                    )
                    next_frontier.add(rn)
        frontier = {m for m in next_frontier if m in name_to_index}
    return steps


def apply_eraser_catalytic_steps(
    steps: Sequence[EraserCatalyticStep],
    edge_set: dict[tuple[str, str], float],
    nodes: set[str],
) -> tuple[int, int, set[str], set[str]]:
    """
    Materialize eraser steps into ``edge_set``.

    Returns (n_steps, n_edges_added, enzymes, products).
    """
    n_edges_before = len(edge_set)
    enzymes: set[str] = set()
    products: set[str] = set()
    for st in steps:
        enzymes.add(st.enzyme)
        products.add(st.product)
        nodes.add(st.substrate)
        nodes.add(st.enzyme)
        nodes.add(st.virtual)
        nodes.add(st.product)
        edge_set.setdefault((st.substrate, st.virtual), st.w_sub_v)
        edge_set.setdefault((st.virtual, st.product), st.w_v_prod)
        edge_set.setdefault((st.enzyme, st.virtual), st.w_enz_v)
    return len(steps), len(edge_set) - n_edges_before, enzymes, products


def load_hme_symbols(node_list_dir: str | Path) -> set[str]:
    """Bare gene symbols for HME enzymes from nodeList_HME (+ layer pool)."""
    node_list_dir = Path(node_list_dir)
    symbols: set[str] = set()
    pools = load_node_list_dir(node_list_dir)
    for x in pools.get("hme", []) or []:
        sym = gin_node_to_gene_symbol(x) or x
        if sym.lower().endswith(".protein"):
            sym = sym[: -len(".protein")]
        symbols.add(sym)
    hme_file = node_list_dir / "nodeList_HME.txt"
    if hme_file.is_file():
        for line in hme_file.read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            sym = gin_node_to_gene_symbol(s) or s
            if sym.lower().endswith(".protein"):
                sym = sym[: -len(".protein")]
            symbols.add(sym)
    return symbols


def reverse_collect_upstream(
    start_names: list[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj_rev: sparse.csr_matrix,
    max_graph_hops: int,
    hme_symbols: set[str],
    tf_symbols: set[str] | None = None,
) -> dict:
    """BFS reverse from histone marks; collect HME/TF proteins reached."""
    tf_symbols = tf_symbols or set()
    starts = [name_to_index[n] for n in start_names if n in name_to_index]
    visited = {i: 0 for i in starts}
    q = deque(starts)
    hme_hits: dict[str, int] = {}
    tf_hits: dict[str, int] = {}
    while q:
        u = q.popleft()
        d = visited[u]
        if d >= max_graph_hops:
            continue
        for v in _outs(adj_rev, int(u)):
            v = int(v)
            if v in visited:
                continue
            visited[v] = d + 1
            q.append(v)
            name = index_to_name[v]
            if is_virtual_gin_node(name):
                continue
            sym = gin_node_to_gene_symbol(name) or name
            if _is_protein(name) or name == sym:
                if sym in hme_symbols:
                    hme_hits[name] = min(hme_hits.get(name, 10**9), d + 1)
                if sym in tf_symbols:
                    tf_hits[name] = min(tf_hits.get(name, 10**9), d + 1)
    return {
        "n_visited": len(visited),
        "hme_proteins": hme_hits,
        "tf_proteins": tf_hits,
    }


def real_hop_from_sources_to_targets(
    sources: set[str],
    targets: list[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
    max_real_hops: int,
) -> dict[str, dict]:
    """Multi-source forward real-hop BFS: sources → targets (via virtuals)."""
    src_idx = [name_to_index[s] for s in sources if s in name_to_index]
    target_idx = {
        name_to_index[t]: t for t in targets if t in name_to_index
    }
    if not src_idx or not target_idx:
        return {
            t: {"reachable": False, "real_hops": None, "from": None, "path": None}
            for t in targets
        }

    # depth counts real hops; graph walk expands through one virtual layer per step
    prev: dict[int, int | None] = {i: None for i in src_idx}
    depth: dict[int, int] = {i: 0 for i in src_idx}
    src_of: dict[int, int] = {i: i for i in src_idx}
    q = deque(src_idx)
    hit: dict[str, dict] = {
        t: {"reachable": False, "real_hops": None, "from": None, "path": None}
        for t in targets
    }
    remaining = set(target_idx)
    while q and remaining:
        u = q.popleft()
        d_real = depth[u]
        if d_real >= max_real_hops:
            continue
        for v in _outs(adj, int(u)):
            v = int(v)
            vn = index_to_name[v]
            if not is_virtual_gin_node(vn):
                continue
            if v in prev:
                continue
            prev[v] = u
            # virtual layer does not increment real hop
            depth[v] = d_real
            src_of[v] = src_of[u]
            for r in _outs(adj, v):
                r = int(r)
                rn = index_to_name[r]
                if is_virtual_gin_node(rn):
                    continue
                if r in prev:
                    continue
                prev[r] = v
                depth[r] = d_real + 1
                src_of[r] = src_of[u]
                if r in remaining:
                    tname = target_idx[r]
                    path_idx = [r]
                    cur: int | None = r
                    while cur is not None and prev[cur] is not None:
                        cur = prev[cur]
                        path_idx.append(cur)  # type: ignore[arg-type]
                    path_idx.reverse()
                    path_names = [index_to_name[i] for i in path_idx]
                    hit[tname] = {
                        "reachable": True,
                        "real_hops": depth[r],
                        "graph_hops": len(path_names) - 1,
                        "from": path_names[0],
                        "path": path_names,
                    }
                    remaining.discard(r)
                q.append(r)
    return hit


def undirected_shortest_from_sources_to_targets(
    sources: set[str],
    targets: list[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
    max_real_hops: int,
) -> dict[str, dict]:
    """Undirected real-hop BFS on the real–real projection; path ordered source→target."""
    adj_rev = adj.T.tocsr()

    def real_neighbors(i: int) -> list[int]:
        out: list[int] = []
        for v in list(_outs(adj, i)) + list(_outs(adj_rev, i)):
            v = int(v)
            if not is_virtual_gin_node(index_to_name[v]):
                continue
            for r in list(_outs(adj, v)) + list(_outs(adj_rev, v)):
                r = int(r)
                if r == i:
                    continue
                if is_virtual_gin_node(index_to_name[r]):
                    continue
                out.append(r)
        return out

    src_idx = [name_to_index[s] for s in sources if s in name_to_index]
    target_idx = {
        name_to_index[t]: t for t in targets if t in name_to_index
    }
    if not src_idx or not target_idx:
        return {
            t: {"reachable": False, "real_hops": None, "from": None, "path": None}
            for t in targets
        }

    start = list(src_idx)
    prev: dict[int, int | None] = {i: None for i in start}
    depth = {i: 0 for i in start}
    src_of = {i: i for i in start}
    q = deque(start)
    hit: dict[str, dict] = {
        t: {"reachable": False, "real_hops": None, "from": None, "path": None}
        for t in targets
    }
    remaining = set(target_idx)
    while q and remaining:
        u = q.popleft()
        if depth[u] >= max_real_hops:
            continue
        for r in real_neighbors(int(u)):
            if r in prev:
                continue
            prev[r] = u
            depth[r] = depth[u] + 1
            src_of[r] = src_of[u]
            if r in remaining:
                tname = target_idx[r]
                path_idx = [r]
                cur: int | None = r
                while cur is not None and prev[cur] is not None:
                    cur = prev[cur]
                    path_idx.append(cur)  # type: ignore[arg-type]
                path_idx.reverse()
                path_names = [index_to_name[i] for i in path_idx]
                hit[tname] = {
                    "reachable": True,
                    "real_hops": len(path_names) - 1,
                    "from": path_names[0],
                    "path": path_names,
                }
                remaining.discard(r)
            q.append(r)
    return hit


def _edge_weight(adj: sparse.csr_matrix, u: int, v: int) -> float:
    """Signed GIN weight; 0 if no arc."""
    return float(adj[u, v])


def _set_edge_from_gin(
    edge_set: dict[tuple[str, str], float],
    a: str,
    b: str,
    *,
    name_to_index: dict[str, int],
    adj: sparse.csr_matrix,
    default: float = 1.0,
) -> bool:
    """Insert (a,b) with GIN sign if missing. Returns True if newly added."""
    if (a, b) in edge_set:
        return False
    if a in name_to_index and b in name_to_index:
        w = _edge_weight(adj, name_to_index[a], name_to_index[b])
        if w == 0.0:
            w = float(default)
        edge_set[(a, b)] = float(w)
    else:
        edge_set[(a, b)] = float(default)
    return True


def _connect_real_pair_via_virtual(
    a: str,
    b: str,
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
    edge_set: dict[tuple[str, str], float],
    nodes: set[str],
) -> int:
    """Add directed a→V→b edges when a shared virtual exists. Returns edges added."""
    if a not in name_to_index or b not in name_to_index:
        return 0
    u = name_to_index[a]
    tidx = name_to_index[b]
    added = 0
    for v in _outs(adj, u):
        v = int(v)
        vn = index_to_name[v]
        if not is_virtual_gin_node(vn):
            continue
        if tidx in {int(x) for x in _outs(adj, v)}:
            if _set_edge_from_gin(
                edge_set, a, vn, name_to_index=name_to_index, adj=adj
            ):
                added += 1
            if _set_edge_from_gin(
                edge_set, vn, b, name_to_index=name_to_index, adj=adj
            ):
                added += 1
            nodes.add(vn)
    return added


def _hme_to_hist_edges(
    hme_nodes: list[str],
    hist_set: set[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
) -> list[tuple[str, str | None, str]]:
    paths: list[tuple[str, str | None, str]] = []
    for h in hme_nodes:
        if h not in name_to_index:
            continue
        i = name_to_index[h]
        for v in _outs(adj, i):
            vn = index_to_name[int(v)]
            if is_virtual_gin_node(vn):
                for r in _outs(adj, int(v)):
                    rn = index_to_name[int(r)]
                    if rn in hist_set:
                        paths.append((h, vn, rn))
            elif vn in hist_set:
                paths.append((h, None, vn))
    return paths


@dataclass
class HmeChainRepairQC:
    direction: str = "sub_PPI → outside_HME → histone"
    n_hist_marks: int = 0
    n_hist_marks_before_backfill: int = 0
    n_sub_proteins: int = 0
    n_hme_upstream: int = 0
    n_hme_in_sub: int = 0
    n_hme_outside: int = 0
    n_hme_bridged: int = 0
    n_hme_to_hist_paths: int = 0
    n_eraser_steps: int = 0
    n_eraser_enzymes: int = 0
    n_eraser_products: int = 0
    n_eraser_edges_added: int = 0
    eraser_product_suffix_counts: dict[str, int] = field(default_factory=dict)
    eraser_enzyme_sample: list[str] = field(default_factory=list)
    n_nodes_added: int = 0
    n_edges_added: int = 0
    max_real_hops: int = 3
    upstream_graph_hops: int = 4
    eraser_max_rounds: int = 3
    repair_eraser_chains: bool = True
    prune_sourceless_histones: bool = True
    backfill_active_hist_marks: bool = True
    backfill_seed_genes_only: bool = True
    drop_nonseed_promoter_active: bool = True
    n_seed_genes_filter: int = 0
    n_nonseed_active_dropped: int = 0
    n_nonseed_state_dropped: int = 0
    n_nonseed_active_edges_removed: int = 0
    dropped_active_sample: list[str] = field(default_factory=list)
    n_active_genes_scanned: int = 0
    n_active_genes_with_backfill: int = 0
    n_hist_marks_backfilled: int = 0
    n_hist_state_edges_backfilled: int = 0
    n_state_nodes_added: int = 0
    hist_backfill_sample: list[str] = field(default_factory=list)
    n_sourceless_hist_removed: int = 0
    sourceless_hist_suffix_counts: dict[str, int] = field(default_factory=dict)
    sourceless_hist_sample: list[str] = field(default_factory=list)
    n_dead_end_hist_removed: int = 0
    n_island_nodes_removed: int = 0
    n_island_edges_removed: int = 0
    island_node_sample: list[str] = field(default_factory=list)
    n_catalysis_stubs_removed: int = 0
    catalysis_stub_sample: list[str] = field(default_factory=list)
    n_reprune_nodes_removed: int = 0
    n_reprune_edges_removed: int = 0
    reprune_node_sample: list[str] = field(default_factory=list)
    directed_by_real_hops: dict[str, int] = field(default_factory=dict)
    undirected_by_real_hops: dict[str, int] = field(default_factory=dict)
    bridged_hme_sample: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "direction": self.direction,
            "n_hist_marks": self.n_hist_marks,
            "n_hist_marks_before_backfill": self.n_hist_marks_before_backfill,
            "backfill_active_hist_marks": self.backfill_active_hist_marks,
            "backfill_seed_genes_only": self.backfill_seed_genes_only,
            "drop_nonseed_promoter_active": self.drop_nonseed_promoter_active,
            "n_seed_genes_filter": self.n_seed_genes_filter,
            "n_nonseed_active_dropped": self.n_nonseed_active_dropped,
            "n_nonseed_state_dropped": self.n_nonseed_state_dropped,
            "n_nonseed_active_edges_removed": self.n_nonseed_active_edges_removed,
            "dropped_active_sample": self.dropped_active_sample,
            "n_active_genes_scanned": self.n_active_genes_scanned,
            "n_active_genes_with_backfill": self.n_active_genes_with_backfill,
            "n_hist_marks_backfilled": self.n_hist_marks_backfilled,
            "n_hist_state_edges_backfilled": self.n_hist_state_edges_backfilled,
            "n_state_nodes_added": self.n_state_nodes_added,
            "hist_backfill_sample": self.hist_backfill_sample,
            "n_sub_proteins": self.n_sub_proteins,
            "n_hme_upstream": self.n_hme_upstream,
            "n_hme_in_sub": self.n_hme_in_sub,
            "n_hme_outside": self.n_hme_outside,
            "n_hme_bridged": self.n_hme_bridged,
            "n_hme_to_hist_paths": self.n_hme_to_hist_paths,
            "repair_eraser_chains": self.repair_eraser_chains,
            "eraser_max_rounds": self.eraser_max_rounds,
            "n_eraser_steps": self.n_eraser_steps,
            "n_eraser_enzymes": self.n_eraser_enzymes,
            "n_eraser_products": self.n_eraser_products,
            "n_eraser_edges_added": self.n_eraser_edges_added,
            "eraser_product_suffix_counts": self.eraser_product_suffix_counts,
            "eraser_enzyme_sample": self.eraser_enzyme_sample,
            "prune_sourceless_histones": self.prune_sourceless_histones,
            "n_sourceless_hist_removed": self.n_sourceless_hist_removed,
            "sourceless_hist_suffix_counts": self.sourceless_hist_suffix_counts,
            "sourceless_hist_sample": self.sourceless_hist_sample,
            "n_dead_end_hist_removed": self.n_dead_end_hist_removed,
            "n_island_nodes_removed": self.n_island_nodes_removed,
            "n_island_edges_removed": self.n_island_edges_removed,
            "island_node_sample": self.island_node_sample,
            "n_catalysis_stubs_removed": self.n_catalysis_stubs_removed,
            "catalysis_stub_sample": self.catalysis_stub_sample,
            "n_reprune_nodes_removed": self.n_reprune_nodes_removed,
            "n_reprune_edges_removed": self.n_reprune_edges_removed,
            "reprune_node_sample": self.reprune_node_sample,
            "n_nodes_added": self.n_nodes_added,
            "n_edges_added": self.n_edges_added,
            "max_real_hops": self.max_real_hops,
            "upstream_graph_hops": self.upstream_graph_hops,
            "directed_by_real_hops": self.directed_by_real_hops,
            "undirected_by_real_hops": self.undirected_by_real_hops,
            "bridged_hme_sample": self.bridged_hme_sample,
        }


def _gin_in_degree(
    name: str,
    *,
    name_to_index: dict[str, int],
    adj_rev: sparse.csr_matrix,
) -> int:
    i = name_to_index.get(name)
    if i is None:
        return 0
    return int(adj_rev.indptr[i + 1] - adj_rev.indptr[i])


def _is_hme_protein_node(name: str, hme_symbols: set[str]) -> bool:
    if not _is_protein(name):
        return False
    sym = gin_node_to_gene_symbol(name) or name
    if sym.lower().endswith(".protein"):
        sym = sym[: -len(".protein")]
    return sym in hme_symbols or is_eraser_hme_symbol(sym)


def _histone_links_promoter(name: str, outs: dict[str, set[str]]) -> bool:
    for nxt in outs.get(name, ()):
        low = nxt.lower()
        if ";state" in low or low.endswith(".promoter.active"):
            return True
    return False


def _is_island_keep_seed(
    name: str,
    *,
    hme_symbols: set[str],
    outs: dict[str, set[str]],
) -> bool:
    """Nodes that may anchor a retained undirected component after hist prune."""
    if is_virtual_gin_node(name):
        return False
    if _is_histone_mark(name):
        # Dead-end eraser products (bare H3K*) must not keep a lone island.
        return _histone_links_promoter(name, outs)
    if _is_hme_protein_node(name, hme_symbols):
        return False
    return True


def _is_hist_catalysis_stub(
    name: str,
    *,
    inn: dict[str, set[str]],
    outs: dict[str, set[str]],
) -> bool:
    """Virtual ``MARK;ENZYME`` with no real histone place neighbor left."""
    if not is_virtual_gin_node(name):
        return False
    u = name.upper()
    if "H3K" not in u and "H4K" not in u and "H2AK" not in u and "H2BK" not in u:
        return False
    nbrs = inn.get(name, set()) | outs.get(name, set())
    return not any(
        _is_histone_mark(x) and not is_virtual_gin_node(x) for x in nbrs
    )


def _reverse_ancestor_keep(
    edge_set: dict[tuple[str, str], float],
    anchors: set[str],
) -> set[str]:
    """Nodes with a directed path *to* an anchor (ancestors ∪ anchors)."""
    nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    seeds = {a for a in anchors if a in nodes}
    if not seeds:
        return set(nodes)
    pred: dict[str, set[str]] = {n: set() for n in nodes}
    for a, b in edge_set:
        pred[b].add(a)
    keep = set(seeds)
    stack = list(seeds)
    while stack:
        v = stack.pop()
        for u in pred.get(v, ()):
            if u not in keep:
                keep.add(u)
                stack.append(u)
    return keep


def prune_sourceless_histones_and_islands(
    edge_set: dict[tuple[str, str], float],
    *,
    name_to_index: dict[str, int],
    adj_rev: sparse.csr_matrix,
    hme_symbols: set[str],
    anchor_nodes: set[str] | None = None,
) -> tuple[dict[tuple[str, str], float], dict[str, object]]:
    """
    After HME/eraser repair:

    1. Drop histone marks that have **no incoming edge in the full GIN**
       (e.g. many ``*.promoter.H3K9ac`` are GIN sources-only) and all incident
       subgraph edges / now-dangling virtuals.
    2. Drop undirected components that no longer touch a keep-seed
       (gene/RNA/promoter/metabolite/writer-backed histone) — typically orphan
       HME enzymes and bare H3K intermediates left behind by (1).
    3. Drop hist-catalysis virtual stubs (no real histone neighbor).
    4. Reverse-ancestor prune from ``anchor_nodes`` (drop HDAC etc. that no
       longer lie on any path into a seed).
    """
    stats: dict[str, object] = {
        "n_sourceless_hist_removed": 0,
        "sourceless_hist_suffix_counts": {},
        "sourceless_hist_sample": [],
        "n_dead_end_hist_removed": 0,
        "n_island_nodes_removed": 0,
        "n_island_edges_removed": 0,
        "island_node_sample": [],
        "n_catalysis_stubs_removed": 0,
        "catalysis_stub_sample": [],
        "n_reprune_nodes_removed": 0,
        "n_reprune_edges_removed": 0,
        "reprune_node_sample": [],
    }
    if not edge_set:
        return edge_set, stats

    nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    hist = [
        n
        for n in nodes
        if _is_histone_mark(n) and not is_virtual_gin_node(n)
    ]
    sourceless = sorted(
        h
        for h in hist
        if _gin_in_degree(h, name_to_index=name_to_index, adj_rev=adj_rev) <= 0
    )
    stats["n_sourceless_hist_removed"] = len(sourceless)
    stats["sourceless_hist_sample"] = sourceless[:40]
    suf: Counter[str] = Counter(h.split(".")[-1] for h in sourceless)
    stats["sourceless_hist_suffix_counts"] = dict(suf.most_common())

    drop: set[str] = set(sourceless)

    def _strip_nodes(to_drop: set[str]) -> None:
        nonlocal edge_set
        if not to_drop:
            return
        edge_set = {
            (a, b): w
            for (a, b), w in edge_set.items()
            if a not in to_drop and b not in to_drop
        }
        # Drop now-dangling virtuals (degree 0).
        for _ in range(4):
            cur = {a for a, _ in edge_set} | {b for _, b in edge_set}
            deg: Counter[str] = Counter()
            for a, b in edge_set:
                deg[a] += 1
                deg[b] += 1
            doomed = {
                n for n in cur if is_virtual_gin_node(n) and deg.get(n, 0) == 0
            }
            if not doomed:
                break
            edge_set = {
                (a, b): w
                for (a, b), w in edge_set.items()
                if a not in doomed and b not in doomed
            }

    _strip_nodes(drop)

    # Eraser products (bare H3K*) have no →promoter;state in GIN; after their
    # sourceless substrates are gone they only feed token piles — drop them too.
    nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    outs_tmp: dict[str, set[str]] = {n: set() for n in nodes}
    for a, b in edge_set:
        outs_tmp[a].add(b)
    dead_end_hist = sorted(
        n
        for n in nodes
        if _is_histone_mark(n)
        and not is_virtual_gin_node(n)
        and not _histone_links_promoter(n, outs_tmp)
    )
    if dead_end_hist:
        _strip_nodes(set(dead_end_hist))
        stats["n_dead_end_hist_removed"] = len(dead_end_hist)

    # Island prune: keep undirected components that contain a keep-seed.
    nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    und: dict[str, set[str]] = {n: set() for n in nodes}
    outs: dict[str, set[str]] = {n: set() for n in nodes}
    for a, b in edge_set:
        und[a].add(b)
        und[b].add(a)
        outs[a].add(b)

    keep_seeds = {
        n
        for n in nodes
        if _is_island_keep_seed(n, hme_symbols=hme_symbols, outs=outs)
    }
    # Also treat degree-0 keep-seeds (isolated RNA etc.) as retained singletons.
    retained: set[str] = set()
    seen: set[str] = set()
    for seed in keep_seeds:
        if seed in seen:
            continue
        stack = [seed]
        comp: list[str] = []
        seen.add(seed)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in und.get(u, ()):
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        retained.update(comp)

    # Isolated nodes with no edges: keep only if keep-seed.
    for n in nodes:
        if n not in seen and n in keep_seeds:
            retained.add(n)

    island_nodes = sorted(nodes - retained)
    if island_nodes:
        before_e = len(edge_set)
        edge_set = {
            (a, b): w
            for (a, b), w in edge_set.items()
            if a in retained and b in retained
        }
        stats["n_island_nodes_removed"] = len(island_nodes)
        stats["n_island_edges_removed"] = before_e - len(edge_set)
        stats["island_node_sample"] = island_nodes[:40]
    else:
        # Still drop explicit degree-0 non-keep nodes left after hist removal.
        nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
        deg = Counter()
        for a, b in edge_set:
            deg[a] += 1
            deg[b] += 1
        orphan = sorted(
            n
            for n in nodes
            if deg.get(n, 0) == 0
            and not _is_island_keep_seed(n, hme_symbols=hme_symbols, outs=outs)
        )
        if orphan:
            drop_o = set(orphan)
            edge_set = {
                (a, b): w
                for (a, b), w in edge_set.items()
                if a not in drop_o and b not in drop_o
            }
            # degree-0 nodes aren't in edge_set keys as endpoints of edges —
            # they disappear once we rebuild from edge_set; track count.
            stats["n_island_nodes_removed"] = len(orphan)
            stats["island_node_sample"] = orphan[:40]

    # Catalysis stubs: MARK;ENZYME virtuals whose real histone ends were pruned.
    nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    inn2: dict[str, set[str]] = {n: set() for n in nodes}
    outs2: dict[str, set[str]] = {n: set() for n in nodes}
    for a, b in edge_set:
        outs2[a].add(b)
        inn2[b].add(a)
    stubs = sorted(
        n
        for n in nodes
        if _is_hist_catalysis_stub(n, inn=inn2, outs=outs2)
    )
    if stubs:
        _strip_nodes(set(stubs))
        stats["n_catalysis_stubs_removed"] = len(stubs)
        stats["catalysis_stub_sample"] = stubs[:40]

    # Reverse ancestor prune from seed anchors (drop orphan HDAC etc.).
    anchors = set(anchor_nodes or ())
    if anchors:
        before_nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
        before_e = len(edge_set)
        keep = _reverse_ancestor_keep(edge_set, anchors)
        edge_set = {
            (a, b): w
            for (a, b), w in edge_set.items()
            if a in keep and b in keep
        }
        after_nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
        dropped = sorted(before_nodes - after_nodes)
        stats["n_reprune_nodes_removed"] = len(dropped)
        stats["n_reprune_edges_removed"] = before_e - len(edge_set)
        stats["reprune_node_sample"] = dropped[:40]

    return edge_set, stats


def repair_hme_chains_on_edges(
    edges: Sequence[tuple[str, str, float]],
    ed: IndexedEdges,
    *,
    hme_symbols: set[str],
    tf_symbols: set[str] | None = None,
    max_real_hops: int = 3,
    upstream_graph_hops: int = 4,
    prefer_directed: bool = True,
    repair_eraser_chains: bool = True,
    eraser_max_rounds: int = 3,
    prune_sourceless_histones: bool = True,
    backfill_active_hist_marks: bool = True,
    backfill_seed_genes: set[str] | None = None,
    do_drop_nonseed_promoter_active: bool = True,
    anchor_nodes: set[str] | None = None,
    n_jobs: int = 1,
) -> tuple[list[tuple[str, str, float]], HmeChainRepairQC]:
    """
    Augment subgraph edges with PPI→HME bridges, writer HME→hist edges, and
    eraser catalytic chains (substrate→product intermediates).

    Pre-step (default on): for gene pathway seeds with ``promoter.active``,
    backfill same-gene histone marks from GIN that feed ``promoter;state``,
    then attach enzymes. Non-seed ``promoter.active`` / ``;state`` are dropped
    by default (RNA/protein kept).

    Bridge paths prefer directed real-hop; unreachable targets fall back to
    undirected shortest real-hop (same as the fork experiment).

    Post-repair cleanup (default on): drop GIN-sourceless histones, dead-end
    intermediates, islands, catalysis stubs, then reverse-ancestor prune from
    ``anchor_nodes``.

    ``n_jobs``: when ``>1``, run directed and undirected hop BFS in parallel.
    """
    seed_filter = set(backfill_seed_genes) if backfill_seed_genes is not None else None
    # Default policy: seed-restricted when a seed gene set is supplied.
    seed_only = seed_filter is not None
    qc = HmeChainRepairQC(
        max_real_hops=max_real_hops,
        upstream_graph_hops=upstream_graph_hops,
        repair_eraser_chains=repair_eraser_chains,
        eraser_max_rounds=eraser_max_rounds,
        prune_sourceless_histones=prune_sourceless_histones,
        backfill_active_hist_marks=backfill_active_hist_marks,
        backfill_seed_genes_only=seed_only,
        drop_nonseed_promoter_active=bool(do_drop_nonseed_promoter_active and seed_only),
        n_seed_genes_filter=len(seed_filter) if seed_filter is not None else 0,
    )
    edge_set = {(a, b): float(w) for a, b, w in edges}
    n_edges_before = len(edge_set)
    nodes = {a for a, _, _ in edges} | {b for _, b, _ in edges}
    before_nodes = set(nodes)

    adj = directed_csr_from_indexed(ed)
    adj_rev = adj.T.tocsr()
    name_to = ed.name_to_index
    idx_to = ed.index_to_name

    hist_before = sorted(
        n for n in nodes if _is_histone_mark(n) and not is_virtual_gin_node(n)
    )
    qc.n_hist_marks_before_backfill = len(hist_before)

    if seed_only and do_drop_nonseed_promoter_active:
        drop_st = drop_nonseed_promoter_active(
            edge_set, nodes, seed_genes=seed_filter or set()
        )
        qc.n_nonseed_active_dropped = int(drop_st["n_nonseed_active_dropped"])
        qc.n_nonseed_state_dropped = int(drop_st["n_nonseed_state_dropped"])
        qc.n_nonseed_active_edges_removed = int(drop_st["n_edges_removed"])
        qc.dropped_active_sample = list(drop_st.get("dropped_active_sample") or [])

    if backfill_active_hist_marks:
        bf = backfill_active_promoter_hist_marks(
            edge_set,
            nodes,
            name_to_index=name_to,
            index_to_name=idx_to,
            adj=adj,
            adj_rev=adj_rev,
            seed_genes=seed_filter,
        )
        qc.n_active_genes_scanned = int(bf["n_active_genes_scanned"])
        qc.n_active_genes_with_backfill = int(bf["n_active_genes_with_backfill"])
        qc.n_hist_marks_backfilled = int(bf["n_hist_marks_backfilled"])
        qc.n_hist_state_edges_backfilled = int(bf["n_hist_state_edges_backfilled"])
        qc.n_state_nodes_added = int(bf["n_state_nodes_added"])
        qc.hist_backfill_sample = list(bf["hist_backfill_sample"] or [])

    hist_marks = sorted(
        n for n in nodes if _is_histone_mark(n) and not is_virtual_gin_node(n)
    )
    sub_proteins = {n for n in nodes if _is_protein(n)}
    qc.n_hist_marks = len(hist_marks)
    qc.n_sub_proteins = len(sub_proteins)
    if not hist_marks or not sub_proteins:
        after_nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
        qc.n_nodes_added = len(after_nodes - before_nodes)
        qc.n_edges_added = len(edge_set) - n_edges_before
        out = [(a, b, w) for (a, b), w in sorted(edge_set.items())]
        return out, qc

    up = reverse_collect_upstream(
        hist_marks,
        name_to_index=name_to,
        index_to_name=idx_to,
        adj_rev=adj_rev,
        max_graph_hops=upstream_graph_hops,
        hme_symbols=hme_symbols,
        tf_symbols=tf_symbols,
    )
    hme_prots = sorted(up["hme_proteins"])
    qc.n_hme_upstream = len(hme_prots)

    eraser_syms = eraser_symbols_from_hme(hme_symbols)
    eraser_steps: list[EraserCatalyticStep] = []
    if repair_eraser_chains and eraser_syms:
        eraser_steps = collect_eraser_catalytic_steps(
            hist_marks,
            name_to_index=name_to,
            index_to_name=idx_to,
            adj=adj,
            eraser_symbols=eraser_syms,
            max_rounds=eraser_max_rounds,
        )

    # Eraser enzymes required by catalytic steps (may extend reverse-collect set)
    eraser_enzymes = {st.enzyme for st in eraser_steps}
    hme_union = sorted(set(hme_prots) | eraser_enzymes)
    hme_in_sub = [h for h in hme_union if h in nodes]
    hme_out = [h for h in hme_union if h not in nodes]
    qc.n_hme_in_sub = len(hme_in_sub)
    qc.n_hme_outside = len(hme_out)

    expanded: list[str] = []
    for h in hme_out:
        if h in name_to:
            expanded.append(h)
        if not h.endswith(".protein"):
            cand = f"{h}.protein"
            if cand in name_to:
                expanded.append(cand)
    expanded = sorted({n for n in expanded if n in name_to})

    from gin_tool.subgraph_pipeline.parallel_util import run_parallel_pair

    hop_dist, sp_dist = run_parallel_pair(
        lambda: real_hop_from_sources_to_targets(
            sub_proteins,
            expanded,
            name_to_index=name_to,
            index_to_name=idx_to,
            adj=adj,
            max_real_hops=max_real_hops,
        ),
        lambda: undirected_shortest_from_sources_to_targets(
            sub_proteins,
            expanded,
            name_to_index=name_to,
            index_to_name=idx_to,
            adj=adj,
            max_real_hops=max_real_hops,
        ),
        n_jobs=n_jobs,
    )

    def _hop_hist(dist_map: dict[str, dict]) -> dict[str, int]:
        c: Counter[str] = Counter()
        for rec in dist_map.values():
            if rec.get("reachable"):
                c[str(int(rec["real_hops"]))] += 1
            else:
                c["unreachable"] += 1
        return dict(c)

    qc.directed_by_real_hops = _hop_hist(hop_dist)
    qc.undirected_by_real_hops = _hop_hist(sp_dist)

    bridge_dist = hop_dist if prefer_directed else sp_dist
    bridged: list[str] = []
    for hme, rec in bridge_dist.items():
        path = rec.get("path") or []
        if len(path) < 2:
            # directed miss → undirected fallback
            urec = sp_dist.get(hme) or {}
            path = urec.get("path") or []
            if len(path) < 2:
                continue
        bridged.append(hme)
        for n in path:
            nodes.add(n)
        for i in range(len(path) - 1):
            a, b = path[i], path[i + 1]
            if is_virtual_gin_node(a) or is_virtual_gin_node(b):
                _set_edge_from_gin(
                    edge_set, a, b, name_to_index=name_to, adj=adj
                )
                continue
            # real→real: materialize via shared virtual in GIN
            _connect_real_pair_via_virtual(
                a,
                b,
                name_to_index=name_to,
                index_to_name=idx_to,
                adj=adj,
                edge_set=edge_set,
                nodes=nodes,
            )
            # undirected path may require reverse orientation
            _connect_real_pair_via_virtual(
                b,
                a,
                name_to_index=name_to,
                index_to_name=idx_to,
                adj=adj,
                edge_set=edge_set,
                nodes=nodes,
            )

    # Always attach writer edges HME → histone (including HME already in-sub)
    all_hme_for_hist = sorted(set(hme_union) | set(bridged) | set(hme_in_sub))
    # normalize to nodes present in GIN
    hme_norm: list[str] = []
    for h in all_hme_for_hist:
        if h in name_to:
            hme_norm.append(h)
        cand = h if h.endswith(".protein") else f"{h}.protein"
        if cand in name_to:
            hme_norm.append(cand)
    hme_norm = sorted(set(hme_norm))
    # Include all histone marks currently in the subgraph (after backfill).
    hist_set = {
        n for n in nodes if _is_histone_mark(n) and not is_virtual_gin_node(n)
    }
    hme_hist = _hme_to_hist_edges(
        hme_norm,
        hist_set,
        name_to_index=name_to,
        index_to_name=idx_to,
        adj=adj,
    )
    qc.n_hme_to_hist_paths = len(hme_hist)
    for a, v, c in hme_hist:
        nodes.add(a)
        nodes.add(c)
        if v:
            nodes.add(v)
            _set_edge_from_gin(edge_set, a, v, name_to_index=name_to, adj=adj)
            _set_edge_from_gin(edge_set, v, c, name_to_index=name_to, adj=adj)
        else:
            _set_edge_from_gin(edge_set, a, c, name_to_index=name_to, adj=adj)

    # Eraser chains: pull demethylated/deacetylated products + enzyme catalysis
    # even when products are absent from the context hist set (no →state in GIN).
    if eraser_steps:
        n_steps, n_e_add, enz_set, prod_set = apply_eraser_catalytic_steps(
            eraser_steps, edge_set, nodes
        )
        qc.n_eraser_steps = n_steps
        qc.n_eraser_enzymes = len(enz_set)
        qc.n_eraser_products = len(prod_set)
        qc.n_eraser_edges_added = n_e_add
        qc.eraser_enzyme_sample = sorted(enz_set)[:40]
        suf_c: Counter[str] = Counter()
        for p in prod_set:
            suf_c[p.split(".")[-1]] += 1
        qc.eraser_product_suffix_counts = dict(suf_c.most_common())

    if prune_sourceless_histones:
        edge_set, prune_stats = prune_sourceless_histones_and_islands(
            edge_set,
            name_to_index=name_to,
            adj_rev=adj_rev,
            hme_symbols=hme_symbols,
            anchor_nodes=anchor_nodes,
        )
        qc.n_sourceless_hist_removed = int(prune_stats["n_sourceless_hist_removed"])
        qc.sourceless_hist_suffix_counts = dict(
            prune_stats["sourceless_hist_suffix_counts"] or {}
        )
        qc.sourceless_hist_sample = list(prune_stats["sourceless_hist_sample"] or [])
        qc.n_dead_end_hist_removed = int(prune_stats["n_dead_end_hist_removed"])
        qc.n_island_nodes_removed = int(prune_stats["n_island_nodes_removed"])
        qc.n_island_edges_removed = int(prune_stats["n_island_edges_removed"])
        qc.island_node_sample = list(prune_stats["island_node_sample"] or [])
        qc.n_catalysis_stubs_removed = int(prune_stats["n_catalysis_stubs_removed"])
        qc.catalysis_stub_sample = list(prune_stats["catalysis_stub_sample"] or [])
        qc.n_reprune_nodes_removed = int(prune_stats["n_reprune_nodes_removed"])
        qc.n_reprune_edges_removed = int(prune_stats["n_reprune_edges_removed"])
        qc.reprune_node_sample = list(prune_stats["reprune_node_sample"] or [])

    after_nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    qc.n_nodes_added = len(after_nodes - before_nodes)
    qc.n_edges_added = len(edge_set) - n_edges_before
    qc.n_hme_bridged = len(bridged)
    qc.bridged_hme_sample = bridged[:40]
    qc.n_hist_marks = sum(
        1
        for n in after_nodes
        if _is_histone_mark(n) and not is_virtual_gin_node(n)
    )

    out = [(a, b, w) for (a, b), w in sorted(edge_set.items())]
    return out, qc


def repair_hme_chains_with_node_lists(
    edges: Sequence[tuple[str, str, float]],
    ed: IndexedEdges,
    node_list_dir: str | Path,
    *,
    max_real_hops: int = 3,
    upstream_graph_hops: int = 4,
    repair_eraser_chains: bool = True,
    eraser_max_rounds: int = 3,
    prune_sourceless_histones: bool = True,
    backfill_active_hist_marks: bool = True,
    backfill_seed_genes: set[str] | None = None,
    do_drop_nonseed_promoter_active: bool = True,
    anchor_nodes: set[str] | None = None,
    n_jobs: int = 1,
) -> tuple[list[tuple[str, str, float]], HmeChainRepairQC]:
    hme = load_hme_symbols(node_list_dir)
    tf = load_tf_symbol_pool(Path(node_list_dir))
    return repair_hme_chains_on_edges(
        edges,
        ed,
        hme_symbols=hme,
        tf_symbols=tf,
        max_real_hops=max_real_hops,
        upstream_graph_hops=upstream_graph_hops,
        repair_eraser_chains=repair_eraser_chains,
        eraser_max_rounds=eraser_max_rounds,
        prune_sourceless_histones=prune_sourceless_histones,
        backfill_active_hist_marks=backfill_active_hist_marks,
        backfill_seed_genes=backfill_seed_genes,
        do_drop_nonseed_promoter_active=do_drop_nonseed_promoter_active,
        anchor_nodes=anchor_nodes,
        n_jobs=n_jobs,
    )

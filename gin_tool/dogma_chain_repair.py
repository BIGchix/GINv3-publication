"""Post-prune central-dogma chain repair for context subgraphs.

After ancestor-only prune, many ``*.protein`` nodes remain without their
``.RNA`` / ``.promoter.active`` chain. This module:

1. Restores each incomplete gene's transcription/translation chain from the
   full GIN when those nodes exist.
2. Repairs genes that already have promoter/RNA/protein stubs but are missing
   the ``promoter.active → *;transcription → .RNA`` link (common after prune).
3. Optionally restores chains for ``promoter.active``-only genes when the full
   dogma exists in GIN (``repair_promoter_only``, **default off**): many
   orphan promoters are miRNA-host loci kept only because
   ``promoter.active → miRNA → seed``, not via the host gene's own RNA/protein.
4. From ``GENE.promoter.active``, walks **reverse real-hop 1** (graph-hop 2 on
   bipartite GIN) and, if any reached TF is already in the subgraph, keeps the
   connecting path nodes/edges (``TF → … → promoter.active``).

Chains are kept even when no subgraph TF is reached, so promoter/RNA places can
act as token sources under ``mark_mode=ancestor_protein`` / ``source_protein``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy import sparse

from gin_tool.active_tf import gin_node_to_gene_symbol
from gin_tool.io import IndexedEdges, directed_csr_from_indexed
from gin_tool.mapping import is_virtual_gin_node


def _outs(adj: sparse.csr_matrix, i: int) -> list[int]:
    lo, hi = adj.indptr[i], adj.indptr[i + 1]
    return [int(x) for x in adj.indices[lo:hi]]


def _gene_from_protein(name: str) -> str | None:
    if not name.lower().endswith(".protein"):
        return None
    if ";" in name:
        return None
    return name[: -len(".protein")]


def _gene_from_promoter_active(name: str) -> str | None:
    if not name.endswith(".promoter.active"):
        return None
    return name[: -len(".promoter.active")]


def _dogma_candidate_names(gene: str) -> list[str]:
    return [
        f"{gene}.protein",
        f"{gene}.RNA",
        f"{gene}.promoter.active",
        f"{gene}.promoter;state",
        f"{gene}.protein;transcription",
        f"{gene}.protein;translation",
    ]


def incomplete_protein_genes(subgraph_nodes: set[str]) -> list[str]:
    """Genes with ``.protein`` in the subgraph but missing ``.RNA`` or promoter."""
    low = {n.lower() for n in subgraph_nodes}
    genes: list[str] = []
    seen: set[str] = set()
    for n in subgraph_nodes:
        gene = _gene_from_protein(n)
        if not gene or gene in seen:
            continue
        seen.add(gene)
        has_rna = f"{gene}.rna".lower() in low
        has_pa = f"{gene}.promoter.active".lower() in low
        if not (has_rna and has_pa):
            genes.append(gene)
    return sorted(genes)


def promoter_only_genes(subgraph_nodes: set[str]) -> list[str]:
    """Genes with ``.promoter.active`` but neither ``.RNA`` nor ``.protein``."""
    genes: list[str] = []
    for n in subgraph_nodes:
        gene = _gene_from_promoter_active(n)
        if not gene:
            continue
        if f"{gene}.RNA" in subgraph_nodes or f"{gene}.protein" in subgraph_nodes:
            continue
        genes.append(gene)
    return sorted(set(genes))


def _adjacency_from_edges(
    edges: Sequence[tuple[str, str, float]] | dict[tuple[str, str], float],
) -> dict[str, list[str]]:
    outs: dict[str, list[str]] = {}
    if isinstance(edges, dict):
        items = edges.items()
        for (a, b), _w in items:
            outs.setdefault(a, []).append(b)
    else:
        for a, b, _w in edges:
            outs.setdefault(a, []).append(b)
    return outs


def has_transcription_path(
    gene: str,
    subgraph_nodes: set[str],
    edge_outs: dict[str, list[str]],
) -> bool:
    """True if ``promoter.active → (transcription virtual) → .RNA`` exists."""
    pa = f"{gene}.promoter.active"
    rna = f"{gene}.RNA"
    if pa not in subgraph_nodes or rna not in subgraph_nodes:
        return False
    for mid in edge_outs.get(pa, []):
        if mid == rna:
            return True
        if not is_virtual_gin_node(mid):
            continue
        if "transcription" not in mid.lower():
            continue
        if rna in edge_outs.get(mid, []):
            return True
    return False


def genes_missing_transcription(
    subgraph_nodes: set[str],
    edges: Sequence[tuple[str, str, float]] | dict[tuple[str, str], float],
) -> list[str]:
    """
    Genes that have ``promoter.active`` and ``.RNA`` (and usually protein) but
    lack the transcription virtual path connecting them.
    """
    outs = _adjacency_from_edges(edges)
    genes: list[str] = []
    for n in subgraph_nodes:
        gene = _gene_from_promoter_active(n)
        if not gene:
            continue
        if f"{gene}.RNA" not in subgraph_nodes:
            continue
        if has_transcription_path(gene, subgraph_nodes, outs):
            continue
        genes.append(gene)
    return sorted(set(genes))


def _collect_dogma_nodes(
    gene: str,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
    adj_rev: sparse.csr_matrix,
) -> set[str]:
    """Strict central-dogma names only (no PPI / complex virtuals).

    Allowlist: ``.protein`` / ``.RNA`` / ``.promoter.active`` / ``.promoter;state``
    plus ``.protein;transcription`` / ``.protein;translation``. Optionally pick up
    a transcription/translation virtual that sits between those reals if named
    differently but still matches ``{gene}.*;(transcription|translation)``.
    """
    _ = index_to_name, adj, adj_rev  # kept for API stability / future use
    nodes = {n for n in _dogma_candidate_names(gene) if n in name_to_index}
    # Tight extra: only gene-prefixed transcription/translation virtuals.
    for n in _dogma_candidate_names(gene):
        if n not in name_to_index:
            continue
        i = name_to_index[n]
        for v in _outs(adj, i) + _outs(adj_rev, i):
            vn = index_to_name[v]
            if not is_virtual_gin_node(vn):
                continue
            lower = vn.lower()
            if not vn.startswith(f"{gene}."):
                continue
            if "transcription" in lower or "translation" in lower:
                nodes.add(vn)
    return nodes


def _reverse_real_hop1_tf_paths(
    promoter_active: str,
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj_rev: sparse.csr_matrix,
    subgraph_tf_nodes: set[str],
) -> tuple[set[str], list[str]]:
    """
    Reverse graph-hop ≤2 from ``promoter.active``; keep paths that hit a
    subgraph TF. Returns (path_nodes, hit_tf_names).
    """
    if promoter_active not in name_to_index:
        return set(), []
    start = name_to_index[promoter_active]
    path_nodes: set[str] = {promoter_active}
    hits: list[str] = []
    for v in _outs(adj_rev, start):
        vn = index_to_name[v]
        for r in _outs(adj_rev, v):
            rn = index_to_name[r]
            if is_virtual_gin_node(rn):
                continue
            if rn in subgraph_tf_nodes:
                path_nodes.add(vn)
                path_nodes.add(rn)
                hits.append(rn)
    return path_nodes, sorted(set(hits))


def _induced_named_edges(
    node_names: set[str],
    *,
    name_to_index: dict[str, int],
    index_to_name: Sequence[str],
    adj: sparse.csr_matrix,
) -> list[tuple[str, str, float]]:
    """Directed edges with both ends in ``node_names`` (CSR neighborhood scan)."""
    if not node_names:
        return []
    idx = [name_to_index[n] for n in node_names if n in name_to_index]
    if not idx:
        return []
    idx_set = set(idx)
    out: list[tuple[str, str, float]] = []
    seen: set[tuple[int, int]] = set()
    for si in idx:
        for ti in _outs(adj, si):
            if ti not in idx_set:
                continue
            key = (si, ti)
            if key in seen:
                continue
            seen.add(key)
            out.append((index_to_name[si], index_to_name[ti], float(adj[si, ti])))
    return out


def subgraph_tf_nodes(subgraph_nodes: set[str], tf_pool: set[str]) -> set[str]:
    """Real nodes in the subgraph whose gene symbol is in the TF pool."""
    out: set[str] = set()
    for n in subgraph_nodes:
        if is_virtual_gin_node(n):
            continue
        sym = gin_node_to_gene_symbol(n)
        if sym and sym in tf_pool:
            out.add(n)
    return out


@dataclass
class DogmaChainRepairQC:
    n_incomplete_proteins: int = 0
    n_missing_transcription: int = 0
    n_promoter_only: int = 0
    n_repaired_genes: int = 0
    n_repaired_transcription: int = 0
    n_repaired_promoter_only: int = 0
    n_genes_no_dogma_in_gin: int = 0
    n_with_tf_link: int = 0
    n_without_tf_link: int = 0
    n_nodes_added: int = 0
    n_edges_added: int = 0
    repaired_genes: list[str] = field(default_factory=list)
    transcription_repaired_genes: list[str] = field(default_factory=list)
    promoter_only_repaired_genes: list[str] = field(default_factory=list)
    tf_linked_genes: list[str] = field(default_factory=list)
    sample_tf_hits: dict[str, list[str]] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "n_incomplete_proteins": self.n_incomplete_proteins,
            "n_missing_transcription": self.n_missing_transcription,
            "n_promoter_only": self.n_promoter_only,
            "n_repaired_genes": self.n_repaired_genes,
            "n_repaired_transcription": self.n_repaired_transcription,
            "n_repaired_promoter_only": self.n_repaired_promoter_only,
            "n_genes_no_dogma_in_gin": self.n_genes_no_dogma_in_gin,
            "n_with_tf_link": self.n_with_tf_link,
            "n_without_tf_link": self.n_without_tf_link,
            "n_nodes_added": self.n_nodes_added,
            "n_edges_added": self.n_edges_added,
            "repaired_genes": self.repaired_genes,
            "transcription_repaired_genes": self.transcription_repaired_genes,
            "promoter_only_repaired_genes": self.promoter_only_repaired_genes,
            "tf_linked_genes": self.tf_linked_genes,
            "sample_tf_hits": self.sample_tf_hits,
        }


def _repair_one_gene_payload(
    gene: str,
    tags: frozenset[str],
    *,
    name_to: dict[str, int],
    idx_to: Sequence[str],
    adj: sparse.csr_matrix,
    adj_rev: sparse.csr_matrix,
    sub_tfs: set[str],
    nodes: set[str],
) -> dict:
    """Return edges/keep/hits for one gene (serial or parallel worker body)."""
    dogma = _collect_dogma_nodes(gene, name_to, idx_to, adj, adj_rev)
    has_chain = any(
        n in dogma
        for n in (f"{gene}.RNA", f"{gene}.promoter.active", f"{gene}.promoter;state")
    )
    tx = f"{gene}.protein;transcription"
    if not has_chain and "missing_transcription" not in tags:
        return {
            "gene": gene,
            "tags": set(tags),
            "skip_no_dogma": True,
            "edges": [],
            "hits": [],
        }
    if "missing_transcription" in tags and tx not in name_to and not has_chain:
        return {
            "gene": gene,
            "tags": set(tags),
            "skip_no_dogma": True,
            "edges": [],
            "hits": [],
        }

    pa = f"{gene}.promoter.active"
    path_nodes: set[str] = set()
    hits: list[str] = []
    if pa in name_to:
        path_nodes, hits = _reverse_real_hop1_tf_paths(
            pa,
            name_to_index=name_to,
            index_to_name=idx_to,
            adj_rev=adj_rev,
            subgraph_tf_nodes=sub_tfs,
        )

    keep = set(dogma) | path_nodes
    prot = f"{gene}.protein"
    if prot in nodes or prot in name_to:
        keep.add(prot)
    rna = f"{gene}.RNA"
    if rna in nodes or rna in name_to:
        keep.add(rna)
    if pa in nodes or pa in name_to:
        keep.add(pa)

    edges = _induced_named_edges(
        keep,
        name_to_index=name_to,
        index_to_name=idx_to,
        adj=adj,
    )
    return {
        "gene": gene,
        "tags": set(tags),
        "skip_no_dogma": False,
        "edges": edges,
        "hits": hits,
    }


# Process-pool worker state (set by initializer; one CSR rebuild per worker).
_DOGMA_POOL: dict = {}


def _dogma_pool_init(
    ed: IndexedEdges,
    sub_tfs: frozenset[str],
    nodes: frozenset[str],
) -> None:
    adj = directed_csr_from_indexed(ed)
    _DOGMA_POOL.clear()
    _DOGMA_POOL.update(
        {
            "ed": ed,
            "adj": adj,
            "adj_rev": adj.T.tocsr(),
            "sub_tfs": set(sub_tfs),
            "nodes": set(nodes),
        }
    )


def _dogma_pool_worker(item: tuple[str, frozenset[str]]) -> dict:
    gene, tags = item
    ed: IndexedEdges = _DOGMA_POOL["ed"]
    return _repair_one_gene_payload(
        gene,
        tags,
        name_to=ed.name_to_index,
        idx_to=ed.index_to_name,
        adj=_DOGMA_POOL["adj"],
        adj_rev=_DOGMA_POOL["adj_rev"],
        sub_tfs=_DOGMA_POOL["sub_tfs"],
        nodes=_DOGMA_POOL["nodes"],
    )


def repair_dogma_chains_on_edges(
    edges: Sequence[tuple[str, str, float]],
    ed: IndexedEdges,
    *,
    tf_pool: set[str],
    real_hops: int = 1,
    repair_missing_transcription: bool = True,
    repair_promoter_only: bool = False,
    n_jobs: int = 1,
) -> tuple[list[tuple[str, str, float]], DogmaChainRepairQC]:
    """
    Augment a directed subgraph edge list with missing transcription chains.

    ``real_hops`` currently only supports 1 (reverse graph-hop = 2).
    ``repair_promoter_only`` defaults to False: orphan promoters are often
    miRNA hosts whose own protein is not an ancestor of seeds.

    ``n_jobs``: CPU workers for per-gene repair (``1`` = serial; ``>1`` process pool).
    """
    if real_hops != 1:
        raise ValueError("dogma_chain_repair currently supports real_hops=1 only")

    from gin_tool.subgraph_pipeline.parallel_util import map_processes, resolve_n_jobs

    qc = DogmaChainRepairQC()
    edge_set = {(a, b): float(w) for a, b, w in edges}
    n_edges_before = len(edge_set)
    nodes = {a for a, _, _ in edges} | {b for _, b, _ in edges}
    before_nodes = set(nodes)

    incomplete = incomplete_protein_genes(nodes)
    missing_tx = (
        genes_missing_transcription(nodes, edge_set) if repair_missing_transcription else []
    )
    prom_only = promoter_only_genes(nodes) if repair_promoter_only else []
    qc.n_incomplete_proteins = len(incomplete)
    qc.n_missing_transcription = len(missing_tx)
    qc.n_promoter_only = len(prom_only)

    reason: dict[str, set[str]] = {}
    for g in incomplete:
        reason.setdefault(g, set()).add("incomplete_protein")
    for g in missing_tx:
        reason.setdefault(g, set()).add("missing_transcription")
    for g in prom_only:
        reason.setdefault(g, set()).add("promoter_only")

    targets = sorted(reason)
    if not targets:
        return [(a, b, w) for (a, b), w in sorted(edge_set.items())], qc

    sub_tfs = subgraph_tf_nodes(nodes, tf_pool)
    jobs = resolve_n_jobs(n_jobs)
    payloads = [(g, frozenset(reason[g])) for g in targets]

    if jobs <= 1:
        adj = directed_csr_from_indexed(ed)
        adj_rev = adj.T.tocsr()
        name_to = ed.name_to_index
        idx_to = ed.index_to_name
        results = [
            _repair_one_gene_payload(
                g,
                tags,
                name_to=name_to,
                idx_to=idx_to,
                adj=adj,
                adj_rev=adj_rev,
                sub_tfs=sub_tfs,
                nodes=nodes,
            )
            for g, tags in payloads
        ]
    else:
        results = map_processes(
            _dogma_pool_worker,
            payloads,
            n_jobs=jobs,
            initializer=_dogma_pool_init,
            initargs=(ed, frozenset(sub_tfs), frozenset(nodes)),
        )

    for res in results:
        gene = res["gene"]
        tags = res["tags"]
        if res["skip_no_dogma"]:
            qc.n_genes_no_dogma_in_gin += 1
            continue
        for e_a, e_b, e_w in res["edges"]:
            edge_set.setdefault((e_a, e_b), e_w)
            nodes.add(e_a)
            nodes.add(e_b)

        qc.n_repaired_genes += 1
        qc.repaired_genes.append(gene)
        if "missing_transcription" in tags:
            qc.n_repaired_transcription += 1
            qc.transcription_repaired_genes.append(gene)
        if "promoter_only" in tags:
            qc.n_repaired_promoter_only += 1
            qc.promoter_only_repaired_genes.append(gene)
        hits = res["hits"]
        if hits:
            qc.n_with_tf_link += 1
            qc.tf_linked_genes.append(gene)
            if len(qc.sample_tf_hits) < 12:
                qc.sample_tf_hits[gene] = hits[:8]
        else:
            qc.n_without_tf_link += 1

    after_nodes = {a for a, _ in edge_set} | {b for _, b in edge_set}
    qc.n_nodes_added = len(after_nodes - before_nodes)
    qc.n_edges_added = len(edge_set) - n_edges_before

    out = [(a, b, w) for (a, b), w in sorted(edge_set.items())]
    return out, qc

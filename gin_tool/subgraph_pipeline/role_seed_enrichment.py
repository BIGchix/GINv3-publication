"""Context v4: early role soft-seeds (**HME protein + RNA** by default).

Soft seeds join primary seeds for context BFS / ancestor prune so HME enzymes
and their cognate transcripts enter the expansion universe before topological
qualification.

Metabolic hyperedge enzymes are **not** soft-seeded here by default; they are
recovered after LCC from in-subgraph hyperedges (``met_enzyme_recovery``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from gin_tool.hme_chain_repair import load_hme_symbols
from gin_tool.io import IndexedEdges, directed_csr_from_indexed
from gin_tool.mapping import is_virtual_gin_node
from gin_tool.spn_export import _is_metabolite

_DEFAULT_HUB_METS = frozenset(
    {
        "SID8145663",
        "SID3303",
        "SID315674483",
        "SID3380",
        "SID3305",
        "SID8148096",
        "SID315675471",
        "SID85646635",
        "SID3308",
        "SID8145742",
    }
)


def _outs(adj, i: int) -> list[int]:
    return [int(j) for j in adj.indices[adj.indptr[i] : adj.indptr[i + 1]]]


def proteins_from_hyperedge_name(name: str) -> set[str]:
    """Parse ``*.protein`` tokens embedded in a ``;``-joined virtual name."""
    out: set[str] = set()
    for part in name.split(";"):
        p = part.strip()
        if p.endswith(".protein"):
            out.add(p)
    return out


def _neighbor_indices(adj, adj_rev, i: int) -> set[int]:
    return set(_outs(adj, i)) | set(_outs(adj_rev, i))


@dataclass
class RoleSeedEnrichmentQC:
    n_primary_seeds: int = 0
    n_hme_soft: int = 0
    n_met_enzyme_soft: int = 0
    n_hub_mets_skipped: int = 0
    hme_soft_sample: list[str] = field(default_factory=list)
    met_enzyme_soft_sample: list[str] = field(default_factory=list)
    hub_mets_skipped: list[str] = field(default_factory=list)
    soft_seeds: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_primary_seeds": self.n_primary_seeds,
            "n_hme_soft": self.n_hme_soft,
            "n_met_enzyme_soft": self.n_met_enzyme_soft,
            "n_hub_mets_skipped": self.n_hub_mets_skipped,
            "hme_soft_sample": self.hme_soft_sample,
            "met_enzyme_soft_sample": self.met_enzyme_soft_sample,
            "hub_mets_skipped": self.hub_mets_skipped,
            "n_soft_seeds": len(self.soft_seeds),
            "soft_seeds_sample": self.soft_seeds[:40],
        }


def collect_hme_soft_seeds(
    ed: IndexedEdges,
    *,
    node_list_dir: str | Path,
) -> list[str]:
    """HME ``.protein`` and ``.RNA`` places present in GIN.

    Soft-seeding both keeps RNA→protein→hist paths available for later dogma
    wiring and RNAseq place-level injection.
    """
    symbols = load_hme_symbols(node_list_dir)
    out: list[str] = []
    for sym in sorted(symbols):
        base = sym[: -len(".protein")] if sym.endswith(".protein") else sym
        prot = sym if sym.endswith(".protein") else f"{sym}.protein"
        rna = f"{base}.RNA"
        if prot in ed.name_to_index:
            out.append(prot)
        if rna in ed.name_to_index:
            out.append(rna)
    return out


def collect_met_enzyme_soft_seeds(
    ed: IndexedEdges,
    seed_names: Sequence[str],
    *,
    hub_mets: Iterable[str] | None = None,
    max_enzymes_per_met: int = 40,
    max_total_enzymes: int = 400,
    n_jobs: int = 1,
) -> tuple[list[str], list[str]]:
    """Enzymes from hyperedges adjacent to non-hub seed metabolites.

    Returns ``(enzyme_proteins, skipped_hub_mets)``.
    """
    from gin_tool.subgraph_pipeline.parallel_util import map_threads, resolve_n_jobs

    hubs = frozenset(hub_mets) if hub_mets is not None else _DEFAULT_HUB_METS
    adj = directed_csr_from_indexed(ed)
    adj_rev = adj.T.tocsr()
    skipped: list[str] = []
    enzymes: set[str] = set()

    seed_mets = [
        n
        for n in seed_names
        if _is_metabolite(n) and n in ed.name_to_index
    ]
    work_mets: list[str] = []
    for met in seed_mets:
        if met in hubs:
            skipped.append(met)
        else:
            work_mets.append(met)

    def _one(met: str) -> list[str]:
        i = ed.name_to_index[met]
        per_met: set[str] = set()
        for j in _neighbor_indices(adj, adj_rev, i):
            vname = ed.index_to_name[j]
            if is_virtual_gin_node(vname) or ";" in vname:
                per_met |= proteins_from_hyperedge_name(vname)
                for k in _neighbor_indices(adj, adj_rev, j):
                    kn = ed.index_to_name[k]
                    if kn.endswith(".protein") and ";" not in kn:
                        per_met.add(kn)
            elif vname.endswith(".protein") and ";" not in vname:
                per_met.add(vname)
        ranked = sorted(p for p in per_met if p in ed.name_to_index)
        return ranked[: max(0, int(max_enzymes_per_met))]

    jobs = resolve_n_jobs(n_jobs)
    parts = map_threads(_one, work_mets, n_jobs=jobs) if work_mets else []
    for ranked in parts:
        for p in ranked:
            enzymes.add(p)
            if len(enzymes) >= max_total_enzymes:
                return sorted(enzymes), skipped
    return sorted(enzymes), skipped


def enrich_role_soft_seeds(
    ed: IndexedEdges,
    primary_seeds: Sequence[str],
    *,
    node_list_dir: str | Path,
    include_hme: bool = True,
    include_met_enzymes: bool = False,
    hub_mets: Iterable[str] | None = None,
    max_enzymes_per_met: int = 40,
    max_total_enzymes: int = 400,
    n_jobs: int = 1,
) -> tuple[list[str], RoleSeedEnrichmentQC]:
    """Return ``primary ∪ soft`` (stable order) and QC."""
    primary = list(dict.fromkeys(primary_seeds))
    qc = RoleSeedEnrichmentQC(n_primary_seeds=len(primary))
    soft: list[str] = []
    primary_set = set(primary)

    if include_hme:
        hme = collect_hme_soft_seeds(ed, node_list_dir=node_list_dir)
        for p in hme:
            if p not in primary_set:
                soft.append(p)
        qc.n_hme_soft = len(hme)
        qc.hme_soft_sample = hme[:40]

    if include_met_enzymes:
        enz, skipped = collect_met_enzyme_soft_seeds(
            ed,
            primary,
            hub_mets=hub_mets,
            max_enzymes_per_met=max_enzymes_per_met,
            max_total_enzymes=max_total_enzymes,
            n_jobs=n_jobs,
        )
        for p in enz:
            if p not in primary_set and p not in soft:
                soft.append(p)
        qc.n_met_enzyme_soft = len(enz)
        qc.met_enzyme_soft_sample = enz[:40]
        qc.n_hub_mets_skipped = len(skipped)
        qc.hub_mets_skipped = skipped

    qc.soft_seeds = list(soft)
    merged = primary + soft
    return merged, qc

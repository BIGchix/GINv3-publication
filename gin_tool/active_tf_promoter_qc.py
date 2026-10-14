"""QC: how many seed promoters are reachable from the active TF list (full GIN)."""

from __future__ import annotations

import json
import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

from gin_tool.io import iter_sifi_edges


WARN_SUGGESTION = (
    "Manually expand the active TF list (--active-tfs / --active-tfs-file / "
    "--active-tf-file) with additional TFs that regulate the uncovered seed "
    "promoters, then re-run."
)


@dataclass
class PromoterCoverageHopQC:
    max_graph_hops: int
    n_seed_genes: int
    n_active_tfs: int
    n_genes_covered: int
    uncovered_genes: list[str] = field(default_factory=list)
    per_gene_n_tfs: dict[str, int] = field(default_factory=dict)
    n_tfs_covering_ge1: int = 0

    @property
    def frac_covered(self) -> float:
        if self.n_seed_genes <= 0:
            return 1.0
        return self.n_genes_covered / self.n_seed_genes

    def to_dict(self) -> dict:
        return {
            "tag": f"graph-hop={self.max_graph_hops}",
            "n_seed_genes": self.n_seed_genes,
            "n_active_tfs": self.n_active_tfs,
            "n_genes_covered": self.n_genes_covered,
            "frac_covered": self.frac_covered,
            "n_uncovered": len(self.uncovered_genes),
            "uncovered_genes": list(self.uncovered_genes),
            "n_tfs_covering_ge1": self.n_tfs_covering_ge1,
            "per_gene_n_tfs": dict(self.per_gene_n_tfs),
        }


@dataclass
class ActiveTFSeedPromoterQC:
    """Full-GIN reachability of seed promoters from active TFs."""

    n_active_tfs: int
    n_seed_genes: int
    seed_genes: list[str]
    active_tf_symbols: list[str]
    missing_promoter_nodes: list[str]
    graphhop1: PromoterCoverageHopQC
    graphhop2: PromoterCoverageHopQC
    warn_max_graph_hops: int = 2

    @property
    def primary(self) -> PromoterCoverageHopQC:
        return self.graphhop2 if self.warn_max_graph_hops >= 2 else self.graphhop1

    @property
    def has_uncovered(self) -> bool:
        return bool(self.primary.uncovered_genes)

    def warning_message(self) -> str | None:
        hop = self.primary
        if not hop.uncovered_genes:
            return None
        sample = ", ".join(hop.uncovered_genes[:12])
        more = ""
        if len(hop.uncovered_genes) > 12:
            more = f" (+{len(hop.uncovered_genes) - 12} more)"
        return (
            f"[WARN] active-TF→seed-promoter QC: "
            f"{len(hop.uncovered_genes)}/{hop.n_seed_genes} seed promoters "
            f"unreachable from {self.n_active_tfs} active TFs at "
            f"graph-hop≤{hop.max_graph_hops} "
            f"(uncovered: {sample}{more}). {WARN_SUGGESTION}"
        )

    def to_dict(self) -> dict:
        return {
            "n_active_tfs": self.n_active_tfs,
            "n_seed_genes": self.n_seed_genes,
            "active_tf_symbols": list(self.active_tf_symbols),
            "seed_genes": list(self.seed_genes),
            "missing_promoter_nodes": list(self.missing_promoter_nodes),
            "warn_max_graph_hops": self.warn_max_graph_hops,
            "graphhop1": self.graphhop1.to_dict(),
            "graphhop2": self.graphhop2.to_dict(),
            "has_uncovered": self.has_uncovered,
            "warning": self.warning_message(),
            "suggestion": WARN_SUGGESTION if self.has_uncovered else None,
        }


def _normalize_tf_sources(active_tfs: Iterable[str]) -> tuple[list[str], set[str]]:
    """Return (sorted symbols, source node names used as BFS roots)."""
    symbols: list[str] = []
    sources: set[str] = set()
    for raw in active_tfs:
        t = str(raw).strip()
        if not t or t.startswith("#"):
            continue
        if t.endswith(".protein"):
            sym = t[: -len(".protein")]
            symbols.append(sym)
            sources.add(t)
            sources.add(sym)
        else:
            symbols.append(t)
            sources.add(t)
            sources.add(f"{t}.protein")
    # unique preserve order
    seen: set[str] = set()
    uniq_syms: list[str] = []
    for s in symbols:
        if s not in seen:
            seen.add(s)
            uniq_syms.append(s)
    return uniq_syms, sources


def _seed_promoter_targets(seed_genes: Sequence[str]) -> dict[str, set[str]]:
    """Map gene → promoter nodes that count as 'connected'."""
    out: dict[str, set[str]] = {}
    for g in seed_genes:
        g = str(g).strip()
        if not g or g.startswith("#"):
            continue
        out[g] = {f"{g}.promoter.active", f"{g}.promoter;state"}
    return out


def qc_active_tf_seed_promoter_connectivity(
    sifi_path: str | Path,
    active_tfs: Sequence[str],
    seed_genes: Sequence[str],
    *,
    max_lines: int | None = None,
    warn_max_graph_hops: int = 2,
) -> ActiveTFSeedPromoterQC:
    """
    Count seed promoters reachable from active TFs on the full directed GIN.

    Path of interest is typically:
      TF(.protein) → Gene.promoter;state → Gene.promoter.active
    so graph-hop≤1 covers ``;state`` and ≤2 covers ``.promoter.active``.
    """
    tf_syms, tf_sources = _normalize_tf_sources(active_tfs)
    targets_by_gene = _seed_promoter_targets(seed_genes)
    seed_gene_list = sorted(targets_by_gene)
    all_targets = set().union(*targets_by_gene.values()) if targets_by_gene else set()

    out_adj: dict[str, set[str]] = defaultdict(set)
    seen_promoter_nodes: set[str] = set()
    for s, t, _w in iter_sifi_edges(str(sifi_path), max_lines=max_lines):
        if s in tf_sources or s in out_adj:
            out_adj[s].add(t)
        if t in all_targets:
            seen_promoter_nodes.add(t)
        if s in all_targets:
            seen_promoter_nodes.add(s)

    # Second pass: expand hop-1 neighborhood outbound (for hop=2)
    hop1_nodes: set[str] = set()
    for src in list(tf_sources):
        hop1_nodes |= out_adj.get(src, set())
    if hop1_nodes:
        for s, t, _w in iter_sifi_edges(str(sifi_path), max_lines=max_lines):
            if s in hop1_nodes:
                out_adj[s].add(t)

    def cover_at(max_hops: int) -> PromoterCoverageHopQC:
        gene_tfs: dict[str, set[str]] = {g: set() for g in seed_gene_list}
        # BFS from each TF source up to max_hops
        for src in tf_sources:
            # map source back to symbol
            sym = src[: -len(".protein")] if src.endswith(".protein") else src
            if sym not in tf_syms and src.endswith(".protein"):
                sym = src[: -len(".protein")]
            frontier = {src}
            visited = {src}
            for _depth in range(max_hops):
                nxt: set[str] = set()
                for u in frontier:
                    for v in out_adj.get(u, ()):
                        if v in all_targets:
                            for g, tset in targets_by_gene.items():
                                if v in tset:
                                    gene_tfs[g].add(sym)
                        if v not in visited:
                            visited.add(v)
                            nxt.add(v)
                frontier = nxt
                if not frontier:
                    break
        per_gene = {g: len(gene_tfs[g]) for g in seed_gene_list}
        uncovered = sorted(g for g in seed_gene_list if per_gene[g] == 0)
        n_tfs_ge1 = len({t for gs in gene_tfs.values() for t in gs})
        return PromoterCoverageHopQC(
            max_graph_hops=max_hops,
            n_seed_genes=len(seed_gene_list),
            n_active_tfs=len(tf_syms),
            n_genes_covered=len(seed_gene_list) - len(uncovered),
            uncovered_genes=uncovered,
            per_gene_n_tfs=per_gene,
            n_tfs_covering_ge1=n_tfs_ge1,
        )

    expected_active = {f"{g}.promoter.active" for g in seed_gene_list}
    missing = sorted(expected_active - seen_promoter_nodes)

    return ActiveTFSeedPromoterQC(
        n_active_tfs=len(tf_syms),
        n_seed_genes=len(seed_gene_list),
        seed_genes=seed_gene_list,
        active_tf_symbols=tf_syms,
        missing_promoter_nodes=missing,
        graphhop1=cover_at(1),
        graphhop2=cover_at(2),
        warn_max_graph_hops=warn_max_graph_hops,
    )


def emit_active_tf_promoter_qc(
    qc: ActiveTFSeedPromoterQC,
    *,
    out_json: Path | None = None,
    print_warning: bool = True,
    use_warnings_module: bool = False,
) -> list[str]:
    """
    Persist QC JSON and emit an immediate warning if any seed promoter is uncovered.

    Returns warning strings (empty if all covered). Does not abort the pipeline.
    """
    if out_json is not None:
        out_json = Path(out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(
            json.dumps(qc.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    msg = qc.warning_message()
    if not msg:
        ok = (
            f"[QC] active-TF→seed-promoter: "
            f"{qc.primary.n_genes_covered}/{qc.primary.n_seed_genes} covered "
            f"at graph-hop≤{qc.primary.max_graph_hops} "
            f"(active_tfs={qc.n_active_tfs})"
        )
        if print_warning:
            print(ok, flush=True)
        return []
    if print_warning:
        print(msg, flush=True)
    if use_warnings_module:
        warnings.warn(msg, UserWarning, stacklevel=2)
    return [msg]


def seed_genes_from_gin_nodes(seed_gin_nodes: Sequence[str]) -> list[str]:
    """Extract gene symbols that have a ``.promoter`` seed anchor among GIN nodes."""
    genes: list[str] = []
    seen: set[str] = set()
    for n in seed_gin_nodes:
        lower = n.lower()
        if ".promoter" not in lower:
            continue
        g = n.split(".promoter", 1)[0]
        if g and g not in seen:
            seen.add(g)
            genes.append(g)
    return genes

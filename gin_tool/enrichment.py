"""Lightweight ORA-style enrichment against a custom universe (e.g. subgraph nodes)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from scipy.stats import fisher_exact


@dataclass
class EnrichmentRow:
    term_id: str
    k: int  # hits in query
    K: int  # query size
    n: int  # term size in universe
    N: int  # universe size
    odds_ratio: float
    pvalue: float


def ora_fisher(
    query_genes: set[str],
    term_genes: set[str],
    universe: set[str],
) -> EnrichmentRow | None:
    """
    2x2 Fisher exact test for one term.

    Contingency (inside universe U):
    a = |Q∩T|, b = |Q\\T|, c = |T\\Q|, d = |U \\ (Q∪T)|
    """
    q = query_genes & universe
    t = term_genes & universe
    if not q or not t:
        return None
    a = len(q & t)
    b = len(q - t)
    c = len(t - q)
    d = len(universe) - len(q | t)
    oddsratio, p = fisher_exact([[a, b], [c, d]])
    return EnrichmentRow(
        term_id="term",
        k=a,
        K=len(q),
        n=len(t),
        N=len(universe),
        odds_ratio=float(oddsratio),
        pvalue=float(p),
    )


def enrich_many_terms(
    query_genes: set[str],
    term_to_genes: dict[str, set[str]],
    universe: set[str],
) -> list[EnrichmentRow]:
    rows: list[EnrichmentRow] = []
    q = query_genes & universe
    for tid, genes in term_to_genes.items():
        t = genes & universe
        if not t:
            continue
        a = len(q & t)
        if a == 0:
            continue
        b = len(q - t)
        c = len(t - q)
        d = len(universe) - len(q | t)
        oddsratio, p = fisher_exact([[a, b], [c, d]])
        rows.append(
            EnrichmentRow(
                term_id=tid,
                k=a,
                K=len(q),
                n=len(t),
                N=len(universe),
                odds_ratio=float(oddsratio),
                pvalue=float(p),
            )
        )
    rows.sort(key=lambda r: r.pvalue)
    return rows


def universe_from_subgraph_nodes(node_names: Iterable[str]) -> set[str]:
    """Strip to gene-like tokens is not done here; universe is literal GIN node names."""
    return {n for n in node_names if n}

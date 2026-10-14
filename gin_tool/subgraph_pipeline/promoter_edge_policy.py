"""Histone / promoter hub edge policy: HME must not act as TFs on ``promoter;state``.

Canonical rule (standard pipeline default):
  drop ``HME.protein → *.promoter;state``; HME only affects hubs via histone marks
  (``HME → mark → ;state``, with mark signs already in GIN).

Optional:
  - scale ``histone_mark → ;state`` weights
  - force DNA-methylation-axis repressors (seed list) ``→ ;state`` weight sign to −1
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

DEFAULT_REPRESSIVE_CHROMATIN_TSV = (
    Path(__file__).resolve().parents[2]
    / "seedLists"
    / "chromatin"
    / "repressive_chromatin_enzymes.tsv"
)


def default_hme_enzyme_genes(
    node_list_dir: str | Path | None = None,
) -> frozenset[str]:
    """Bare gene symbols from ``nodeList/nodeList_HME.txt`` (+ HME pool)."""
    from gin_tool.hme_chain_repair import load_hme_symbols

    if node_list_dir is None:
        # repo_root/nodeList — gin_tool/subgraph_pipeline → parents[2]
        node_list_dir = Path(__file__).resolve().parents[2] / "nodeList"
    return frozenset(load_hme_symbols(node_list_dir))


# Backward-compatible alias (= nodeList_HME at import time).
CHROMATIN_ENZYME_GENES: frozenset[str] = default_hme_enzyme_genes()


def load_repressive_chromatin_genes(
    tsv: str | Path | None = None,
) -> dict[str, float]:
    """gene → suggested_sign (usually −1) from repressive_chromatin_enzymes.tsv."""
    path = Path(tsv) if tsv else DEFAULT_REPRESSIVE_CHROMATIN_TSV
    out: dict[str, float] = {}
    if not path.is_file():
        return out
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        if not line.strip() or line.startswith("#"):
            continue
        if i == 0 and line.lower().startswith("gene"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        gene = parts[0].strip()
        try:
            sign = float(parts[2])
        except ValueError:
            sign = -1.0
        if gene:
            out[gene] = sign if sign != 0 else -1.0
    return out


def _is_promoter_state_target(name: str) -> bool:
    low = name.lower()
    return ".promoter;state" in low or low.endswith("promoter;state")


def _is_promoter_active_target(name: str) -> bool:
    low = name.lower()
    return ".promoter.active" in low or low.endswith("promoter.active")


def _is_promoter_incoming_target(name: str) -> bool:
    return _is_promoter_state_target(name) or _is_promoter_active_target(name)


def _is_histone_mark_place(name: str) -> bool:
    from gin_tool.spn_export import _is_histone_mark

    return bool(_is_histone_mark(name))


def _gene_from_protein_place(name: str) -> str | None:
    if not name.endswith(".protein"):
        return None
    return name[: -len(".protein")]


@dataclass
class PromoterEdgePolicyStats:
    n_in: int = 0
    n_out: int = 0
    n_dropped_enzyme_to_state: int = 0
    n_scaled_hist_to_state: int = 0
    n_repressive_signed: int = 0
    hist_scale: float = 1.0
    incoming_promoter_scale: float = 1.0
    n_scaled_incoming_promoter: int = 0
    dropped_enzyme_genes: list[str] | None = None
    repressive_signed_genes: list[str] | None = None
    enzyme_list_source: str = "nodeList_HME"
    n_enzyme_genes_in_list: int = 0
    sign_repressive_to_state: bool = False

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


def apply_promoter_edge_policy(
    gin_edges: Sequence[tuple[str, str, float]],
    *,
    drop_chromatin_enzyme_to_state: bool = True,
    histone_to_state_scale: float = 1.0,
    promoter_incoming_edge_scale: float = 1.0,
    enzyme_genes: frozenset[str] | set[str] | None = None,
    node_list_dir: str | Path | None = None,
    sign_repressive_to_state: bool = False,
    repressive_genes_tsv: str | Path | None = None,
) -> tuple[list[tuple[str, str, float]], PromoterEdgePolicyStats]:
    """Filter/scale/sign GIN edges before Petri build.

    - Drop ``HME.protein → *.promoter;state`` (default gene list = nodeList_HME).
    - Optionally multiply ``histone_mark → *.promoter;state`` by ``histone_to_state_scale``.
    - Optionally multiply **all** incoming edges to ``*.promoter;state`` / ``*.promoter.active``
      by ``promoter_incoming_edge_scale`` (after hist/repressive transforms).
    - Optionally set DNA-methylation-axis repressor ``protein → ;state`` weight to
      ``sign(suggested) * abs(w)`` (typically −1).
    """
    if enzyme_genes is not None:
        genes = frozenset(enzyme_genes)
        src_label = "caller"
    else:
        genes = default_hme_enzyme_genes(node_list_dir)
        src_label = "nodeList_HME"
    scale = float(histone_to_state_scale)
    incoming_scale = float(promoter_incoming_edge_scale)
    repressive = (
        load_repressive_chromatin_genes(repressive_genes_tsv)
        if sign_repressive_to_state
        else {}
    )
    stats = PromoterEdgePolicyStats(
        n_in=len(gin_edges),
        hist_scale=scale,
        incoming_promoter_scale=incoming_scale,
        enzyme_list_source=src_label,
        n_enzyme_genes_in_list=len(genes),
        sign_repressive_to_state=bool(sign_repressive_to_state),
    )
    dropped_genes: set[str] = set()
    signed_genes: set[str] = set()
    out: list[tuple[str, str, float]] = []

    for src, tgt, w in gin_edges:
        if drop_chromatin_enzyme_to_state and _is_promoter_state_target(tgt):
            g = _gene_from_protein_place(src)
            if g is not None and g in genes:
                stats.n_dropped_enzyme_to_state += 1
                dropped_genes.add(g)
                continue
        wf = float(w)
        if scale != 1.0 and _is_promoter_state_target(tgt) and _is_histone_mark_place(src):
            wf *= scale
            stats.n_scaled_hist_to_state += 1
        if sign_repressive_to_state and _is_promoter_state_target(tgt):
            g = _gene_from_protein_place(src)
            if g is not None and g in repressive:
                sug = float(repressive[g])
                mag = abs(wf) if wf != 0.0 else 1.0
                new_w = (-mag) if sug < 0 else mag
                if new_w != wf:
                    stats.n_repressive_signed += 1
                    signed_genes.add(g)
                wf = new_w
        if incoming_scale != 1.0 and _is_promoter_incoming_target(tgt):
            wf *= incoming_scale
            stats.n_scaled_incoming_promoter += 1
        out.append((src, tgt, wf))

    stats.n_out = len(out)
    stats.dropped_enzyme_genes = sorted(dropped_genes)
    stats.repressive_signed_genes = sorted(signed_genes)
    return out, stats

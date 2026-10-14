"""Map seed inputs (genes + metabolites) to GIN walk anchors."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import re
from typing import Any, Callable, Literal

from gin_tool.gin_defaults import REPO_ROOT
from gin_tool.subgraph_pipeline.manifest import CentralDogmaMappingQC
from gin_tool.subgraph_pipeline.metabolite_id_resolve import (
    GinMetaboliteIndex,
    PubChemMetaboliteResolver,
    resolve_metabolite_input,
)

_DEFAULT_PUBCHEM_CACHE = REPO_ROOT / "runtime" / "cache" / "metabolite_pubchem.json"


SeedKind = Literal["gene", "metabolite"]


@dataclass
class CentralDogmaSeed:
    """One seed gene and GIN anchors along DNA→RNA→protein."""

    symbol: str
    kind: SeedKind = "gene"
    rna_node: str | None = None
    protein_node: str | None = None
    promoter_node: str | None = None  # real: GENE.promoter.active
    metabolite_node: str | None = None

    @property
    def walk_anchor_nodes(self) -> list[str]:
        if self.kind == "metabolite":
            return [self.metabolite_node] if self.metabolite_node else []
        # DNA → RNA → protein (all present; no prefer / first-hit collapse)
        out: list[str] = []
        if self.promoter_node:
            out.append(self.promoter_node)
        if self.rna_node:
            out.append(self.rna_node)
        if self.protein_node:
            out.append(self.protein_node)
        return out

    @property
    def regulatory_anchor_nodes(self) -> list[str]:
        return self.walk_anchor_nodes.copy()

    @property
    def triplet_complete(self) -> bool:
        if self.kind == "metabolite":
            return bool(self.metabolite_node)
        return bool(self.rna_node and self.protein_node and self.promoter_node)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SeedMappingQC(CentralDogmaMappingQC):
    n_metabolite: int = 0
    n_gene: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def looks_like_metabolite_id(name: str) -> bool:
    s = name.strip()
    if not s:
        return False
    upper = s.upper()
    if re.match(r"^CAS\d", upper):
        return True
    if re.match(r"^(CHEBI|SID|HMDB|KEGGC)\d", upper):
        return True
    if re.match(r"^(?:CPD:|KEGG:)?C\d{5}$", upper):
        return True
    lower = s.lower()
    return lower.endswith(".met") or "chebi" in lower


def resolve_metabolite_gin_node(
    symbol: str,
    gin_is_member: Callable[[str], bool],
    metabolite_pool: set[str] | None = None,
    *,
    metabolite_index: GinMetaboliteIndex | None = None,
    pubchem: PubChemMetaboliteResolver | None = None,
    allow_online: bool = True,
) -> str | None:
    """Resolve one metabolite identifier to a GIN node name."""
    raw = symbol.strip()
    if not raw:
        return None
    pool = metabolite_pool or set()
    index = metabolite_index or GinMetaboliteIndex.from_pool(pool)

    result = resolve_metabolite_input(
        raw,
        index,
        pubchem=pubchem,
        allow_online=allow_online and pubchem is not None,
    )
    if result.gin_node and gin_is_member(result.gin_node):
        return result.gin_node

    # legacy direct membership
    if raw in pool and gin_is_member(raw):
        return raw
    upper = raw.upper()
    for cand in (raw, upper):
        if cand in pool and gin_is_member(cand):
            return cand
    if gin_is_member(raw):
        return raw
    if upper != raw and gin_is_member(upper):
        return upper
    return None


def gin_central_dogma_node_names(symbol: str) -> dict[str, str]:
    """
    Candidate GIN names for one HGNC-like symbol.

    Promoter uses real node ``GENE.promoter.active`` (not virtual ``;state``).
    """
    sym = symbol.strip()
    return {
        "rna": f"{sym}.RNA",
        "protein": f"{sym}.protein",
        "protein_legacy": sym,
        "promoter": f"{sym}.promoter.active",
    }


def resolve_protein_node(candidates: dict[str, str], gin_is_member: Callable[[str], bool]) -> str | None:
    if gin_is_member(candidates["protein"]):
        return candidates["protein"]
    legacy = candidates.get("protein_legacy")
    if legacy and gin_is_member(legacy):
        return legacy
    return None


def map_union_genes_central_dogma(
    union_genes: list[str],
    gin_is_member: Callable[[str], bool],
) -> tuple[list[CentralDogmaSeed], CentralDogmaMappingQC]:
    records, qc = map_union_seeds(union_genes, gin_is_member)
    return records, qc


def map_union_seeds(
    union_inputs: list[str],
    gin_is_member: Callable[[str], bool],
    *,
    seed_kinds: dict[str, str] | None = None,
    metabolite_pool: set[str] | None = None,
    metabolite_index: GinMetaboliteIndex | None = None,
    pubchem: PubChemMetaboliteResolver | None = None,
    allow_online_metabolite: bool = True,
    allow_gene_metabolite_fallback: bool = False,
) -> tuple[list[CentralDogmaSeed], SeedMappingQC]:
    """
    Map seed tokens to GIN anchors.

    When ``seed_kinds`` is provided (preferred), each token is processed strictly
    as ``gene`` or ``metabolite``. Heuristic pool / id-shape detection is only used
    when ``seed_kinds`` is omitted (legacy mixed lists).
    """
    records: list[CentralDogmaSeed] = []
    unmapped: list[str] = []
    pool = metabolite_pool or set()
    index = metabolite_index or GinMetaboliteIndex.from_pool(pool)
    resolver = pubchem
    if allow_online_metabolite and resolver is None and pool:
        resolver = PubChemMetaboliteResolver(cache_path=_DEFAULT_PUBCHEM_CACHE)

    for raw in union_inputs:
        symbol = str(raw).strip()
        if not symbol:
            continue

        kind: SeedKind | None = None
        if seed_kinds is not None:
            k = seed_kinds.get(symbol) or seed_kinds.get(symbol.upper())
            if k in ("gene", "metabolite"):
                kind = k  # type: ignore[assignment]

        if kind is None:
            in_pool = symbol in pool or symbol.upper() in pool
            looks_met = looks_like_metabolite_id(symbol)
            kind = "metabolite" if (looks_met or in_pool) else "gene"

        if kind == "metabolite":
            node = resolve_metabolite_gin_node(
                symbol,
                gin_is_member,
                pool,
                metabolite_index=index,
                pubchem=resolver,
                allow_online=allow_online_metabolite,
            )
            rec = CentralDogmaSeed(symbol=symbol, kind="metabolite", metabolite_node=node)
            records.append(rec)
            if not node:
                unmapped.append(symbol)
            continue

        gene_sym = symbol.split(".")[0].strip().upper()
        cands = gin_central_dogma_node_names(gene_sym)
        rec = CentralDogmaSeed(
            symbol=gene_sym,
            kind="gene",
            rna_node=cands["rna"] if gin_is_member(cands["rna"]) else None,
            protein_node=resolve_protein_node(cands, gin_is_member),
            promoter_node=cands["promoter"] if gin_is_member(cands["promoter"]) else None,
        )
        if rec.walk_anchor_nodes:
            records.append(rec)
            continue

        # Optional legacy fallback: treat unmapped gene-channel token as metabolite name.
        if allow_gene_metabolite_fallback or seed_kinds is None:
            node = resolve_metabolite_gin_node(
                symbol,
                gin_is_member,
                pool,
                metabolite_index=index,
                pubchem=resolver,
                allow_online=allow_online_metabolite,
            )
            if node:
                records.append(
                    CentralDogmaSeed(symbol=symbol, kind="metabolite", metabolite_node=node)
                )
                continue

        records.append(rec)
        unmapped.append(gene_sym)

    qc = SeedMappingQC(
        n_input=len(union_inputs),
        n_with_rna=sum(1 for r in records if r.rna_node),
        n_with_protein=sum(1 for r in records if r.protein_node),
        n_with_promoter=sum(1 for r in records if r.promoter_node),
        n_triplet_complete=sum(1 for r in records if r.kind == "gene" and r.triplet_complete),
        n_unmapped=len(unmapped),
        unmapped_symbols=unmapped[:200],
        n_metabolite=sum(1 for r in records if r.kind == "metabolite" and r.metabolite_node),
        n_gene=sum(1 for r in records if r.kind == "gene"),
    )
    return records, qc


def collect_walk_input_nodes(records: list[CentralDogmaSeed]) -> set[str]:
    nodes: set[str] = set()
    for rec in records:
        nodes.update(rec.walk_anchor_nodes)
    return nodes


def collect_regulatory_seed_nodes(records: list[CentralDogmaSeed]) -> set[str]:
    return collect_walk_input_nodes(records)

"""Build seed manifests from gene lists, metabolite lists, or AnnData (scRNA mode)."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from gin_tool.gin_defaults import DEFAULT_SYMBOL_BLACKLIST
from gin_tool.mapping import load_symbol_blacklist
from gin_tool.subgraph_pipeline.manifest import SeedManifest


def _normalize_gene_symbol(name: str) -> str:
    return str(name).split(".")[0].strip().upper()


def _read_lines(path: str | Path) -> list[str]:
    return Path(path).read_text(encoding="utf-8").splitlines()


def _iter_nonempty(lines: Iterable[str]) -> list[str]:
    out: list[str] = []
    for raw in lines:
        line = str(raw).split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def build_seed_manifest_from_inputs(
    *,
    genes: Iterable[str] | None = None,
    metabolites: Iterable[str] | None = None,
    run_id: str = "seed_inputs",
    blacklist_path: str | Path | None = None,
) -> SeedManifest:
    """
    Build a seed manifest with **explicit** gene vs metabolite channels.

    Genes are uppercased HGNC-like symbols (blacklist applied).
    Metabolites keep user casing / id form (CAS / CHEBI / CID / SID / …);
    they are never run through the gene blacklist or uppercasing.
    """
    bl_path = blacklist_path if blacklist_path is not None else DEFAULT_SYMBOL_BLACKLIST
    blacklist = load_symbol_blacklist(bl_path)

    seen: set[str] = set()
    union: list[str] = []
    seed_kinds: dict[str, str] = {}
    blocked: list[str] = []

    for raw in _iter_nonempty(genes or []):
        sym = _normalize_gene_symbol(raw)
        if not sym:
            continue
        if sym in blacklist:
            blocked.append(sym)
            continue
        if sym not in seen:
            seen.add(sym)
            union.append(sym)
            seed_kinds[sym] = "gene"

    for raw in _iter_nonempty(metabolites or []):
        mid = raw.strip()
        if not mid:
            continue
        # Prefer original token; avoid colliding with an already-added gene symbol.
        if mid in seen:
            # Same string already present as gene — keep gene unless user only sent met.
            continue
        key = mid
        if key not in seen:
            seen.add(key)
            union.append(key)
            seed_kinds[key] = "metabolite"

    if not union:
        raise ValueError("no seed inputs remain after filtering")

    n_gene = sum(1 for k in seed_kinds.values() if k == "gene")
    n_met = sum(1 for k in seed_kinds.values() if k == "metabolite")
    return SeedManifest(
        run_id=run_id,
        source="seed_inputs",
        union_genes=union,
        seed_kinds=seed_kinds,
        stats={
            "n_input": n_gene + n_met + len(blocked),
            "n_after_blacklist": len(union),
            "n_blacklisted": len(blocked),
            "n_gene": n_gene,
            "n_metabolite": n_met,
            "blacklisted_sample": blocked[:30],
            "blacklist_path": str(bl_path),
            "channels": {"genes": n_gene, "metabolites": n_met},
        },
    )


def build_seed_manifest_from_gene_list(
    genes: Iterable[str],
    *,
    metabolites: Iterable[str] | None = None,
    run_id: str = "gene_list",
    blacklist_path: str | Path | None = None,
    node_list_dir: str | Path | None = None,  # kept for API compat; unused
) -> SeedManifest:
    """
    Gene-channel constructor (optionally with a separate metabolite iterable).

    Metabolites must be passed via ``metabolites=`` — they are not inferred from
    the gene list (avoids CASP12↔metabolite-pool collisions).
    """
    del node_list_dir  # no longer used for heuristic typing
    manifest = build_seed_manifest_from_inputs(
        genes=genes,
        metabolites=metabolites,
        run_id=run_id,
        blacklist_path=blacklist_path,
    )
    # Preserve legacy source label when only genes are supplied.
    if not metabolites:
        manifest.source = "gene_list"
    return manifest


def build_seed_manifest_from_file(
    path: str | Path,
    *,
    metabolite_file: str | Path | None = None,
    run_id: str | None = None,
    blacklist_path: str | Path | None = None,
) -> SeedManifest:
    """Load gene symbols from ``path``; optional metabolite ids from ``metabolite_file``."""
    p = Path(path)
    rid = run_id or p.stem
    mets = _read_lines(metabolite_file) if metabolite_file else None
    return build_seed_manifest_from_gene_list(
        _read_lines(p),
        metabolites=mets,
        run_id=rid,
        blacklist_path=blacklist_path,
    )


def build_seed_manifest_from_files(
    *,
    gene_file: str | Path | None = None,
    metabolite_file: str | Path | None = None,
    run_id: str | None = None,
    blacklist_path: str | Path | None = None,
) -> SeedManifest:
    """Explicit dual-file constructor; at least one file required."""
    if not gene_file and not metabolite_file:
        raise ValueError("need gene_file and/or metabolite_file")
    genes = _read_lines(gene_file) if gene_file else []
    mets = _read_lines(metabolite_file) if metabolite_file else []
    rid = run_id
    if rid is None:
        if gene_file:
            rid = Path(gene_file).stem
        else:
            rid = Path(metabolite_file).stem  # type: ignore[arg-type]
    return build_seed_manifest_from_inputs(
        genes=genes,
        metabolites=mets,
        run_id=rid,
        blacklist_path=blacklist_path,
    )

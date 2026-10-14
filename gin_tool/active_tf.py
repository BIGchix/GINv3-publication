"""Infer active transcription factors for context-aware subgraph expansion."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

from gin_tool.enrichment import enrich_many_terms
from gin_tool.io import iter_sifi_edges
from gin_tool.mapping import map_ids_to_gin_nodes
from gin_tool.node_lists import load_node_list_dir


def bh_fdr(pvalues: Sequence[float]) -> list[float]:
    """Benjamini–Hochberg adjusted p-values (returns list aligned to input)."""
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvalues[i])
    adj = [1.0] * m
    prev = 1.0
    for rank_rev, idx in enumerate(reversed(order)):
        rank = m - rank_rev
        p = pvalues[idx]
        q = min(prev, p * m / rank)
        adj[idx] = q
        prev = q
    return adj


def gin_node_to_gene_symbol(name: str) -> str | None:
    """
    Map a GIN entity node to an HGNC-like gene symbol when possible.

    Returns None for virtual transition nodes and unrecognized patterns.
    """
    if ";" in name and ".promoter" not in name.lower():
        return None
    lower = name.lower()
    if lower.endswith(".rna"):
        return name.rsplit(".", 1)[0]
    if lower.endswith(".protein"):
        return name.rsplit(".", 1)[0]
    if ".promoter" in lower:
        return name.split(".promoter", 1)[0]
    if lower.endswith(".state"):
        return name[: -len(".state")]
    if "." not in name:
        return name
    return None


def promoter_target_to_gene(target: str) -> str | None:
    """Extract gene symbol from a TF regulatory target such as ``TP53.promoter;state``."""
    if ".promoter" not in target.lower():
        return None
    gene = target.split(".promoter", 1)[0]
    if not gene or ";" in gene:
        return None
    return gene


def load_tf_gin_node_map(node_list_dir: Path) -> dict[str, str]:
    """Map bare HGNC symbol -> qualified GIN TF node name (e.g. ``TP53.protein``)."""
    pools = load_node_list_dir(node_list_dir)
    out: dict[str, str] = {}
    for n in pools.get("tf", []):
        n = n.strip()
        if not n:
            continue
        sym = gin_node_to_gene_symbol(n) or n
        out[sym] = n
    return out


def load_tf_symbol_pool(node_list_dir: Path) -> set[str]:
    """TF gene symbols (bare HGNC) from ``nodeList_TF.txt``."""
    return set(load_tf_gin_node_map(node_list_dir).keys())


def parse_id_list(spec: str | None, *, path: Path | None = None) -> list[str]:
    """Comma-separated IDs or one ID per line from an optional file."""
    out: list[str] = []
    if path is not None:
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line)
    if spec:
        out.extend(s.strip() for s in spec.split(",") if s.strip())
    seen: set[str] = set()
    unique: list[str] = []
    for x in out:
        if x not in seen:
            seen.add(x)
            unique.append(x)
    return unique


def build_tf_regulon_map(
    sifi_path: str,
    tf_symbols: set[str],
    *,
    tf_gin_by_bare: Mapping[str, str] | None = None,
    max_lines: int | None = None,
) -> tuple[dict[str, set[str]], set[str]]:
    """
    Stream ``.sifi`` and collect TF → target gene symbols.

    Returns (regulon_map, universe) keyed by **bare** TF symbols.
    """
    regulons: dict[str, set[str]] = {tf: set() for tf in tf_symbols}
    source_to_bare: dict[str, str] = {}
    if tf_gin_by_bare:
        for bare, gin in tf_gin_by_bare.items():
            if bare in tf_symbols:
                source_to_bare[gin] = bare
                source_to_bare[bare] = bare
    else:
        for tf in tf_symbols:
            source_to_bare[tf] = tf

    universe: set[str] = set()
    for source, target, _w in iter_sifi_edges(sifi_path, max_lines=max_lines):
        bare = source_to_bare.get(source)
        if bare is None:
            bare_sym = gin_node_to_gene_symbol(source)
            if bare_sym and bare_sym in regulons:
                bare = bare_sym
            elif source in regulons:
                bare = source
            else:
                continue
        if bare not in regulons:
            continue
        gene = promoter_target_to_gene(target)
        if not gene:
            continue
        regulons[bare].add(gene)
        universe.add(gene)
    return regulons, universe


def seed_nodes_to_gene_symbols(seed_gin_nodes: Iterable[str]) -> set[str]:
    symbols: set[str] = set()
    for node in seed_gin_nodes:
        sym = gin_node_to_gene_symbol(node)
        if sym:
            symbols.add(sym)
    return symbols


@dataclass
class ActiveTFRecord:
    gin_node: str
    symbol: str
    overlap_k: int
    overlap_genes: list[str]
    pvalue: float
    p_adj: float
    odds_ratio: float
    regulon_size: int
    method: str
    evidence_tier: str = "moderate"

    def to_dict(self) -> dict:
        return {
            "gin_node": self.gin_node,
            "symbol": self.symbol,
            "overlap_k": self.overlap_k,
            "overlap_genes": self.overlap_genes,
            "pvalue": self.pvalue,
            "p_adj": self.p_adj,
            "odds_ratio": self.odds_ratio,
            "regulon_size": self.regulon_size,
            "method": self.method,
            "evidence_tier": self.evidence_tier,
        }


@dataclass
class ActiveTFRequest:
    """Inputs for active-TF inference."""

    seed_gin_nodes: Sequence[str]
    sifi_path: str
    node_list_dir: Path
    user_active_tfs: Sequence[str] | None = None
    gin_node_is_member: Callable[[str], bool] | None = None
    top_n: int = 30
    fdr_alpha: float = 0.05
    min_overlap: int = 2
    low_power_threshold: int = 5
    max_lines: int | None = None
    # Pipeline default is rwr_cohesion (B); rwr_rewire (C) optional; fisher=legacy.
    method: str = "rwr_cohesion"
    rwr_work_dir: Path | None = None
    rwr_n_null: int = 1000
    rwr_n_jobs: int = 8
    rwr_emp_p_cutoff: float = 0.05
    rwr_rewire_sweeps: float = 0.1
    rwr_real_hops: int = 1
    rwr_seed: int = 42


@dataclass
class ActiveTFResult:
    records: list[ActiveTFRecord]
    method: str
    seed_symbols: set[str]
    universe_size: int
    tf_pool_size: int
    warnings: list[str] = field(default_factory=list)
    skipped: bool = False
    skip_reason: str | None = None

    def active_gin_nodes(self) -> list[str]:
        return [r.gin_node for r in self.records]

    def active_symbols(self) -> list[str]:
        return [r.symbol for r in self.records]

    def to_qc_dict(self) -> dict:
        return {
            "method": self.method,
            "skipped": self.skipped,
            "skip_reason": self.skip_reason,
            "seed_symbol_count": len(self.seed_symbols),
            "seed_symbols_sample": sorted(self.seed_symbols)[:50],
            "universe_size": self.universe_size,
            "tf_pool_size": self.tf_pool_size,
            "active_tf_count": len(self.records),
            "warnings": self.warnings,
            "active_tfs": [r.to_dict() for r in self.records],
        }


def _map_user_tfs(
    user_ids: Sequence[str],
    gin_node_is_member: Callable[[str], bool],
    tf_pool: set[str],
) -> tuple[list[ActiveTFRecord], list[str]]:
    mapping, report = map_ids_to_gin_nodes(
        user_ids,
        gin_node_is_member,
        missing_strategy="warn",
    )
    warnings = list(report.warnings)
    records: list[ActiveTFRecord] = []
    for uid in user_ids:
        uid = uid.strip()
        if not uid:
            continue
        gnode = mapping.get(uid)
        if not gnode:
            warnings.append(f"active TF not mapped: {uid!r}")
            continue
        sym = gin_node_to_gene_symbol(gnode) or uid
        records.append(
            ActiveTFRecord(
                gin_node=gnode,
                symbol=sym,
                overlap_k=0,
                overlap_genes=[],
                pvalue=0.0,
                p_adj=0.0,
                odds_ratio=0.0,
                regulon_size=0,
                method="user_supplied",
                evidence_tier="user",
            )
        )
    _ = tf_pool  # reserved for future TF-pool validation warnings
    return records, warnings


def _regulon_enrich(
    seed_symbols: set[str],
    regulons: Mapping[str, set[str]],
    universe: set[str],
    tf_pool: set[str],
    *,
    top_n: int,
    fdr_alpha: float,
    min_overlap: int,
    low_power: bool,
) -> tuple[list[ActiveTFRecord], list[str]]:
    warnings: list[str] = []
    if not seed_symbols:
        return [], ["empty seed gene symbol set"]
    if not universe:
        return [], ["empty TF target universe from GIN"]

    term_to_genes = {tf: genes for tf, genes in regulons.items() if genes}
    rows = enrich_many_terms(seed_symbols, term_to_genes, universe)
    if not rows:
        if low_power:
            return _connectivity_fallback(seed_symbols, regulons, tf_pool), [
                "low_power: |S|<threshold; using connectivity fallback (overlap>=1)"
            ]
        return [], ["no TF regulon overlapped seeds"]

    p_adj = bh_fdr([r.pvalue for r in rows])
    candidates: list[ActiveTFRecord] = []
    for row, q in zip(rows, p_adj):
        if row.k < min_overlap and not low_power:
            continue
        overlap_genes = sorted(seed_symbols & term_to_genes[row.term_id])
        candidates.append(
            ActiveTFRecord(
                gin_node=row.term_id,
                symbol=row.term_id,
                overlap_k=row.k,
                overlap_genes=overlap_genes,
                pvalue=row.pvalue,
                p_adj=q,
                odds_ratio=row.odds_ratio,
                regulon_size=row.n,
                method="regulon_enrich",
                evidence_tier="low" if low_power else "moderate",
            )
        )

    # Force-include seed TFs present in pool
    forced = seed_symbols & tf_pool
    by_sym = {r.symbol: r for r in candidates}
    for sym in sorted(forced):
        if sym in by_sym:
            by_sym[sym].method = "regulon_enrich_forced_seed_tf"
            continue
        overlap_genes = sorted(seed_symbols & regulons.get(sym, set()))
        by_sym[sym] = ActiveTFRecord(
            gin_node=sym,
            symbol=sym,
            overlap_k=len(overlap_genes),
            overlap_genes=overlap_genes,
            pvalue=0.0,
            p_adj=0.0,
            odds_ratio=0.0,
            regulon_size=len(regulons.get(sym, set())),
            method="regulon_enrich_forced_seed_tf",
            evidence_tier="moderate",
        )
    candidates = list(by_sym.values())

    if low_power:
        candidates.sort(key=lambda r: (-r.overlap_k, r.pvalue))
        return candidates[:top_n], warnings

    selected = [r for r in candidates if r.p_adj <= fdr_alpha]
    if not selected:
        selected = sorted(candidates, key=lambda r: r.pvalue)[:top_n]
        warnings.append(
            f"no TF passed FDR<={fdr_alpha}; returning top {min(top_n, len(selected))} by p-value"
        )
    else:
        selected.sort(key=lambda r: (r.p_adj, -r.overlap_k))
        selected = selected[:top_n]
    return selected, warnings


def _connectivity_fallback(
    seed_symbols: set[str],
    regulons: Mapping[str, set[str]],
    tf_pool: set[str],
) -> list[ActiveTFRecord]:
    out: list[ActiveTFRecord] = []
    for tf in sorted(tf_pool):
        overlap = seed_symbols & regulons.get(tf, set())
        if not overlap:
            continue
        out.append(
            ActiveTFRecord(
                gin_node=tf,
                symbol=tf,
                overlap_k=len(overlap),
                overlap_genes=sorted(overlap),
                pvalue=1.0,
                p_adj=1.0,
                odds_ratio=0.0,
                regulon_size=len(regulons.get(tf, set())),
                method="connectivity_fallback",
                evidence_tier="low",
            )
        )
    out.sort(key=lambda r: (-r.overlap_k, r.symbol))
    return out


def _records_from_rwr(
    selected: list[dict],
    tf_gin_map: Mapping[str, str],
    *,
    method_tag: str,
) -> list[ActiveTFRecord]:
    out: list[ActiveTFRecord] = []
    for row in selected:
        sym = str(row["symbol"])
        gin = tf_gin_map.get(sym, sym)
        emp = float(row.get("rwr_emp_p", 1.0))
        q = float(row.get("rwr_emp_q", emp))
        score = float(row.get("rwr_score", 0.0))
        out.append(
            ActiveTFRecord(
                gin_node=gin,
                symbol=sym,
                overlap_k=0,
                overlap_genes=[],
                pvalue=emp,
                p_adj=q,
                odds_ratio=score,
                regulon_size=0,
                method=method_tag,
                evidence_tier="rwr_null",
            )
        )
    return out


def infer_active_tfs(request: ActiveTFRequest) -> ActiveTFResult:
    """
    Infer or accept active TF nodes for downstream subgraph filtering.

    If ``user_active_tfs`` is non-empty, maps user IDs and returns without ranking.
    Else ``method``:
      - ``rwr_cohesion`` (B, default): directed RWR + cohesion null
      - ``rwr_rewire`` (C): directed RWR + edge-rewire null (slower, stricter)
      - ``fisher``: legacy regulon Fisher enrich
    """
    tf_gin_map = load_tf_gin_node_map(request.node_list_dir)
    tf_pool = set(tf_gin_map.keys())
    seed_symbols = seed_nodes_to_gene_symbols(request.seed_gin_nodes)

    if request.user_active_tfs:
        if request.gin_node_is_member is None:
            raise ValueError("gin_node_is_member required when user_active_tfs is set")
        records, warnings = _map_user_tfs(
            request.user_active_tfs,
            request.gin_node_is_member,
            tf_pool,
        )
        return ActiveTFResult(
            records=records,
            method="user_supplied",
            seed_symbols=seed_symbols,
            universe_size=0,
            tf_pool_size=len(tf_pool),
            warnings=warnings,
        )

    if not seed_symbols:
        return ActiveTFResult(
            records=[],
            method="skipped",
            seed_symbols=set(),
            universe_size=0,
            tf_pool_size=len(tf_pool),
            warnings=[],
            skipped=True,
            skip_reason="no_gene_symbols_in_seeds",
        )

    method = (request.method or "rwr_cohesion").strip().lower()
    if method in ("rwr_cohesion", "cohesion", "b", "rwr_rewire", "rewire", "c"):
        from gin_tool.active_tf_rwr import resolve_null_strength, run_rwr_null_active_tfs

        strength = resolve_null_strength(method)
        method_tag = f"rwr_{strength}"
        work = Path(request.rwr_work_dir or Path("gin_out") / "active_tf_rwr_work")
        selected, _summary = run_rwr_null_active_tfs(
            sorted(seed_symbols),
            method=method,
            work_dir=work,
            n_rwr_null=request.rwr_n_null,
            n_jobs=request.rwr_n_jobs,
            emp_p_cutoff=request.rwr_emp_p_cutoff,
            rewire_sweeps=request.rwr_rewire_sweeps,
            real_hops=request.rwr_real_hops,
            seed=request.rwr_seed,
        )
        records = _records_from_rwr(selected, tf_gin_map, method_tag=method_tag)
        warnings: list[str] = []
        if not records:
            warnings.append(
                f"no active TFs with rwr_emp_p<={request.rwr_emp_p_cutoff} ({method_tag})"
            )
        return ActiveTFResult(
            records=records,
            method=method_tag,
            seed_symbols=seed_symbols,
            universe_size=0,
            tf_pool_size=len(tf_pool),
            warnings=warnings,
        )

    if method not in ("fisher", "regulon_enrich"):
        raise ValueError(
            f"unknown active-TF method {request.method!r}; "
            "use rwr_cohesion, rwr_rewire, fisher, or pass user_active_tfs"
        )

    regulons, universe = build_tf_regulon_map(
        request.sifi_path,
        tf_pool,
        tf_gin_by_bare=tf_gin_map,
        max_lines=request.max_lines,
    )
    low_power = len(seed_symbols) < request.low_power_threshold
    records, warnings = _regulon_enrich(
        seed_symbols,
        regulons,
        universe,
        tf_pool,
        top_n=request.top_n,
        fdr_alpha=request.fdr_alpha,
        min_overlap=request.min_overlap,
        low_power=low_power,
    )
    if low_power:
        warnings.insert(
            0,
            f"|S|={len(seed_symbols)} < {request.low_power_threshold}; enrichment downgraded",
        )
    out_method = (
        "connectivity_fallback"
        if records and all(r.method == "connectivity_fallback" for r in records)
        else "regulon_enrich"
    )
    if not records:
        warnings.append("no active TFs selected")
    return ActiveTFResult(
        records=records,
        method=out_method,
        seed_symbols=seed_symbols,
        universe_size=len(universe),
        tf_pool_size=len(tf_pool),
        warnings=warnings,
    )


def write_active_tfs_tsv(path: str | Path, result: ActiveTFResult) -> None:
    path = Path(path)
    header = (
        "gin_node\tsymbol\toverlap_k\tpvalue\tp_adj\todds_ratio\t"
        "regulon_size\tmethod\tevidence_tier\toverlap_genes\n"
    )
    lines = [header]
    for r in result.records:
        genes = ";".join(r.overlap_genes)
        lines.append(
            f"{r.gin_node}\t{r.symbol}\t{r.overlap_k}\t{r.pvalue:.6g}\t{r.p_adj:.6g}\t"
            f"{r.odds_ratio:.6g}\t{r.regulon_size}\t{r.method}\t{r.evidence_tier}\t{genes}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


def write_active_tf_qc_json(path: str | Path, result: ActiveTFResult) -> None:
    Path(path).write_text(
        json.dumps(result.to_qc_dict(), indent=2),
        encoding="utf-8",
    )

"""GRN allowlist prune: TF→promoter;state edges filtered by ChIP evidence.

Rules (aligned with ``p53_cell_grn_prune_v4_seedhist`` / ``wnt_cell_grn_prune_v4_seedhist``):

- **Unmeasured TF** in a cell: keep all TF→``promoter;state`` edges (atlas incomplete ≠ inactive).
- **Measured TF**: keep edge iff ``(TF, gene)`` is in ``allowlist_union.tsv``.
- **Pathway seed pairs** (optional, global): when both TF and target gene are pathway
  seed symbols, keep the edge even for measured TFs (curated pathway transcription).
- All other edges are untouched. Degree-0 nodes are dropped.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from gin_tool.active_tf import load_tf_symbol_pool
from gin_tool.subgraph import write_subgraph_tsv
from gin_tool.subgraph_pipeline.met_source_seed import load_subgraph_edges_tsv


def _gene_from_promoter_state(tgt: str) -> str | None:
    if tgt.endswith(".promoter;state"):
        return tgt.split(".", 1)[0]
    return None


def _tf_gene_from_protein(src: str) -> str | None:
    if src.endswith(".protein") and ";" not in src:
        return src[: -len(".protein")]
    return None


def is_tf_to_promoter_state_edge(
    src: str, tgt: str, *, tf_symbols: set[str]
) -> bool:
    gene = _tf_gene_from_protein(src)
    return gene is not None and gene in tf_symbols and _gene_from_promoter_state(tgt) is not None


def load_allowlist_pairs(path: str | Path) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    with Path(path).open(encoding="utf-8") as fh:
        header: list[str] | None = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [c.lower() for c in parts]
                continue
            row = dict(zip(header, parts))
            tf = row.get("tf", "").strip()
            gene = row.get("gene", "").strip()
            if tf and gene:
                pairs.add((tf, gene))
    return pairs


def load_measured_tfs(path: str | Path) -> set[str]:
    out: set[str] = set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.add(line)
    return out


def drop_degree_zero_nodes(
    edges: list[tuple[str, str, float]],
) -> tuple[list[tuple[str, str, float]], int]:
    if not edges:
        return [], 0
    nodes: set[str] = set()
    for src, tgt, _ in edges:
        nodes.add(src)
        nodes.add(tgt)
    incident: dict[str, int] = {n: 0 for n in nodes}
    for src, tgt, _ in edges:
        incident[src] += 1
        incident[tgt] += 1
    dead = {n for n, d in incident.items() if d == 0}
    if not dead:
        return edges, 0
    kept = [(s, t, w) for s, t, w in edges if s not in dead and t not in dead]
    return kept, len(dead)


@dataclass
class GrnAllowlistPruneResult:
    cell: str
    source_edges: str
    allowlist: str
    peaks_dir: str
    edges: list[tuple[str, str, float]]
    n_edges_in: int
    n_edges_out: int
    n_edges_removed_total: int
    n_tf_to_state_in: int = 0
    n_tf_to_state_kept: int = 0
    n_tf_to_state_removed: int = 0
    n_tf_to_state_spared_unmeasured: int = 0
    n_tf_to_state_kept_allowlisted: int = 0
    n_tf_to_state_kept_seed_pair: int = 0
    n_measured_tfs: int = 0
    n_spared_tfs: int = 0
    n_pruned_tfs: int = 0
    spared_tfs_sample: list[str] = field(default_factory=list)
    pruned_tfs_sample: list[str] = field(default_factory=list)
    n_nodes_in: int = 0
    n_nodes_out: int = 0
    n_nodes_removed: int = 0
    n_allowlist_pairs: int = 0
    removed_sample: list[str] = field(default_factory=list)
    measured_tfs: list[str] = field(default_factory=list)
    seed_genes: list[str] = field(default_factory=list)
    keep_seed_tf_to_seed_gene: bool = False

    def to_qc_dict(self) -> dict:
        rules = {
            "target": "TF.protein → GENE.promoter;state (TF ∈ nodeList_TF)",
            "if_tf_unmeasured_in_cell": (
                "keep all original TF→state edges (atlas incomplete ≠ inactive)"
            ),
            "if_tf_measured_in_cell": "keep edge iff (TF, gene) ∈ allowlist_union",
            "measured_definition": (
                "TF label in cell peak BEDs (GTRD/ReMap/Atlas), mapped to nodeList_TF"
            ),
            "untouched": "non TF→promoter;state edges",
            "no_new_edges": True,
            "drop_degree0_nodes": True,
        }
        if self.keep_seed_tf_to_seed_gene:
            rules["keep_seed_tf_to_seed_gene"] = (
                "measured TF→seed-gene promoter;state kept when both ends are pathway seeds"
            )
        return {
            "cell": self.cell,
            "source_edges": self.source_edges,
            "allowlist": self.allowlist,
            "peaks_dir": self.peaks_dir,
            "n_edges_in": self.n_edges_in,
            "n_edges_out": self.n_edges_out,
            "n_edges_removed_total": self.n_edges_removed_total,
            "n_tf_to_state_in": self.n_tf_to_state_in,
            "n_tf_to_state_kept": self.n_tf_to_state_kept,
            "n_tf_to_state_removed": self.n_tf_to_state_removed,
            "n_tf_to_state_spared_unmeasured": self.n_tf_to_state_spared_unmeasured,
            "n_tf_to_state_kept_allowlisted": self.n_tf_to_state_kept_allowlisted,
            "n_tf_to_state_kept_seed_pair": self.n_tf_to_state_kept_seed_pair,
            "n_measured_tfs": self.n_measured_tfs,
            "n_spared_tfs": self.n_spared_tfs,
            "n_pruned_tfs": self.n_pruned_tfs,
            "spared_tfs_sample": self.spared_tfs_sample,
            "pruned_tfs_sample": self.pruned_tfs_sample,
            "n_nodes_in": self.n_nodes_in,
            "n_nodes_out": self.n_nodes_out,
            "n_nodes_removed": self.n_nodes_removed,
            "n_allowlist_pairs": self.n_allowlist_pairs,
            "removed_sample": self.removed_sample,
            "rules": rules,
            "seed_genes": self.seed_genes[:20] if len(self.seed_genes) > 20 else self.seed_genes,
            "n_seed_genes": len(self.seed_genes),
            "keep_seed_tf_to_seed_gene": self.keep_seed_tf_to_seed_gene,
        }


def prune_grn_allowlist(
    edges: list[tuple[str, str, float]],
    *,
    cell: str,
    source_edges: str,
    allowlist_path: str | Path,
    peaks_dir: str | Path,
    measured_tfs: set[str],
    tf_symbols: set[str],
    seed_genes: set[str] | None = None,
    keep_seed_tf_to_seed_gene: bool = False,
) -> GrnAllowlistPruneResult:
    """Filter TF→promoter;state edges; return pruned edge list + QC stats."""
    allowlist = load_allowlist_pairs(allowlist_path)
    seed_genes = seed_genes or set()
    nodes_in = {n for s, t, _ in edges for n in (s, t)}

    kept: list[tuple[str, str, float]] = []
    removed_tf_state: list[tuple[str, str, float]] = []
    n_tf_in = 0
    n_kept = 0
    n_removed = 0
    n_spared_unmeasured = 0
    n_kept_allow = 0
    n_kept_seed = 0
    tf_had_removed: dict[str, bool] = {}
    tf_had_kept: dict[str, bool] = {}

    for src, tgt, w in edges:
        if not is_tf_to_promoter_state_edge(src, tgt, tf_symbols=tf_symbols):
            kept.append((src, tgt, w))
            continue
        n_tf_in += 1
        tf = _tf_gene_from_protein(src)
        gene = _gene_from_promoter_state(tgt)
        assert tf is not None and gene is not None

        keep = False
        if tf not in measured_tfs:
            keep = True
            n_spared_unmeasured += 1
        elif (tf, gene) in allowlist:
            keep = True
            n_kept_allow += 1
        elif keep_seed_tf_to_seed_gene and tf in seed_genes and gene in seed_genes:
            keep = True
            n_kept_seed += 1

        if keep:
            kept.append((src, tgt, w))
            n_kept += 1
            tf_had_kept[tf] = True
        else:
            removed_tf_state.append((src, tgt, w))
            n_removed += 1
            tf_had_removed[tf] = True

    spared_tfs = sorted(tf for tf in tf_symbols if tf not in measured_tfs)
    pruned_tfs = sorted(tf for tf in measured_tfs if tf_had_removed.get(tf))

    kept2, _n_dead = drop_degree_zero_nodes(kept)
    nodes_out = {n for s, t, _ in kept2 for n in (s, t)}

    removed_sample = [
        f"{s}\t{t}\t{w:g}" for s, t, w in removed_tf_state[:30]
    ]

    return GrnAllowlistPruneResult(
        cell=cell,
        source_edges=source_edges,
        allowlist=str(allowlist_path),
        peaks_dir=str(peaks_dir),
        edges=kept2,
        n_edges_in=len(edges),
        n_edges_out=len(kept2),
        n_edges_removed_total=len(edges) - len(kept2),
        n_tf_to_state_in=n_tf_in,
        n_tf_to_state_kept=n_kept,
        n_tf_to_state_removed=n_removed,
        n_tf_to_state_spared_unmeasured=n_spared_unmeasured,
        n_tf_to_state_kept_allowlisted=n_kept_allow,
        n_tf_to_state_kept_seed_pair=n_kept_seed,
        n_measured_tfs=len(measured_tfs & tf_symbols),
        n_spared_tfs=len(spared_tfs),
        n_pruned_tfs=len(pruned_tfs),
        spared_tfs_sample=spared_tfs[:40],
        pruned_tfs_sample=pruned_tfs[:40],
        n_nodes_in=len(nodes_in),
        n_nodes_out=len(nodes_out),
        n_nodes_removed=len(nodes_in) - len(nodes_out),
        n_allowlist_pairs=len(allowlist),
        removed_sample=removed_sample,
        measured_tfs=sorted(measured_tfs),
        seed_genes=sorted(seed_genes),
        keep_seed_tf_to_seed_gene=keep_seed_tf_to_seed_gene,
    )


def write_pruned_cell_dir(
    out_dir: str | Path,
    result: GrnAllowlistPruneResult,
) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_subgraph_tsv(str(out / "subgraph_edges.tsv"), result.edges)
    (out / "prune_qc.json").write_text(
        json.dumps(result.to_qc_dict(), indent=2) + "\n", encoding="utf-8"
    )
    (out / "measured_tfs.txt").write_text(
        "\n".join(result.measured_tfs) + ("\n" if result.measured_tfs else ""),
        encoding="utf-8",
    )
    return out


def default_chip_align_root() -> Path:
    return Path(__file__).resolve().parents[2] / "runtime" / "chip_gin_align"


def chip_align_cell_dir(chip_root: str | Path, cell: str) -> Path:
    return Path(chip_root) / cell


def load_seed_genes(seed_file: str | Path) -> set[str]:
    genes: set[str] = set()
    for line in Path(seed_file).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.add(line)
    return genes

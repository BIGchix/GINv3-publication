"""JSON-serializable manifests for subgraph pipeline seed genes."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class SeedManifest:
    """Seed gene set for subgraph extraction."""

    run_id: str
    union_genes: list[str]
    source: str = "gene_list"  # gene_list | scrna
    cell_type_key: str = ""
    n_hvg: int = 0
    n_markers_per_type: int = 0
    min_cells_for_markers: int = 0
    hvg_genes: list[str] = field(default_factory=list)
    markers_by_cell_type: dict[str, list[str]] = field(default_factory=dict)
    marker_genes_flat: list[str] = field(default_factory=list)
    overlap_hvg_markers: list[str] = field(default_factory=list)
    stats: dict[str, Any] = field(default_factory=dict)
    gin_anchor_nodes: list[str] = field(default_factory=list)
    seed_kinds: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

    @classmethod
    def load(cls, path: str | Path) -> SeedManifest:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(**data)


@dataclass
class NodeCountStats:
    n_real_nodes: int = 0
    n_virtual_nodes: int = 0
    n_gin_edges: int = 0
    n_compressed_edges: int = 0
    n_seeds_no_loop: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ConnectivityQC:
    n_components_before: int = 0
    n_components_after: int = 0
    largest_component_frac: float = 0.0
    nodes_removed_l1: int = 0
    bridge_edges_added: int = 0
    seeds_in_real_nodes: int = 0
    seeds_total: int = 0
    node_stats: NodeCountStats = field(default_factory=NodeCountStats)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        return d


@dataclass
class SeedSelfLoopQC:
    n_seeds: int = 0
    n_triplet_complete: int = 0
    n_rna_self_loop: int = 0
    n_protein_self_loop: int = 0
    n_promoter_self_loop: int = 0
    n_either_self_loop: int = 0
    n_failed: int = 0
    max_real_hops: int = 0
    failed_symbols_sample: list[str] = field(default_factory=list)
    node_stats: NodeCountStats = field(default_factory=NodeCountStats)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class EdgeTopologyQC:
    n_edges: int = 0
    canonical_forward: int = 0
    cross_layer_forward: int = 0
    reverse: int = 0
    lateral: int = 0
    involving_other: int = 0
    by_layer_pair: dict[str, int] = field(default_factory=dict)
    modeling: dict[str, Any] = field(default_factory=dict)
    node_stats: NodeCountStats = field(default_factory=NodeCountStats)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CentralDogmaMappingQC:
    n_input: int = 0
    n_with_rna: int = 0
    n_with_protein: int = 0
    n_with_promoter: int = 0
    n_triplet_complete: int = 0
    n_unmapped: int = 0
    unmapped_symbols: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SubgraphExtractResult:
    run_id: str
    seed_manifest: SeedManifest
    walk: Any
    active_tf_symbols: list[str]
    gin_seed_nodes: list[str]
    central_dogma_qc: CentralDogmaMappingQC
    seed_self_loop_qc: SeedSelfLoopQC
    connectivity_qc: ConnectivityQC
    edge_topology_qc: EdgeTopologyQC
    prefilter_stats: Any | None = None
    warnings: list[str] = field(default_factory=list)

    def node_stats(self) -> NodeCountStats:
        from gin_tool.subgraph_pipeline.stats import node_count_stats

        n_no_loop = self.seed_self_loop_qc.n_failed
        return node_count_stats(self.walk, n_seeds_no_loop=n_no_loop)

    def save_summary(self, path: str | Path) -> None:
        summary = {
            "run_id": self.run_id,
            "node_stats_final": self.node_stats().to_dict(),
            "central_dogma": self.central_dogma_qc.to_dict(),
            "seed_self_loop": self.seed_self_loop_qc.to_dict(),
            "connectivity": self.connectivity_qc.to_dict(),
            "edge_topology": self.edge_topology_qc.to_dict(),
            "active_tf_count": len(self.active_tf_symbols),
            "warnings": self.warnings,
        }
        Path(path).write_text(
            json.dumps(summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

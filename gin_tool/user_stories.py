"""Core user stories for the GIN-based multi-omics tool (product scope)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import FrozenSet


@dataclass(frozen=True)
class UserStory:
    """One end-to-end workflow the tool optimizes for."""

    id: str
    title: str
    inputs: FrozenSet[str]
    outputs: FrozenSet[str]
    description: str


# US1: differential / candidate lists -> mechanistic neighborhood
US1 = UserStory(
    id="US1",
    title="From omics hits to GIN neighborhood and report",
    inputs=frozenset(
        {
            "gene_or_entity_ids",  # symbols or GIN node suffixes, e.g. TP53.protein
            "optional_weights",  # effect sizes / -log10 p
            "gin_edge_list_path",  # .sifi TSV
        }
    ),
    outputs=frozenset(
        {
            "subgraph_edges_tsv",
            "propagation_scores_tsv",
            "qc_mapping_json",
        }
    ),
    description=(
        "User supplies RNA/protein/metabolite identifiers (after mapping). "
        "Tool extracts a k-hop bipartite neighborhood in GIN, runs degree-aware "
        "diffusion from seeds, and writes an edge list plus ranked nodes for "
        "interpretation (chromatin–TF–target chains appear as paths)."
    ),
)

# US2: cohort matrices + outcomes -> graph-regularized or graph-guided features
US2 = UserStory(
    id="US2",
    title="Cohort multi-omics with GIN prior for benchmarking",
    inputs=frozenset(
        {
            "sample_by_feature_matrices",  # per modality, aligned samples
            "clinical_or_label",  # e.g. survival, response
            "gin_prior_paths",  # edge list or precomputed sparse blocks
        }
    ),
    outputs=frozenset(
        {
            "baseline_vs_gin_metrics_json",
            "optional_propagated_features_npy",
        }
    ),
    description=(
        "User runs a fixed benchmark protocol (see gin_tool.benchmark): "
        "e.g. univariate + SNF-style fusion vs seed-diffusion features on the "
        "same labels, reporting AUROC/C-index deltas. GIN is optional prior, "
        "not a replacement for standard fusion."
    ),
)

USER_STORIES: tuple[UserStory, ...] = (US1, US2)

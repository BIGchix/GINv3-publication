"""Benchmark protocol: public multi-omics cohorts vs GIN-augmented features (scaffold)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any


@dataclass
class BenchmarkProtocol:
    """
    Documented plan for comparing **baseline** vs **GIN-prior** models (US2).

    Recommended public resources (no download in this package):
    - TCGA multi-omics: RNA-seq + clinical (optional methylation) via GDC; use
      held-out AUROC for binary outcomes or Harrell's C-index for survival.
    - Cancer cell lines: CCLE/GDSC drug response + expression (single-omics
      baseline; extend with metabolomics subsets where available).
    """

    name: str
    cohort: str
    modalities: tuple[str, ...]
    label: str
    baseline_methods: tuple[str, ...] = ("elastic_net", "snf_fusion_stub")
    gin_methods: tuple[str, ...] = ("rwr_top_features",)
    splits: str = "5fold_cv"
    metrics: tuple[str, ...] = ("auroc", "c_index")
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "cohort": self.cohort,
            "modalities": list(self.modalities),
            "label": self.label,
            "baseline_methods": list(self.baseline_methods),
            "gin_methods": list(self.gin_methods),
            "splits": self.splits,
            "metrics": list(self.metrics),
            "notes": self.notes,
        }


DEFAULT_TCGA_STYLE = BenchmarkProtocol(
    name="tcga_style_binary",
    cohort="TCGA (example: BRCA stage I vs II, user-defined)",
    modalities=("rna", "clinical"),
    label="binary_stage_or_subtype",
    baseline_methods=("logistic_l1", "similarity_network_fusion"),
    gin_methods=("diffusion_smoothed_expression", "rwr_seeds_from_de"),
    notes=(
        "Preprocess RNA to gene symbols; map to GIN .rna nodes. Baseline: "
        "penalized logistic on top variance genes. GIN: replace or augment "
        "features with RWR-smoothed scores on k-hop subgraph union over DE seeds. "
        "Compare mean AUROC +/- SE across folds; report effect size (Cliff's delta)."
    ),
)


@dataclass
class BenchmarkResult:
    protocol_name: str
    metric: str
    baseline_value: float
    gin_value: float
    details: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> str:
        return json.dumps(
            {
                "protocol_name": self.protocol_name,
                "metric": self.metric,
                "baseline_value": self.baseline_value,
                "gin_value": self.gin_value,
                "delta": self.gin_value - self.baseline_value,
                "details": self.details,
            },
            indent=2,
        )


def compare_scalar_metrics(
    protocol: BenchmarkProtocol,
    metric: str,
    baseline_value: float,
    gin_value: float,
    **details: Any,
) -> BenchmarkResult:
    """Wrap numeric comparison after user runs sklearn/survival models externally."""
    return BenchmarkResult(
        protocol_name=protocol.name,
        metric=metric,
        baseline_value=baseline_value,
        gin_value=gin_value,
        details={"protocol": protocol.to_dict(), **details},
    )


def toy_demonstration_scores() -> BenchmarkResult:
    """
    Placeholder numbers showing JSON shape; replace with real CV outputs.

    Intended workflow: export CV metrics from a notebook, then call
    ``compare_scalar_metrics``.
    """
    return compare_scalar_metrics(
        DEFAULT_TCGA_STYLE,
        metric="auroc",
        baseline_value=0.72,
        gin_value=0.76,
        note="Synthetic illustration only",
    )

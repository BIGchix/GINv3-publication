"""QC mapping dashboard."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from gin_tool.mapping import MappingQCReport


def render_qc_dashboard(viz_dir: Path, report: MappingQCReport) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # Mapped vs unmapped
    ax = axes[0]
    mapped = report.mapped_count
    unmapped = report.input_count - report.mapped_count
    ax.bar(
        ["mapped", "unmapped"],
        [mapped, unmapped],
        color=["#27ae60", "#e74c3c"],
    )
    ax.set_title("Seed mapping")
    ax.set_ylabel("count")

    # Layer counts (mapped seeds only)
    ax = axes[1]
    layers = report.layer_counts
    if layers:
        ax.pie(
            list(layers.values()),
            labels=list(layers.keys()),
            autopct="%1.0f%%",
            startangle=90,
        )
    ax.set_title("Mapped seed layers")

    # Unmapped reasons
    ax = axes[2]
    reason_counts: Counter[str] = Counter()
    for row in report.unmapped_detail:
        reason_counts[row.get("reason", "unknown")] += 1
    if report.ambiguous_mappings:
        reason_counts["ambiguous_resolved"] = len(report.ambiguous_mappings)
    if reason_counts:
        ax.bar(list(reason_counts.keys()), list(reason_counts.values()), color="#3498db")
        ax.tick_params(axis="x", rotation=30)
    else:
        ax.text(0.5, 0.5, "no failures", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("Mapping issues")

    fig.suptitle(f"GIN mapping QC ({report.gin_version})", fontsize=12)
    fig.tight_layout()
    fig.savefig(viz_dir / "qc_dashboard.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

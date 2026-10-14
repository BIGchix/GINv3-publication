"""Locked pipeline version strings for context subgraph + SPN showcase."""

from __future__ import annotations

# v3 (fallback): Directed export + promoter reverse BFS + ancestor-only prune +
# dogma then HME/eraser + histprune.
CONTEXT_PIPELINE_VERSION_V3 = "context_v3_ancprune_dogma_hme_histprune"

# v4 (default): HME soft-seeds before extract; LCC → met-enzyme recovery
# (post-topology only) → dogma → HME (seed-only active-hist backfill + drop
# non-seed promoter.active/;state, keep RNA/protein + eraser/histprune).
# Dogma precedes HME so seed promoter.active backfill sees the full active set.
# Metabolic hyperedge enzymes are NOT soft-seeded at 0b (avoids BFS blow-up).
CONTEXT_PIPELINE_VERSION_V4 = (
    "context_v4_hme_soft_metenz3a_dogma_hme_seedhist_dropnsact_histprune"
)

CONTEXT_PIPELINE_VERSION = CONTEXT_PIPELINE_VERSION_V4

CONTEXT_PIPELINE_ALIASES: dict[str, str] = {
    "v3": CONTEXT_PIPELINE_VERSION_V3,
    "v4": CONTEXT_PIPELINE_VERSION_V4,
    CONTEXT_PIPELINE_VERSION_V3: CONTEXT_PIPELINE_VERSION_V3,
    CONTEXT_PIPELINE_VERSION_V4: CONTEXT_PIPELINE_VERSION_V4,
}

# Re-export SPN policy versions for a single import surface.
from gin_tool.subgraph_pipeline.spn_policy import (  # noqa: E402
    SPN_PIPELINE_VERSION,
    SPN_PIPELINE_VERSION_V2,
    SPN_PIPELINE_VERSION_V8,
)

__all__ = [
    "CONTEXT_PIPELINE_VERSION",
    "CONTEXT_PIPELINE_VERSION_V3",
    "CONTEXT_PIPELINE_VERSION_V4",
    "CONTEXT_PIPELINE_ALIASES",
    "SPN_PIPELINE_VERSION",
    "SPN_PIPELINE_VERSION_V2",
    "SPN_PIPELINE_VERSION_V8",
]

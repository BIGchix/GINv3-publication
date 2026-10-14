"""Resolve metabolic stoichiometry coefficients for SPN arcs (phase C).

When ``stoich=kegg``, look up BioPAX/KEGGR ν for substrate/product places.
Unresolved edges fall back to 1:1 (conserving transfer). Coverage is reported
for QC honesty — typical histprune subgraphs have ~11% KEGGR-tagged met edges.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping


@dataclass
class StoichCoverage:
    n_met_edges: int = 0
    n_keggr_tagged: int = 0
    n_resolved: int = 0
    n_fallback_one_to_one: int = 0
    unresolved_ids: list[str] = field(default_factory=list)
    cache_path: str | None = None

    def to_dict(self) -> dict[str, Any]:
        frac = (self.n_resolved / self.n_met_edges) if self.n_met_edges else 0.0
        return {
            "n_met_edges": self.n_met_edges,
            "n_keggr_tagged": self.n_keggr_tagged,
            "n_resolved": self.n_resolved,
            "n_fallback_one_to_one": self.n_fallback_one_to_one,
            "resolved_fraction": frac,
            "unresolved_ids_sample": self.unresolved_ids[:20],
            "cache_path": self.cache_path,
            "status": "stub" if self.n_resolved == 0 else "ok",
        }


def load_biopax_stoich_cache(path: str | Path | None) -> dict[str, dict[str, float]]:
    """Load optional JSON cache: ``{KEGGRxxxxx: {place: nu, ...}}``.

    Returns empty dict when missing — callers fall back to 1:1.
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.is_file():
        return {}
    import json

    data = json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return {}
    out: dict[str, dict[str, float]] = {}
    for rid, coeffs in data.items():
        if isinstance(coeffs, dict):
            out[str(rid)] = {str(k): float(v) for k, v in coeffs.items()}
    return out


def resolve_reaction_coeffs(
    transition_id: str,
    *,
    cache: Mapping[str, Mapping[str, float]] | None = None,
) -> dict[str, float] | None:
    """Return place→ν for a KEGGR transition, or None → use 1:1."""
    if not cache:
        return None
    tid = transition_id
    # Allow ids embedded in longer virtual names.
    for key, coeffs in cache.items():
        if key in tid:
            return dict(coeffs)
    return None


def report_stoich_coverage(
    transition_ids: list[str],
    *,
    met_transition_ids: list[str] | None = None,
    cache: Mapping[str, Mapping[str, float]] | None = None,
    cache_path: str | None = None,
) -> StoichCoverage:
    """Coverage QC for metabolic transitions (KEGGR resolve vs 1:1 fallback)."""
    mets = list(met_transition_ids) if met_transition_ids is not None else list(transition_ids)
    cov = StoichCoverage(n_met_edges=len(mets), cache_path=cache_path)
    cache = cache or {}
    for tid in mets:
        tagged = "KEGGR" in tid.upper() or "keggr" in tid.lower()
        if tagged:
            cov.n_keggr_tagged += 1
        coeffs = resolve_reaction_coeffs(tid, cache=cache)
        if coeffs:
            cov.n_resolved += 1
        else:
            cov.n_fallback_one_to_one += 1
            if tagged:
                cov.unresolved_ids.append(tid)
    return cov

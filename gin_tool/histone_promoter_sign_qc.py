"""Audit signed histone → promoter edges (activation vs inhibition)."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


def _is_histone(name: str) -> bool:
    u = name.upper()
    return any(h in u for h in ("H3K", "H4K", "H2AK", "H2BK", "H3AC", "H4AC"))


def _is_promoter_state(name: str) -> bool:
    return ".promoter;state" in name.lower()


def _is_promoter_active(name: str) -> bool:
    return name.lower().endswith(".promoter.active")


def _mark_class(name: str) -> str:
    u = name.upper()
    if "K27ME" in u or "K9ME" in u:
        return "repressive_like"
    if "K4ME" in u or "K36ME" in u or "K79ME" in u:
        return "active_like"
    if "AC" in u:
        return "acetyl_active_like"
    return "other"


@dataclass
class HistonePromoterSignQC:
    n_hist_to_promoter_state: int = 0
    n_inhibitory: int = 0
    n_activating: int = 0
    n_zero_or_other: int = 0
    by_mark_class: dict[str, dict[str, int]] = field(default_factory=dict)
    n_paths_to_promoter_active: int = 0
    n_inhibitory_paths_to_active: int = 0
    sample_inhibitory: list[dict] = field(default_factory=list)
    sample_activating: list[dict] = field(default_factory=list)
    note: str = (
        "GIN stores histone → GENE.promoter;state (signed), then "
        "promoter;state → promoter.active (+1). SPN maps w_in*w_out<0 to inhibition."
    )

    def to_dict(self) -> dict:
        return {
            "n_hist_to_promoter_state": self.n_hist_to_promoter_state,
            "n_inhibitory": self.n_inhibitory,
            "n_activating": self.n_activating,
            "n_zero_or_other": self.n_zero_or_other,
            "by_mark_class": self.by_mark_class,
            "n_paths_to_promoter_active": self.n_paths_to_promoter_active,
            "n_inhibitory_paths_to_active": self.n_inhibitory_paths_to_active,
            "sample_inhibitory": self.sample_inhibitory,
            "sample_activating": self.sample_activating,
            "note": self.note,
        }


def audit_histone_promoter_signs(
    edges: Iterable[tuple[str, str, float]],
) -> HistonePromoterSignQC:
    edge_list = list(edges)
    out: dict[str, list[tuple[str, float]]] = {}
    for a, b, w in edge_list:
        out.setdefault(a, []).append((b, float(w)))

    qc = HistonePromoterSignQC()
    class_counts: dict[str, Counter[str]] = {}

    for a, b, w in edge_list:
        if not (_is_histone(a) and _is_promoter_state(b)):
            continue
        qc.n_hist_to_promoter_state += 1
        wf = float(w)
        cls = _mark_class(a)
        bucket = class_counts.setdefault(cls, Counter())
        if wf < 0:
            qc.n_inhibitory += 1
            bucket["inhibitory"] += 1
            if len(qc.sample_inhibitory) < 12:
                qc.sample_inhibitory.append({"source": a, "target": b, "weight": wf})
        elif wf > 0:
            qc.n_activating += 1
            bucket["activating"] += 1
            if len(qc.sample_activating) < 12:
                qc.sample_activating.append({"source": a, "target": b, "weight": wf})
        else:
            qc.n_zero_or_other += 1
            bucket["other"] += 1

        for c, w2 in out.get(b, []):
            if not _is_promoter_active(c):
                continue
            qc.n_paths_to_promoter_active += 1
            if wf * float(w2) < 0:
                qc.n_inhibitory_paths_to_active += 1

    qc.by_mark_class = {k: dict(v) for k, v in sorted(class_counts.items())}
    return qc


def audit_histone_promoter_signs_from_tsv(path: str | Path) -> HistonePromoterSignQC:
    rows: list[tuple[str, str, float]] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        w = float(parts[2]) if len(parts) > 2 else 1.0
        rows.append((parts[0], parts[1], w))
    return audit_histone_promoter_signs(rows)


def write_histone_promoter_sign_qc(path: str | Path, qc: HistonePromoterSignQC) -> Path:
    out = Path(path)
    out.write_text(json.dumps(qc.to_dict(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out

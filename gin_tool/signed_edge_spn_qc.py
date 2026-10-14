"""QC: signed GIN edges → SPN activation/inhibition kinds (all layers, not only histone)."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from gin_tool.mapping import is_virtual_gin_node


def _is_histone(name: str) -> bool:
    u = name.upper()
    return any(h in u for h in ("H3K", "H4K", "H2AK", "H2BK", "H3AC", "H4AC"))


@dataclass
class SignedEdgeSpnQC:
    n_neg_gin_edges: int = 0
    n_neg_arcs_in_petri: int = 0
    n_neg_dropped: int = 0
    n_spn_inhibition: int = 0
    n_spn_activation: int = 0
    n_expected_rvr_inhibition: int = 0
    n_got_rvr_inhibition: int = 0
    n_both_negative_rvr: int = 0
    inhibition_by_source_kind: dict[str, int] = field(default_factory=dict)
    sample_both_neg: list[dict] = field(default_factory=list)
    sample_missed: list[dict] = field(default_factory=list)
    note: str = (
        "SPN inhibition if w_in<0 or w_out<0 on R→V→R (covers histone XOR signs "
        "and protein both-negative inhibition encoding)."
    )

    def to_dict(self) -> dict:
        return {
            "n_neg_gin_edges": self.n_neg_gin_edges,
            "n_neg_arcs_in_petri": self.n_neg_arcs_in_petri,
            "n_neg_dropped": self.n_neg_dropped,
            "n_spn_inhibition": self.n_spn_inhibition,
            "n_spn_activation": self.n_spn_activation,
            "n_expected_rvr_inhibition": self.n_expected_rvr_inhibition,
            "n_got_rvr_inhibition": self.n_got_rvr_inhibition,
            "n_both_negative_rvr": self.n_both_negative_rvr,
            "inhibition_by_source_kind": self.inhibition_by_source_kind,
            "sample_both_neg": self.sample_both_neg,
            "sample_missed": self.sample_missed,
            "ok": self.n_expected_rvr_inhibition == self.n_got_rvr_inhibition
            and self.n_neg_dropped == 0,
            "note": self.note,
        }


def _src_kind(name: str) -> str:
    if _is_histone(name):
        return "histone"
    low = name.lower()
    if low.endswith(".protein"):
        return "protein"
    if ".promoter" in low:
        return "promoter"
    if low.endswith(".rna"):
        return "rna"
    if name.startswith(("CHEBI", "SID", "CID", "HMDB")):
        return "metabolite"
    return "other"


def audit_signed_edges_in_spn(
    gin_edges: Iterable[tuple[str, str, float]],
    *,
    arc_weights: dict[tuple[str, str], float],
    arcs: Iterable[tuple[str, str]],
    interactions: Iterable[object],
) -> SignedEdgeSpnQC:
    edges = list(gin_edges)
    arc_set = set(arcs)
    qc = SignedEdgeSpnQC()
    neg = [(a, b, float(w)) for a, b, w in edges if float(w) < 0]
    qc.n_neg_gin_edges = len(neg)
    qc.n_neg_arcs_in_petri = sum(1 for w in arc_weights.values() if w < 0)
    qc.n_neg_dropped = sum(1 for a, b, _ in neg if (a, b) not in arc_set)

    inh_set: set[tuple[str, str]] = set()
    src_counts: Counter[str] = Counter()
    n_inh = n_act = 0
    for i in interactions:
        kind = getattr(i, "kind", None)
        if kind == "inhibition":
            n_inh += 1
            inh_set.add((getattr(i, "source"), getattr(i, "target")))
            src_counts[_src_kind(getattr(i, "source"))] += 1
        elif kind == "activation":
            n_act += 1
    qc.n_spn_inhibition = n_inh
    qc.n_spn_activation = n_act
    qc.inhibition_by_source_kind = dict(src_counts)

    outs: dict[str, list[tuple[str, float]]] = {}
    for a, b, w in edges:
        outs.setdefault(a, []).append((b, float(w)))

    for a, b, w in edges:
        if is_virtual_gin_node(a) or not is_virtual_gin_node(b):
            continue
        if (a, b) not in arc_set:
            continue
        for c, w2 in outs.get(b, []):
            if is_virtual_gin_node(c):
                continue
            if not (w < 0 or w2 < 0):
                continue
            qc.n_expected_rvr_inhibition += 1
            if w < 0 and w2 < 0:
                qc.n_both_negative_rvr += 1
                if len(qc.sample_both_neg) < 12:
                    qc.sample_both_neg.append(
                        {"source": a, "via": b, "target": c, "w_in": w, "w_out": w2}
                    )
            if (a, c) in inh_set:
                qc.n_got_rvr_inhibition += 1
            elif len(qc.sample_missed) < 12:
                qc.sample_missed.append(
                    {"source": a, "via": b, "target": c, "w_in": w, "w_out": w2}
                )
    return qc


def write_signed_edge_spn_qc(path: str | Path, qc: SignedEdgeSpnQC) -> Path:
    out = Path(path)
    out.write_text(json.dumps(qc.to_dict(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out

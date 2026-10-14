"""ChromHMM emission → signed mark weights for hist→promoter.active (SPN v8).

Design (validation / non-locus):
  - Load a cell-type emission matrix ``e[state, mark] = P(mark=1 | state)``.
  - Keep only rows annotated as **promoter** compartment.
  - Split those into active vs inactive transcription sets.
  - Signed weight per mark (global, not locus-specific)::

        w_j = mean_{i in A} e_{ij}  -  mean_{i in I} e_{ij}

  - At REGULATORY fire (hist → ``*.promoter.active``): ``δ_signed = δ · w_j``;
    sign(δ_signed) selects activation vs inhibition; magnitude is ``|δ_signed|``.

Emission tables are **not** P(state|mark); the difference above is the agreed
discriminative proxy for this pipeline.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

# Canonical mark token inside GIN place names (e.g. TP53.promoter.H3K27ac).
_MARK_TOKEN_RE = re.compile(
    r"(H3K\d+[A-Za-z0-9]*|H4K\d+[A-Za-z0-9]*|H2A[A-Za-z0-9]+|H2B[A-Za-z0-9]+|"
    r"H3AC[A-Za-z0-9]*|H4AC[A-Za-z0-9]*)",
    re.IGNORECASE,
)

# Compartment values kept for promoter-only filtering.
DEFAULT_PROMOTER_COMPARTMENTS = frozenset(
    {
        "promoter",
        "active_promoter",
        "flanking_promoter",
        "poised_promoter",
        "bivalent_promoter",
        "tss",
        "promoter_flanking",
    }
)

DEFAULT_ACTIVE_LABELS = frozenset(
    {
        "active",
        "active_transcription",
        "tx_active",
        "strong",
    }
)

DEFAULT_INACTIVE_LABELS = frozenset(
    {
        "inactive",
        "inactive_transcription",
        "tx_inactive",
        "repressed",
        "quiescent",
        "poised",  # promoter-poised counted inactive for Tx unless overridden
        "bivalent",
    }
)


@dataclass(frozen=True)
class ChromHmmEmissionTable:
    """Rows = chromatin states; columns = histone marks; values in [0, 1]."""

    marks: tuple[str, ...]
    state_ids: tuple[str, ...]
    mnemonics: tuple[str, ...]
    compartments: tuple[str, ...]
    activities: tuple[str, ...]
    emission: tuple[tuple[float, ...], ...]  # n_states × n_marks

    def to_dict(self) -> dict:
        return {
            "marks": list(self.marks),
            "n_states": len(self.state_ids),
            "state_ids": list(self.state_ids),
            "mnemonics": list(self.mnemonics),
            "compartments": list(self.compartments),
            "activities": list(self.activities),
        }


@dataclass
class ChromHmmWeightResult:
    """Signed mark weights plus QC for one cell-type emission table."""

    weights: dict[str, float]
    n_promoter_states: int = 0
    n_active_states: int = 0
    n_inactive_states: int = 0
    active_state_ids: list[str] = field(default_factory=list)
    inactive_state_ids: list[str] = field(default_factory=list)
    skipped_non_promoter: int = 0
    formula: str = "w_j = mean_A e_ij - mean_I e_ij"
    source_path: str | None = None

    def to_dict(self) -> dict:
        return {
            "formula": self.formula,
            "source_path": self.source_path,
            "n_marks": len(self.weights),
            "weights": dict(sorted(self.weights.items())),
            "n_promoter_states": self.n_promoter_states,
            "n_active_states": self.n_active_states,
            "n_inactive_states": self.n_inactive_states,
            "active_state_ids": list(self.active_state_ids),
            "inactive_state_ids": list(self.inactive_state_ids),
            "skipped_non_promoter": self.skipped_non_promoter,
        }


def normalize_mark_name(name: str) -> str:
    """Uppercase chromHMM / GIN mark token (H3K27ac → H3K27AC)."""
    return str(name).strip().upper().replace("-", "").replace("_", "")


def histone_mark_key_from_place(place: str) -> str | None:
    """Extract chromHMM mark key from a GIN hist place name, or None."""
    if not place or ";" in place:
        return None
    m = _MARK_TOKEN_RE.search(place)
    if not m:
        return None
    return normalize_mark_name(m.group(1))


def _norm_label(s: str) -> str:
    return str(s).strip().lower().replace(" ", "_").replace("-", "_")


def _is_promoter_compartment(
    compartment: str, *, allow: frozenset[str] = DEFAULT_PROMOTER_COMPARTMENTS
) -> bool:
    c = _norm_label(compartment)
    if c in allow:
        return True
    return "promoter" in c or c == "tss"


def _is_active_activity(activity: str, *, allow: frozenset[str] = DEFAULT_ACTIVE_LABELS) -> bool:
    a = _norm_label(activity)
    if a in allow:
        return True
    return a.startswith("active") and "inactive" not in a


def _is_inactive_activity(
    activity: str, *, allow: frozenset[str] = DEFAULT_INACTIVE_LABELS
) -> bool:
    a = _norm_label(activity)
    if a in allow:
        return True
    return (
        a.startswith("inactive")
        or "repress" in a
        or a in {"quiescent", "quiet", "hetero", "heterochromatin"}
    )


def load_chromhmm_emission_tsv(path: str | Path) -> ChromHmmEmissionTable:
    """Load emission TSV.

    Required header columns:
      ``state_id``, ``mnemonic``, ``compartment``, ``activity``, then mark columns.

    ``compartment``: promoter / genebody / enhancer / heterochromatin / …
    ``activity``: active / inactive / other (only promoter rows used for weights).
    Mark columns: numeric P(mark|state) in [0, 1].
    """
    p = Path(path)
    lines = [
        ln
        for ln in p.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]
    if not lines:
        raise ValueError(f"empty chromHMM emission table: {p}")
    header = lines[0].split("\t")
    required = ("state_id", "mnemonic", "compartment", "activity")
    lower = {h.lower(): i for i, h in enumerate(header)}
    for req in required:
        if req not in lower:
            raise ValueError(f"emission TSV missing column {req!r}: {p}")
    meta = {required[i] for i in range(len(required))}
    mark_cols = [(i, h) for i, h in enumerate(header) if h.lower() not in meta]
    if not mark_cols:
        raise ValueError(f"emission TSV has no mark columns: {p}")

    marks = tuple(normalize_mark_name(h) for _, h in mark_cols)
    state_ids: list[str] = []
    mnemonics: list[str] = []
    compartments: list[str] = []
    activities: list[str] = []
    rows: list[tuple[float, ...]] = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) < len(header):
            parts = parts + [""] * (len(header) - len(parts))
        state_ids.append(parts[lower["state_id"]].strip())
        mnemonics.append(parts[lower["mnemonic"]].strip())
        compartments.append(parts[lower["compartment"]].strip())
        activities.append(parts[lower["activity"]].strip())
        vals: list[float] = []
        for i, _h in mark_cols:
            raw = parts[i].strip() if i < len(parts) else ""
            vals.append(float(raw) if raw else 0.0)
        rows.append(tuple(vals))
    return ChromHmmEmissionTable(
        marks=marks,
        state_ids=tuple(state_ids),
        mnemonics=tuple(mnemonics),
        compartments=tuple(compartments),
        activities=tuple(activities),
        emission=tuple(rows),
    )


def signed_mark_weights_from_emission(
    table: ChromHmmEmissionTable,
    *,
    promoter_compartments: frozenset[str] | None = None,
    active_labels: frozenset[str] | None = None,
    inactive_labels: frozenset[str] | None = None,
    source_path: str | None = None,
) -> ChromHmmWeightResult:
    """Compute ``w_j = mean_A e_ij - mean_I e_ij`` on promoter states only."""
    prom_allow = promoter_compartments or DEFAULT_PROMOTER_COMPARTMENTS
    act_allow = active_labels or DEFAULT_ACTIVE_LABELS
    inact_allow = inactive_labels or DEFAULT_INACTIVE_LABELS

    active_idx: list[int] = []
    inactive_idx: list[int] = []
    skipped = 0
    for i, (comp, act) in enumerate(zip(table.compartments, table.activities)):
        if not _is_promoter_compartment(comp, allow=prom_allow):
            skipped += 1
            continue
        if _is_active_activity(act, allow=act_allow):
            active_idx.append(i)
        elif _is_inactive_activity(act, allow=inact_allow):
            inactive_idx.append(i)
        else:
            # Promoter but unlabeled activity — skip for A/I means.
            skipped += 1

    if not active_idx or not inactive_idx:
        raise ValueError(
            "need ≥1 active and ≥1 inactive *promoter* states to form signed weights "
            f"(active={len(active_idx)}, inactive={len(inactive_idx)}, "
            f"skipped_or_nonpromoter={skipped})"
        )

    n_marks = len(table.marks)
    weights: dict[str, float] = {}
    for j, mark in enumerate(table.marks):
        mean_a = sum(table.emission[i][j] for i in active_idx) / len(active_idx)
        mean_i = sum(table.emission[i][j] for i in inactive_idx) / len(inactive_idx)
        weights[mark] = float(mean_a - mean_i)

    return ChromHmmWeightResult(
        weights=weights,
        n_promoter_states=len(active_idx) + len(inactive_idx),
        n_active_states=len(active_idx),
        n_inactive_states=len(inactive_idx),
        active_state_ids=[table.state_ids[i] for i in active_idx],
        inactive_state_ids=[table.state_ids[i] for i in inactive_idx],
        skipped_non_promoter=skipped,
        source_path=source_path,
    )


def load_signed_mark_weights_tsv(path: str | Path) -> dict[str, float]:
    """Load precomputed ``mark\\tweight`` (header optional)."""
    p = Path(path)
    out: dict[str, float] = {}
    for i, ln in enumerate(p.read_text(encoding="utf-8").splitlines()):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < 2:
            continue
        if i == 0 and parts[0].lower() in {"mark", "mark_name", "histone"}:
            continue
        out[normalize_mark_name(parts[0])] = float(parts[1])
    if not out:
        raise ValueError(f"no mark weights in {p}")
    return out


def resolve_chromhmm_mark_weights(
    *,
    emission_tsv: str | Path | None = None,
    weights_tsv: str | Path | None = None,
) -> ChromHmmWeightResult | None:
    """Prefer precomputed weights TSV; else derive from emission TSV."""
    if weights_tsv:
        w = load_signed_mark_weights_tsv(weights_tsv)
        return ChromHmmWeightResult(
            weights=w,
            formula="precomputed mark\\tweight TSV",
            source_path=str(weights_tsv),
        )
    if emission_tsv:
        table = load_chromhmm_emission_tsv(emission_tsv)
        return signed_mark_weights_from_emission(
            table, source_path=str(emission_tsv)
        )
    return None


def lookup_mark_weight(
    weights: Mapping[str, float] | None, place: str
) -> float | None:
    """Return signed w for a hist place, or None if table missing / mark unknown."""
    if not weights:
        return None
    key = histone_mark_key_from_place(place)
    if key is None:
        return None
    if key in weights:
        return float(weights[key])
    # soft match: strip trailing digits variants already normalized
    for k, v in weights.items():
        if key.startswith(k) or k.startswith(key):
            return float(v)
    return None


def apply_signed_chromhmm_drive(
    delta: float,
    *,
    source: str,
    target: str,
    kind: str,
    weights: Mapping[str, float] | None,
) -> tuple[float, str, bool]:
    """Scale REGULATORY δ by signed w for hist→promoter.active.

    Returns ``(delta_abs, kind, applied)``. When no weight applies, returns
    the original ``(delta, kind, False)``.
    """
    from gin_tool.spn_export import _is_histone_mark

    if not weights or delta <= 0:
        return float(delta), kind, False
    if not _is_histone_mark(source):
        return float(delta), kind, False
    tgt = target.lower()
    if ".promoter.active" not in tgt and not tgt.endswith("promoter.active"):
        return float(delta), kind, False
    w = lookup_mark_weight(weights, source)
    if w is None:
        return float(delta), kind, False
    signed = float(delta) * float(w)
    if abs(signed) <= 1e-15:
        return 0.0, kind, True
    new_kind = "activation" if signed > 0.0 else "inhibition"
    return abs(signed), new_kind, True


def write_mark_weights_json(result: ChromHmmWeightResult, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result.to_dict(), indent=2) + "\n", encoding="utf-8")
    return out


def write_mark_weights_tsv(weights: Mapping[str, float], path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["mark\tweight"]
    for k in sorted(weights):
        lines.append(f"{k}\t{float(weights[k]):.6g}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out

"""Assemble per-cell chromHMM emission tables for SPN v8.

Policy (panel):
  - Core marks (H3K4me3/me1, H3K36me3, H3K27me3, H3K9me3): **indiv** 15-state
    emissions per epigenome.
  - **H3K27ac**: from **joint** 18-state emission, looked up by mnemonic
    (not by state index).
  - Align indiv states → joint-15 mnemonics by nearest emission vector
    (same 5 marks), never ``state_i == state_i``.
  - Keep only **promoter** mnemonics; split **active** vs **inactive**.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from gin_tool.subgraph_pipeline.chromhmm_emission_weights import (
    ChromHmmEmissionTable,
    normalize_mark_name,
    signed_mark_weights_from_emission,
    write_mark_weights_json,
    write_mark_weights_tsv,
)

# Joint-15 mnemonic → (compartment, activity) for SPN v8 filtering.
_MNEMONIC_ANNOT: dict[str, tuple[str, str]] = {
    "TSSA": ("promoter", "active"),
    "TSSAFLNK": ("promoter", "active"),
    "TSSFLNK": ("promoter", "active"),
    "TSSFLNKU": ("promoter", "active"),
    "TSSFLNKD": ("promoter", "active"),
    "TSSBIV": ("promoter", "inactive"),
    "BIVFLNK": ("promoter", "inactive"),
    # non-promoter (kept for matching QC; filtered later)
    "TXFLNK": ("genebody", "active"),
    "TX": ("genebody", "active"),
    "TXWK": ("genebody", "active"),
    "ENHG": ("enhancer", "active"),
    "ENHG1": ("enhancer", "active"),
    "ENHG2": ("enhancer", "active"),
    "ENH": ("enhancer", "active"),
    "ENHA1": ("enhancer", "active"),
    "ENHA2": ("enhancer", "active"),
    "ENHWK": ("enhancer", "active"),
    "ENHBIV": ("enhancer", "inactive"),
    "ZNFRPTS": ("other", "inactive"),
    "ZNF/RPTS": ("other", "inactive"),
    "HET": ("heterochromatin", "inactive"),
    "REPRPC": ("other", "inactive"),
    "REPRPCWK": ("other", "inactive"),
    "QUIES": ("other", "inactive"),
}

# Joint-18 mnemonic → which joint-18 rows supply H3K27ac for a joint-15 label.
_K27AC_FROM_JOINT18: dict[str, tuple[str, ...]] = {
    "TSSA": ("TSSA",),
    "TSSAFLNK": ("TSSFLNK", "TSSFLNKU", "TSSFLNKD"),
    "TSSBIV": ("TSSBIV",),
    "BIVFLNK": ("TSSBIV", "ENHBIV"),
}


@dataclass(frozen=True)
class RawEmissionMatrix:
    """Roadmap-style emission file (state id may be ``1`` or ``1_TssA``)."""

    marks: tuple[str, ...]
    state_ids: tuple[str, ...]
    mnemonics: tuple[str, ...]  # may be "" if numeric-only
    emission: np.ndarray  # (n_states, n_marks)


def _parse_state_label(raw: str) -> tuple[str, str]:
    s = raw.strip()
    if "_" in s:
        sid, mnem = s.split("_", 1)
        return sid, mnem
    return s, ""


def _norm_mnemonic(m: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", m.upper())


def load_roadmap_emission_matrix(path: str | Path) -> RawEmissionMatrix:
    p = Path(path)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"empty emission file: {p}")
    header = lines[0].split("\t")
    marks = tuple(normalize_mark_name(h) for h in header[1:])
    state_ids: list[str] = []
    mnemonics: list[str] = []
    rows: list[list[float]] = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        sid, mnem = _parse_state_label(parts[0])
        state_ids.append(sid)
        mnemonics.append(mnem)
        rows.append([float(x) for x in parts[1 : 1 + len(marks)]])
    return RawEmissionMatrix(
        marks=marks,
        state_ids=tuple(state_ids),
        mnemonics=tuple(mnemonics),
        emission=np.asarray(rows, dtype=float),
    )


def load_overlap_tss_enrichment(path: str | Path) -> dict[str, float]:
    """Map state id → TSS fold enrichment from ChromHMM overlap table."""
    p = Path(path)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
    header = lines[0].split("\t")
    # Prefer precise TSS column.
    tss_col = None
    for name in header:
        if name.startswith("TSS.") and "2kbp" not in name:
            tss_col = name
            break
    if tss_col is None:
        for name in header:
            if "TSS" in name:
                tss_col = name
                break
    if tss_col is None:
        raise ValueError(f"no TSS column in overlap table: {p}")
    j = header.index(tss_col)
    out: dict[str, float] = {}
    for ln in lines[1:]:
        if ln.lower().startswith("base"):
            continue
        parts = ln.split("\t")
        sid, _ = _parse_state_label(parts[0])
        out[sid] = float(parts[j])
    return out


def mnemonic_annotation(mnemonic: str) -> tuple[str, str] | None:
    key = _norm_mnemonic(mnemonic)
    if key in _MNEMONIC_ANNOT:
        return _MNEMONIC_ANNOT[key]
    # soft: starts with TSS / BIVFLNK
    if key.startswith("TSSA") or key in {"TSSFLNK", "TSSFLNKU", "TSSFLNKD"}:
        return ("promoter", "active")
    if "TSSBIV" in key or key == "BIVFLNK":
        return ("promoter", "inactive")
    return None


def match_indiv_to_joint15(
    indiv: RawEmissionMatrix,
    joint15: RawEmissionMatrix,
) -> list[tuple[str, str, float]]:
    """For each indiv state, return (joint_state_id, joint_mnemonic, distance)."""
    shared = [m for m in joint15.marks if m in set(indiv.marks)]
    if len(shared) < 3:
        raise ValueError(f"too few shared marks for matching: {shared}")
    j_idx = [joint15.marks.index(m) for m in shared]
    i_idx = [indiv.marks.index(m) for m in shared]
    J = joint15.emission[:, j_idx]
    I = indiv.emission[:, i_idx]
    out: list[tuple[str, str, float]] = []
    for i in range(I.shape[0]):
        d = np.linalg.norm(J - I[i], axis=1)
        j = int(np.argmin(d))
        mnem = joint15.mnemonics[j] or joint15.state_ids[j]
        out.append((joint15.state_ids[j], mnem, float(d[j])))
    return out


def k27ac_for_mnemonic(joint18: RawEmissionMatrix, mnemonic: str) -> float:
    """Mean H3K27ac emission over joint-18 states mapped to this mnemonic family."""
    if "H3K27AC" not in joint18.marks:
        raise ValueError("joint-18 emission missing H3K27ac")
    col = joint18.marks.index("H3K27AC")
    key = _norm_mnemonic(mnemonic)
    targets = _K27AC_FROM_JOINT18.get(key)
    if not targets:
        # try direct match on joint18 mnemonic
        targets = (key,)
    vals: list[float] = []
    for i, m in enumerate(joint18.mnemonics):
        nm = _norm_mnemonic(m)
        if nm in {_norm_mnemonic(t) for t in targets} or nm == key:
            vals.append(float(joint18.emission[i, col]))
    if not vals:
        return 0.0
    return float(sum(vals) / len(vals))


@dataclass
class AssembledCellEmission:
    cell_id: str
    eid: str
    table: ChromHmmEmissionTable
    match_qc: list[dict]


def assemble_cell_emission(
    *,
    cell_id: str,
    eid: str,
    indiv_emission_path: str | Path,
    joint15_emission_path: str | Path,
    joint18_emission_path: str | Path,
    overlap_path: str | Path | None = None,
    min_tss_enrichment: float | None = 5.0,
    core_marks: Sequence[str] | None = None,
) -> AssembledCellEmission:
    """Build promoter-only annotated emission table for one panel cell.

    Columns: core marks from indiv + H3K27ac from joint-18 (by mnemonic).
    Rows: only promoter active/inactive after annotation alignment.
    """
    indiv = load_roadmap_emission_matrix(indiv_emission_path)
    joint15 = load_roadmap_emission_matrix(joint15_emission_path)
    joint18 = load_roadmap_emission_matrix(joint18_emission_path)
    matches = match_indiv_to_joint15(indiv, joint15)
    tss = load_overlap_tss_enrichment(overlap_path) if overlap_path else {}

    want_core = [
        normalize_mark_name(m)
        for m in (
            core_marks
            or ("H3K4me3", "H3K4me1", "H3K36me3", "H3K27me3", "H3K9me3")
        )
    ]
    for m in want_core:
        if m not in indiv.marks:
            raise ValueError(f"{cell_id}: indiv emission missing {m}")

    marks_out = tuple([*want_core, "H3K27AC"])
    state_ids: list[str] = []
    mnemonics: list[str] = []
    compartments: list[str] = []
    activities: list[str] = []
    rows: list[tuple[float, ...]] = []
    match_qc: list[dict] = []

    for i, sid in enumerate(indiv.state_ids):
        j_sid, j_mnem, dist = matches[i]
        ann = mnemonic_annotation(j_mnem)
        tss_v = float(tss.get(sid, float("nan")))
        qc = {
            "indiv_state": sid,
            "matched_joint15": f"{j_sid}_{j_mnem}" if j_mnem else j_sid,
            "distance": dist,
            "tss_enrichment": tss_v,
            "kept": False,
            "reason": "",
        }
        if ann is None:
            qc["reason"] = "no_mnemonic_annotation"
            match_qc.append(qc)
            continue
        comp, act = ann
        # Polycomb / bivalent states with strong TSS overlap → inactive promoters
        # when the cell model did not recover TssBiv/BivFlnk.
        mnem_norm = _norm_mnemonic(j_mnem)
        if (
            comp != "promoter"
            and mnem_norm in {"REPRPC", "REPRPCWK", "ENHBIV"}
            and min_tss_enrichment is not None
            and tss
            and np.isfinite(tss_v)
            and tss_v >= float(min_tss_enrichment)
        ):
            comp, act = "promoter", "inactive"
            qc["reason_override"] = f"{mnem_norm.lower()}_high_tss_as_inactive_promoter"
        if comp != "promoter":
            qc["reason"] = f"non_promoter:{comp}"
            match_qc.append(qc)
            continue
        if min_tss_enrichment is not None and tss and (
            not np.isfinite(tss_v) or tss_v < float(min_tss_enrichment)
        ):
            qc["reason"] = f"tss_below_{min_tss_enrichment}"
            match_qc.append(qc)
            continue
        core_vals = [float(indiv.emission[i, indiv.marks.index(m)]) for m in want_core]
        if mnem_norm in {"REPRPC", "REPRPCWK", "ENHBIV"} and act == "inactive":
            k27 = k27ac_for_mnemonic(joint18, "TssBiv")
        else:
            k27 = k27ac_for_mnemonic(joint18, j_mnem)
        state_ids.append(f"{sid}->{j_mnem}")
        mnemonics.append(j_mnem)
        compartments.append(comp)
        activities.append(act)
        rows.append(tuple([*core_vals, float(k27)]))
        qc["kept"] = True
        qc["reason"] = qc.get("reason_override") or f"promoter_{act}"
        qc["H3K27ac_joint18"] = k27
        match_qc.append(qc)

    if not rows:
        raise ValueError(
            f"{cell_id}/{eid}: no promoter states kept after annotation filter"
        )

    table = ChromHmmEmissionTable(
        marks=marks_out,
        state_ids=tuple(state_ids),
        mnemonics=tuple(mnemonics),
        compartments=tuple(compartments),
        activities=tuple(activities),
        emission=tuple(rows),
    )
    return AssembledCellEmission(
        cell_id=cell_id, eid=eid, table=table, match_qc=match_qc
    )


def write_spn_emission_tsv(table: ChromHmmEmissionTable, path: str | Path) -> Path:
    """Write SPN v8 loader format (state_id/mnemonic/compartment/activity + marks)."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    header = ["state_id", "mnemonic", "compartment", "activity", *table.marks]
    lines = ["\t".join(header)]
    for i, sid in enumerate(table.state_ids):
        row = [
            sid,
            table.mnemonics[i],
            table.compartments[i],
            table.activities[i],
            *[f"{v:.8g}" for v in table.emission[i]],
        ]
        lines.append("\t".join(row))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


DEFAULT_PANEL: dict[str, str] = {
    "h1_wa01": "E003",
    "h9_wa09": "E008",
    "imr90": "E017",
    "keratinocyte": "E127",
    "tcell_cd4": "E043",
}


def build_panel_emissions(
    roadmap_dir: str | Path,
    *,
    out_dir: str | Path,
    panel: Mapping[str, str] | None = None,
    min_tss_enrichment: float | None = 5.0,
) -> dict[str, dict]:
    """Assemble all panel cells; write emission TSVs + signed weight files."""
    root = Path(roadmap_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    joint15 = root / "emissions_15_coreMarks.txt"
    joint18 = root / "emissions_18_core_K27ac.txt"
    panel = dict(panel or DEFAULT_PANEL)
    summary: dict[str, dict] = {}
    for cell, eid in panel.items():
        indiv = root / "indiv_emissions_15" / f"{cell}_{eid}_emissions_15_coreMarks.txt"
        overlap = root / "indiv_overlap" / eid / f"{eid}_15_coreMarks_overlap.txt"
        assembled = assemble_cell_emission(
            cell_id=cell,
            eid=eid,
            indiv_emission_path=indiv,
            joint15_emission_path=joint15,
            joint18_emission_path=joint18,
            overlap_path=overlap if overlap.is_file() else None,
            min_tss_enrichment=min_tss_enrichment,
        )
        em_path = out / f"{cell}_{eid}_promoter_emission.tsv"
        write_spn_emission_tsv(assembled.table, em_path)
        weights = signed_mark_weights_from_emission(
            assembled.table, source_path=str(em_path)
        )
        write_mark_weights_tsv(weights.weights, out / f"{cell}_{eid}_mark_weights.tsv")
        write_mark_weights_json(weights, out / f"{cell}_{eid}_mark_weights.json")
        (out / f"{cell}_{eid}_match_qc.json").write_text(
            __import__("json").dumps(assembled.match_qc, indent=2) + "\n",
            encoding="utf-8",
        )
        summary[cell] = {
            "eid": eid,
            "emission_tsv": str(em_path),
            "n_promoter_rows": len(assembled.table.state_ids),
            "mnemonics": list(assembled.table.mnemonics),
            "activities": list(assembled.table.activities),
            "weights": weights.weights,
            "n_active": weights.n_active_states,
            "n_inactive": weights.n_inactive_states,
        }
    (out / "panel_build_summary.json").write_text(
        __import__("json").dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary

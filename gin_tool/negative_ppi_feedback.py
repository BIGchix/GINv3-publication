"""Curated PPI negative-feedback sign fixes (GINv3.11 + subgraph TSVs).

GIN encodes some E3⊣substrate pairs as COP1-style both-negative R→V→R
(``COP1.protein → COP1.protein;TP53.protein → TP53.protein``, negative),
but several canonical ligase/inhibitor pairs were stored as **positive**
complex edges (e.g. MDM2–TP53). SPN treats any negative endpoint as
inhibition, so the wrong sign removes negative feedback from dynamics.

This module flips those existing RVR legs to ``negative`` / ``-1`` and
adds a missing V→substrate.protein negative leg when the ligase→V edge
already exists (COP1 template). Complex-node outputs such as
``MDM2.protein;TP53.protein → MDM2_TP53`` stay positive (binding).
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

# Ligase/inhibitor → virtual, then virtual → substrate.protein
FLIP_TO_NEGATIVE: tuple[tuple[str, str], ...] = (
    # MDM2 E3 ⊣ TP53 / p63 / p73
    ("MDM2.protein", "MDM2.protein;TP53.protein"),
    ("MDM2.protein;TP53.protein", "TP53.protein"),
    ("MDM2.protein", "MDM2.protein;TP63.protein"),
    ("MDM2.protein;TP63.protein", "TP63.protein"),
    ("MDM2.protein", "MDM2.protein;TP73.protein"),
    ("MDM2.protein;TP73.protein", "TP73.protein"),
    ("MDM2.protein", "MDM2.protein;TP73L.protein"),
    ("MDM2.protein;TP73L.protein", "TP73L.protein"),
    # Pirh2 E3 ⊣ p73
    ("RCHY1.protein", "RCHY1.protein;TP73.protein"),
    ("RCHY1.protein;TP73.protein", "TP73.protein"),
    # Siah1 E3 ⊣ β-catenin
    ("SIAH1.protein", "CTNNB1.protein;SIAH1.protein"),
    ("CTNNB1.protein;SIAH1.protein", "CTNNB1.protein"),
    # β-TrCP E3 ⊣ β-catenin / IκBα
    ("BTRC.protein", "BTRC.protein;CTNNB1.protein"),
    ("BTRC.protein;CTNNB1.protein", "CTNNB1.protein"),
    ("BTRC.protein", "BTRC.protein;NFKBIA.protein"),
    ("BTRC.protein;NFKBIA.protein", "NFKBIA.protein"),
    # p14ARF ⊣ MDM2
    ("CDKN2A.protein", "CDKN2A.protein;MDM2.protein"),
    ("CDKN2A.protein;MDM2.protein", "MDM2.protein"),
)

# Add only when missing (V→substrate). Ligase→V is flipped above if present.
ADD_NEGATIVE_IF_MISSING: tuple[tuple[str, str], ...] = (
    ("MDM2.protein;TP73.protein", "TP73.protein"),
    ("BTRC.protein;CTNNB1.protein", "CTNNB1.protein"),
    ("CDKN2A.protein;MDM2.protein", "MDM2.protein"),
)

# Already COP1-style both-negative in GINv3.11 (do not touch).
ALREADY_NEGATIVE_OK: tuple[tuple[str, str], ...] = (
    ("COP1.protein", "COP1.protein;TP53.protein"),
    ("COP1.protein;TP53.protein", "TP53.protein"),
    ("MDM2.protein", "MDM2.protein;CDKN1A.protein"),
    ("MDM2.protein;CDKN1A.protein", "CDKN1A.protein"),
    ("SKP2.protein", "SKP2.protein;CDKN1A.protein"),
    ("SKP2.protein;CDKN1A.protein", "CDKN1A.protein"),
    ("APC_AXIN1_CTNNB1;GSK3B.protein", "CTNNB1.protein"),
)

FLIP_SET = frozenset(FLIP_TO_NEGATIVE)
ADD_SET = frozenset(ADD_NEGATIVE_IF_MISSING)


def _norm_rel(rel: str) -> str:
    r = rel.strip().lower()
    if r in {"negative", "-1", "-1.0"}:
        return "negative"
    if r in {"positive", "1", "1.0", "+1", "+1.0"}:
        return "positive"
    return rel.strip()


def rewrite_sifi_line(
    line: str,
    *,
    flip: frozenset[tuple[str, str]] = FLIP_SET,
) -> tuple[str, tuple[str, str] | None, str]:
    """Return (written_line, (src, tgt) or None, action).

    ``action`` is ``flip``, ``already``, or ``keep``.
    """
    raw = line.rstrip("\n")
    if not raw.strip():
        return (line if line.endswith("\n") else line + "\n"), None, "keep"
    parts = raw.split("\t")
    if len(parts) < 3:
        return raw + "\n", None, "keep"
    src, tgt, rel = parts[0], parts[1], parts[2]
    key = (src, tgt)
    if key in flip:
        if _norm_rel(rel) != "negative":
            parts[2] = "negative"
            return "\t".join(parts[:3]) + "\n", key, "flip"
        return raw + "\n", key, "already"
    return raw + "\n", key, "keep"


def apply_sifi_text(
    lines: Iterable[str],
    *,
    flip: frozenset[tuple[str, str]] = FLIP_SET,
    add: frozenset[tuple[str, str]] = ADD_SET,
) -> tuple[list[str], dict]:
    """Rewrite sifi lines (``src\\ttgt\\tpositive|negative``)."""
    out: list[str] = []
    seen: set[tuple[str, str]] = set()
    n_flip = n_already = n_keep = 0
    for raw in lines:
        written, key, action = rewrite_sifi_line(raw, flip=flip)
        out.append(written)
        if key:
            seen.add(key)
        if action == "flip":
            n_flip += 1
        elif action == "already":
            n_already += 1
        else:
            n_keep += 1
    n_add = 0
    for src, tgt in sorted(add):
        if (src, tgt) in seen:
            continue
        out.append(f"{src}\t{tgt}\tnegative\n")
        n_add += 1
    stats = {
        "n_flipped": n_flip,
        "n_already_negative": n_already,
        "n_added": n_add,
        "n_unchanged": n_keep,
    }
    return out, stats


def apply_numeric_tsv(
    lines: Iterable[str],
    *,
    flip: frozenset[tuple[str, str]] = FLIP_SET,
    add: frozenset[tuple[str, str]] = ADD_SET,
    present_nodes: set[str] | None = None,
) -> tuple[list[str], dict]:
    """Rewrite subgraph TSV (``src\\ttgt\\t±1``). Skip adds unless both nodes exist."""
    out: list[str] = []
    seen: set[tuple[str, str]] = set()
    nodes: set[str] = set(present_nodes or ())
    n_flip = n_already = n_keep = 0
    for raw in lines:
        line = raw.rstrip("\n")
        if not line.strip():
            out.append(line + "\n")
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            out.append(line + "\n")
            continue
        src, tgt, wraw = parts[0], parts[1], parts[2]
        nodes.add(src)
        nodes.add(tgt)
        key = (src, tgt)
        seen.add(key)
        try:
            w = float(wraw)
        except ValueError:
            n_keep += 1
            out.append(line + "\n")
            continue
        if key in flip:
            if w >= 0:
                parts[2] = "-1"
                n_flip += 1
            else:
                n_already += 1
            out.append("\t".join(parts[:3]) + "\n")
        else:
            n_keep += 1
            out.append(line + "\n")
    n_add = 0
    for src, tgt in sorted(add):
        if (src, tgt) in seen:
            continue
        if src not in nodes or tgt not in nodes:
            continue
        out.append(f"{src}\t{tgt}\t-1\n")
        n_add += 1
    stats = {
        "n_flipped": n_flip,
        "n_already_negative": n_already,
        "n_added": n_add,
        "n_unchanged": n_keep,
    }
    return out, stats

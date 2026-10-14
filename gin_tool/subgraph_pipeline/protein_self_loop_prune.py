"""Drop GIN homodimer protein self-loops (R→V→R) before SPN export.

Pattern (e.g. TP53 / MDM2):
  GENE.protein → GENE.protein;GENE.protein → GENE.protein
These mint/copy protein tokens and distort feedback dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass


def is_protein_homodimer_self_loop_edge(src: str, tgt: str) -> bool:
    """True for both legs of a monomer homodimer R-V-R self-loop."""
    if src.endswith(".protein") and ";" not in src:
        gene = src[: -len(".protein")]
        return tgt == f"{gene}.protein;{gene}.protein"
    if tgt.endswith(".protein") and ";" not in tgt:
        gene = tgt[: -len(".protein")]
        return src == f"{gene}.protein;{gene}.protein"
    return False


@dataclass(frozen=True)
class ProteinSelfLoopPruneStats:
    n_in: int
    n_out: int
    n_removed: int
    removed_pairs: tuple[tuple[str, str], ...]

    def to_dict(self) -> dict:
        return {
            "n_in": self.n_in,
            "n_out": self.n_out,
            "n_removed": self.n_removed,
            "removed_pairs": [list(p) for p in self.removed_pairs],
        }


def drop_protein_self_loops(
    edges: list[tuple[str, str, float]],
) -> tuple[list[tuple[str, str, float]], ProteinSelfLoopPruneStats]:
    """Return filtered edge list and prune stats."""
    kept: list[tuple[str, str, float]] = []
    removed: list[tuple[str, str]] = []
    for src, tgt, w in edges:
        if is_protein_homodimer_self_loop_edge(src, tgt):
            removed.append((src, tgt))
        else:
            kept.append((src, tgt, w))
    stats = ProteinSelfLoopPruneStats(
        n_in=len(edges),
        n_out=len(kept),
        n_removed=len(removed),
        removed_pairs=tuple(removed),
    )
    return kept, stats

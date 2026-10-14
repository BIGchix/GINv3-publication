"""Export subgraph for Cytoscape / external tools."""

from __future__ import annotations

from pathlib import Path

from scipy import sparse

from gin_tool.io import REL_NEG, REL_POS
from gin_tool.mapping import gin_entity_layer


def export_cytoscape_sif(
    path: Path,
    sub_sym: sparse.csr_matrix,
    sub_dir: sparse.csr_matrix,
    names: list[str],
    scores: np.ndarray,
    seed_indices: set[int],
    node_indices: set[int],
    module_by_node: dict[int, int] | None,
    *,
    edge_mode: str = "none",
) -> None:
    """
    Write a simplified interaction file with node attributes as comments.

    Format: tab-separated edge list plus ``#NODE`` attribute lines for Cytoscape import.
    """
    idx_sorted = sorted(node_indices)
    old_to_local = {g: li for li, g in enumerate(idx_sorted)}

    lines: list[str] = [
        "#source\ttarget\tinteraction\tsign\tweight",
    ]
    sym = sub_sym.tocsr()
    for i in idx_sorted:
        lo, hi = sym.indptr[i], sym.indptr[i + 1]
        for j in sym.indices[lo:hi]:
            j = int(j)
            if j not in node_indices or j <= i:
                continue
            w = float(sym[i, j])
            sign = 0
            if edge_mode == "sign":
                sign = _dominant_sign(sub_dir, i, j)
            lines.append(
                f"{names[i]}\t{names[j]}\tundirected\t{sign}\t{w:g}"
            )

    for i in idx_sorted:
        mid = module_by_node.get(i, -1) if module_by_node else -1
        seed = 1 if i in seed_indices else 0
        lines.append(
            f"#NODE\t{names[i]}\trwr_score={scores[i]:g}\tlayer={gin_entity_layer(names[i])}"
            f"\tmodule_id={mid}\tis_seed={seed}"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _dominant_sign(sub_dir: sparse.csr_matrix, i: int, j: int) -> int:
    s = 0
    w_ij = float(sub_dir[i, j]) if sub_dir.shape[0] else 0.0
    w_ji = float(sub_dir[j, i]) if sub_dir.shape[0] else 0.0
    if w_ij > 0:
        s += REL_POS
    elif w_ij < 0:
        s += REL_NEG
    if w_ji > 0:
        s += REL_POS
    elif w_ji < 0:
        s += REL_NEG
    if s > 0:
        return REL_POS
    if s < 0:
        return REL_NEG
    return 0

"""Simplified SVG schematic: modules (rounded boxes) and information-flow polylines."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
from xml.sax.saxutils import escape


def _module_box_center(mi: int, n_mods: int, width: float, height: float) -> tuple[float, float]:
    margin_x, margin_y = 48.0, 40.0
    usable = width - 2 * margin_x
    stride = usable / max(n_mods, 1)
    cx = margin_x + mi * stride + stride / 2.0
    cy = height * 0.74
    return cx, cy


def path_module_sequence(path: Sequence[int], modules: Sequence[frozenset[int]]) -> list[int]:
    """Compress consecutive duplicates of module ids encountered along a vertex path."""
    seq: list[int | None] = []
    for v in path:
        mid: int | None = None
        for i, mod in enumerate(modules):
            if v in mod:
                mid = i
                break
        seq.append(mid)
    out: list[int] = []
    for x in seq:
        if x is None:
            continue
        if not out or out[-1] != x:
            out.append(x)
    return out


def render_schematic_svg(
    out_path: str | Path,
    modules: Sequence[frozenset[int]],
    paths: Sequence[Sequence[int]],
    sub_names: Sequence[str],
    interface_nodes: set[int],
    *,
    width: float = 960.0,
    height: float = 520.0,
    max_flows_drawn: int = 14,
) -> None:
    """
    Draw rounded rectangles for each module and colored polylines for the first
    ``max_flows_drawn`` ranked flows (module-to-module spine). Only **interface**
    real nodes get text labels (others are abstracted into the spine).
    """
    out_path = Path(out_path)
    W, H = width, height
    n_mod = len(modules)
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}">',
        f'<rect width="{W:.0f}" height="{H:.0f}" fill="#fafafa"/>',
        '<text x="24" y="32" font-size="18" font-family="Helvetica,Arial,sans-serif" fill="#222">GIN schematic: modules (boxes) vs information flows (lines)</text>',
    ]

    for mi in range(n_mod):
        cx, cy = _module_box_center(mi, n_mod, W, H)
        bw, bh = 130.0, 64.0
        parts.append(
            f'<rect x="{cx - bw/2:.1f}" y="{cy - bh/2:.1f}" width="{bw:.0f}" height="{bh:.0f}" '
            f'rx="14" fill="#e8f0fe" stroke="#1a5276" stroke-width="1.6"/>'
        )
        parts.append(
            f'<text x="{cx:.1f}" y="{cy + 5:.1f}" text-anchor="middle" font-size="13" '
            f'font-family="Helvetica,Arial,sans-serif" fill="#1a5276">Module {mi}</text>'
        )

    palette = ["#c0392b", "#1e8449", "#7d3c98", "#ca6f1e", "#2471a3", "#117864", "#884ea0"]

    for pi, path in enumerate(paths[:max_flows_drawn]):
        spine = path_module_sequence(path, modules)
        if len(spine) < 2:
            continue
        pts: list[str] = []
        for mid in spine:
            cx, cy = _module_box_center(mid, n_mod, W, H)
            pts.append(f"{cx:.1f},{cy - 110:.1f}")
        color = palette[pi % len(palette)]
        parts.append(
            f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
            f'stroke-width="2.2" stroke-linejoin="round" opacity="0.9"/>'
        )
        parts.append(
            f'<text x="24" y="{58 + pi * 16:.1f}" font-size="11" font-family="Helvetica,Arial,sans-serif" fill="{color}">Flow {pi}</text>'
        )

    lx, ly = W * 0.52, 96.0
    parts.append(
        f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="12" font-weight="bold" '
        f'font-family="Helvetica,Arial,sans-serif" fill="#333">Interface molecules</text>'
    )
    iy = ly + 18.0
    for vi in sorted(interface_nodes):
        label = escape(sub_names[vi][:80])
        parts.append(
            f'<text x="{lx:.1f}" y="{iy:.1f}" font-size="10" font-family="Helvetica,Arial,sans-serif" fill="#444">{label}</text>'
        )
        iy += 14.0
        if iy > H - 24:
            parts.append(
                f'<text x="{lx:.1f}" y="{iy:.1f}" font-size="10" font-family="Helvetica,Arial,sans-serif" fill="#777">…</text>'
            )
            break

    parts.append("</svg>")
    out_path.write_text("\n".join(parts), encoding="utf-8")

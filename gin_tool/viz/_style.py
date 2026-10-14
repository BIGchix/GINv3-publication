"""Shared colors and matplotlib defaults for GIN visualizations."""

from __future__ import annotations

LAYER_COLORS: dict[str, str] = {
    "rna": "#2980b9",
    "protein": "#c0392b",
    "promoter_chromatin": "#8e44ad",
    "metabolite": "#27ae60",
    "virtual": "#95a5a6",
    "other": "#7f8c8d",
}

MODULE_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]

SIGN_EDGE_COLORS = {1: "#27ae60", -1: "#c0392b", 0: "#bdc3c7"}

DPI = 120

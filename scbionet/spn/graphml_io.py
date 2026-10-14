"""Parse mEPN GraphML produced by GINtoSPN (yEd / BioLayout schema)."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

from scbionet.spn.model import SPNInteraction, SPNModel

_NS = {
    "g": "http://graphml.graphdrawing.org/xmlns",
    "y": "http://www.yworks.com/xml/graphml",
}

_TRANSITION_LABELS = {"C", "B", "TR", "TL", "/", "T", "D"}
_TRANSITION_COLORS = {"#9999FF", "#FFFF99", "#99CCFF"}


def _local(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


def _node_label(node_el: ET.Element) -> str:
    for child in node_el.iter():
        loc = _local(child.tag)
        if loc == "NodeLabel" and child.text:
            return child.text.strip()
        if loc == "data" and child.attrib.get("key") == "d0" and child.text:
            return child.text.strip()
    return ""


def _node_shape(node_el: ET.Element) -> str:
    for child in node_el.iter():
        if _local(child.tag) == "Shape":
            return child.attrib.get("type", "")
    return ""


def _node_fill(node_el: ET.Element) -> str:
    for child in node_el.iter():
        if _local(child.tag) == "Fill":
            return child.attrib.get("color", "").upper()
    return ""


def _edge_arrow_target(edge_el: ET.Element) -> str:
    for child in edge_el.iter():
        if _local(child.tag) == "Arrows":
            return child.attrib.get("target", "standard")
    return "standard"


def _edge_weight(edge_el: ET.Element) -> float:
    for child in edge_el.iter():
        if _local(child.tag) == "EdgeLabel" and child.text:
            try:
                return float(child.text.strip())
            except ValueError:
                return 1.0
    return 1.0


def _is_transition(label: str, shape: str, fill: str) -> bool:
    if label in _TRANSITION_LABELS:
        return True
    if ";" in label:
        return True
    if fill in _TRANSITION_COLORS:
        return True
    if shape == "ellipse" and label in _TRANSITION_LABELS:
        return True
    return False


def load_mepn_graphml(path: str | Path, *, add_self_inhibition: bool = True) -> SPNModel:
    """Load GINtoSPN / BioLayout mEPN GraphML into an :class:`SPNModel`."""
    root = ET.parse(path).getroot()
    graph = root.find(".//g:graph", _NS)
    if graph is None:
        graph = root.find(".//{http://graphml.graphdrawing.org/xmlns}graph")
    if graph is None:
        raise ValueError(f"no graph element in {path}")

    nodes: dict[str, dict] = {}
    for node_el in graph.findall("g:node", _NS) or graph.findall(
        "{http://graphml.graphdrawing.org/xmlns}node"
    ):
        nid = node_el.attrib["id"]
        label = _node_label(node_el)
        nodes[nid] = {
            "label": label or nid,
            "is_transition": _is_transition(
                label, _node_shape(node_el), _node_fill(node_el)
            ),
        }

    edges: list[dict] = []
    for edge_el in graph.findall("g:edge", _NS) or graph.findall(
        "{http://graphml.graphdrawing.org/xmlns}edge"
    ):
        edges.append(
            {
                "source": edge_el.attrib["source"],
                "target": edge_el.attrib["target"],
                "inhibitory": _edge_arrow_target(edge_el) == "t_shape",
                "weight": _edge_weight(edge_el),
            }
        )

    transitions = {nid for nid, meta in nodes.items() if meta["is_transition"]}
    places = {nid: meta["label"] for nid, meta in nodes.items() if nid not in transitions}

    interactions: list[SPNInteraction] = []
    for tid in transitions:
        in_places: list[tuple[str, float]] = []
        out_places: list[tuple[str, bool, float]] = []
        for e in edges:
            if e["target"] == tid and e["source"] in places:
                in_places.append((e["source"], e["weight"]))
            elif e["source"] == tid and e["target"] in places:
                out_places.append((e["target"], e["inhibitory"], e["weight"]))

        for src, w_in in in_places:
            src_label = places[src]
            for tgt, inhibitory, w_out in out_places:
                tgt_label = places[tgt]
                kind = "inhibition" if inhibitory else "activation"
                interactions.append(
                    SPNInteraction(
                        source=src_label,
                        target=tgt_label,
                        kind=kind,
                        transition_id=tid,
                        weight=w_in * w_out,
                    )
                )

        # transition with only inputs → treat as inhibition from each input to itself? skip
        if not out_places and len(in_places) == 2:
            # common inhibition topology: regulator → T ← target
            for tgt_pid, _ in in_places:
                for reg_pid, w in in_places:
                    if reg_pid == tgt_pid:
                        continue
                    interactions.append(
                        SPNInteraction(
                            source=places[reg_pid],
                            target=places[tgt_pid],
                            kind="inhibition",
                            transition_id=tid,
                            weight=w,
                        )
                    )

    place_labels = list(dict.fromkeys(interactions[i].source for i in range(len(interactions))))
    place_labels += [e.target for e in interactions if e.target not in place_labels]
    place_labels = sorted(set(place_labels))

    model = SPNModel(
        places=place_labels,
        interactions=interactions,
        place_labels={p: p for p in place_labels},
    )
    if add_self_inhibition:
        model.add_self_inhibition_loops()
    return model


def export_topology_json(model: SPNModel, path: str | Path) -> None:
    """Write place/interaction topology for SPN dashboard consumers."""
    import json

    payload = {
        "places": list(model.places),
        "interactions": [
            {
                "source": e.source,
                "target": e.target,
                "kind": e.kind,
                "weight": e.weight,
                "transition_id": e.transition_id,
            }
            for e in model.interactions
        ],
    }
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")

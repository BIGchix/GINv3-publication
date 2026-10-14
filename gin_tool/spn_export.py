"""Export rehydrated GIN signaling subgraph to mEPN GraphML for BioLayout SPN."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Mapping

from gin_tool.mapping import is_virtual_gin_node
from gin_tool.regulatory_filter import is_dna_element
from gin_tool.spn_rehydrate import RehydratedSubgraph

_GRAPHML_NS = "http://graphml.graphdrawing.org/xmlns"
ET.register_namespace("", _GRAPHML_NS)


class SpnPlaceKind(str, Enum):
    PROTEIN = "Protein"
    PROTEIN_COMPLEX = "Protein_complex"
    RNA = "RNA"
    DNA_REGION = "DNA_region"
    PROMOTER = "Promoter"
    CHROMATIN_LOCUS = "chromatin_locus"
    HISTONE_MARK = "histone_mark"
    METABOLITE = "Metabolite"
    OTHER = "Other"


class SpnTransitionKind(str, Enum):
    BINDING = "Binding"
    CATALYSIS = "Catalysis"
    TRANSCRIPTION = "Transcription"
    TRANSLATION = "Translation"
    GIN_VIRTUAL = "GIN_virtual"
    DEGRADATION = "Degradation"


# mEPN display hints (BioLayout / yEd compatible keys)
_PLACE_STYLE: dict[SpnPlaceKind, dict[str, str]] = {
    SpnPlaceKind.PROTEIN: {"shape": "ellipse", "color": "#FFCCCC"},
    SpnPlaceKind.PROTEIN_COMPLEX: {"shape": "ellipse", "color": "#FFB366"},
    SpnPlaceKind.RNA: {"shape": "ellipse", "color": "#CCCCFF"},
    SpnPlaceKind.DNA_REGION: {"shape": "rectangle", "color": "#CCFFCC"},
    SpnPlaceKind.PROMOTER: {"shape": "rectangle", "color": "#CCFFCC"},
    SpnPlaceKind.CHROMATIN_LOCUS: {"shape": "rectangle", "color": "#E6FFE6"},
    SpnPlaceKind.HISTONE_MARK: {"shape": "rectangle", "color": "#D9FFD9"},
    SpnPlaceKind.METABOLITE: {"shape": "hexagon", "color": "#FFFFCC"},
    SpnPlaceKind.OTHER: {"shape": "ellipse", "color": "#EEEEEE"},
}

_TRANSITION_STYLE: dict[SpnTransitionKind, dict[str, str]] = {
    SpnTransitionKind.BINDING: {"shape": "rectangle", "color": "#999999"},
    SpnTransitionKind.CATALYSIS: {"shape": "rectangle", "color": "#666666"},
    SpnTransitionKind.TRANSCRIPTION: {"shape": "rectangle", "color": "#88AA88"},
    SpnTransitionKind.TRANSLATION: {"shape": "rectangle", "color": "#8888AA"},
    SpnTransitionKind.GIN_VIRTUAL: {"shape": "rectangle", "color": "#AAAAAA"},
    SpnTransitionKind.DEGRADATION: {"shape": "rectangle", "color": "#666666"},
}


# CAS registry IDs are ``CAS\d+-\d+-\d`` (e.g. CAS50-99-7). Do NOT treat gene
# symbols like CASP3 / CAST / CASK as metabolites (prefix collision).
_CAS_REGISTRY_RE = re.compile(r"^CAS\d{2,7}-\d{2}-\d(?:$|[_.;])", re.I)


def _is_metabolite(name: str) -> bool:
    if name.startswith(("CHEBI", "SID", "CID", "HMDB", "KEGGC")):
        return True
    if _CAS_REGISTRY_RE.match(name):
        return True
    lower = name.lower()
    return "chebi" in lower or lower.endswith(".met")


def _is_rna(name: str) -> bool:
    lower = name.lower()
    return lower.endswith(".rna") or "mir-" in lower or "mirna" in lower


def is_mirna_place(name: str) -> bool:
    """True for miRNA / MIR host transcript places (not generic mRNA)."""
    lower = (name or "").lower()
    if "mirna" in lower:
        return True
    # MIR3121.RNA, hsa-miR-338-3p.RNA, hsa-mir-21, …
    if lower.startswith("hsa-mir") or lower.startswith("mir"):
        return True
    if "mir-" in lower or ".mir" in lower:
        return True
    # Gene-style MIR###.RNA (numeric miRNA loci)
    base = lower.split(".", 1)[0]
    if base.startswith("mir") and any(ch.isdigit() for ch in base):
        return True
    return False


def _is_promoter(name: str) -> bool:
    return ".promoter" in name.lower()


def _is_chromatin_locus(name: str) -> bool:
    return name.startswith("chr") or bool(re.match(r"^chr\d", name, re.I))


def _is_histone_mark(name: str) -> bool:
    upper = name.upper()
    return any(h in upper for h in ("H3K", "H4K", "H2AK", "H2BK", "H3AC", "H4AC"))


def _is_protein_complex(name: str) -> bool:
    return "_" in name and ";" not in name and "." not in name and not _is_metabolite(name)


def classify_place(name: str, layer_map: Mapping[str, str] | None = None) -> SpnPlaceKind:
    """Map a real GIN node name to an mEPN place kind."""
    if is_virtual_gin_node(name):
        raise ValueError(f"expected real node, got virtual: {name!r}")

    layer = (layer_map or {}).get(name)
    if layer == "metabolite" or _is_metabolite(name):
        return SpnPlaceKind.METABOLITE
    if layer in ("histone_mark", "hme") or _is_histone_mark(name):
        return SpnPlaceKind.HISTONE_MARK
    if _is_promoter(name):
        return SpnPlaceKind.PROMOTER
    if layer == "promoter_chromatin" or _is_chromatin_locus(name) or is_dna_element(name):
        if _is_chromatin_locus(name):
            return SpnPlaceKind.CHROMATIN_LOCUS
        return SpnPlaceKind.DNA_REGION
    if layer in ("rna", "mirna", "lncrna") or _is_rna(name):
        return SpnPlaceKind.RNA
    if layer == "complex" or _is_protein_complex(name):
        return SpnPlaceKind.PROTEIN_COMPLEX
    return SpnPlaceKind.PROTEIN


def classify_virtual_transition(name: str) -> SpnTransitionKind:
    """Classify a GIN virtual node as an mEPN transition type."""
    lower = name.lower()
    if "translation" in lower:
        return SpnTransitionKind.TRANSLATION
    # Real transcription reactions are ``*;transcription``. Do NOT treat
    # ``GENE.promoter;state`` as transcription — that dropped state→active arcs
    # via allow_petri_arc (transcription may only produce RNA).
    if "transcription" in lower:
        return SpnTransitionKind.TRANSCRIPTION
    if ".promoter;state" in lower or lower.endswith("promoter;state"):
        return SpnTransitionKind.GIN_VIRTUAL
    if _is_metabolite(name) or "chebi" in lower:
        return SpnTransitionKind.CATALYSIS
    return SpnTransitionKind.GIN_VIRTUAL


def is_markable_seed_place(name: str) -> bool:
    """Legacy protein-only mark filter (``*.protein``)."""
    return name.lower().endswith(".protein")


def is_degradation_sink_place(name: str) -> bool:
    low = name.lower()
    return (
        low.startswith("__sink")
        or low.startswith("__spn_sink__")
        or low in {"__sink_rna__", "__sink_protein__", "__sink_metabolite__"}
    )


# Token injection modes for ancestor-pruned context SPNs.
# After prune_mode=ancestors the digraph mostly flows *into* seeds, so marking only
# seed proteins leaves upstream places starved.  "祖先-only" token sources =
# Petri places with in-degree 0 (frontier sources), not every place in the net.
MARK_MODES = (
    "seed_protein",
    "source_protein",
    "ancestor_protein",
    "ancestor_all",  # alias: in-degree-0 places (any type)
    "all_places",  # every non-sink place (legacy full maintain)
)


def _places_with_incoming(net: "SpnPetriNet") -> set[str]:
    return {tgt for src, tgt in net.arcs if tgt in net.places}


def resolve_mark_places(
    net: "SpnPetriNet",
    *,
    seed_gin_names: set[str] | None = None,
    mark_mode: str = "seed_protein",
) -> set[str]:
    """Choose places that receive initial / maintain tokens.

    * ``seed_protein`` — seed anchors ∩ ``*.protein`` (legacy).
    * ``source_protein`` — ``*.protein`` with no incoming Petri arc.
    * ``ancestor_protein`` — every ``*.protein`` place in the net.
    * ``ancestor_all`` — every non-sink place with Petri in-degree 0
      (祖先-only sources; subgraph is already ancestor-pruned).
    * ``all_places`` — every non-sink place (legacy full maintain).
    """
    mode = (mark_mode or "seed_protein").strip().lower()
    if mode not in MARK_MODES:
        raise ValueError(f"unknown mark_mode {mark_mode!r}; expected one of {MARK_MODES}")
    proteins = {p for p in net.places if is_markable_seed_place(p)}
    has_in = _places_with_incoming(net)
    if mode == "all_places":
        return {p for p in net.places if not is_degradation_sink_place(p)}
    if mode == "ancestor_all":
        return {
            p
            for p in net.places
            if p not in has_in and not is_degradation_sink_place(p)
        }
    if mode == "ancestor_protein":
        return proteins
    if mode == "source_protein":
        return {p for p in proteins if p not in has_in}
    seeds = seed_gin_names or set()
    return {p for p in proteins if p in seeds}


def apply_place_tokens(
    net: "SpnPetriNet",
    *,
    seed_gin_names: set[str] | None = None,
    mark_mode: str = "seed_protein",
) -> set[str]:
    """Set ``place_tokens`` according to ``mark_mode``; return marked names."""
    marked = resolve_mark_places(net, seed_gin_names=seed_gin_names, mark_mode=mark_mode)
    for p in net.places:
        net.place_tokens[p] = 1 if p in marked else 0
    return marked


def allow_petri_arc(src: str, tgt: str) -> bool:
    """
    Drop biologically invalid arcs when materializing GIN bipartite edges as Petri arcs.

    Transcription transitions produce RNA, not promoters/DNA (GIN often stores R↔V
    both ways around ``*;transcription``).
    """
    if is_virtual_gin_node(src) and not is_virtual_gin_node(tgt):
        if classify_virtual_transition(src) == SpnTransitionKind.TRANSCRIPTION:
            return _is_rna(tgt)
    return True


def infer_synthetic_transition(src: str, tgt: str, layer_map: Mapping[str, str] | None = None) -> SpnTransitionKind:
    """GINtoSPN Fig.5D-style synthetic transition when no GIN virtual intermediate exists."""
    sk = classify_place(src, layer_map)
    tk = classify_place(tgt, layer_map)
    if sk in (
        SpnPlaceKind.DNA_REGION,
        SpnPlaceKind.PROMOTER,
        SpnPlaceKind.CHROMATIN_LOCUS,
        SpnPlaceKind.HISTONE_MARK,
    ) and tk == SpnPlaceKind.RNA:
        return SpnTransitionKind.TRANSCRIPTION
    if sk == SpnPlaceKind.RNA and tk in (SpnPlaceKind.PROTEIN, SpnPlaceKind.PROTEIN_COMPLEX):
        return SpnTransitionKind.TRANSLATION
    if sk in (SpnPlaceKind.PROTEIN, SpnPlaceKind.PROTEIN_COMPLEX) and tk in (
        SpnPlaceKind.PROTEIN,
        SpnPlaceKind.PROTEIN_COMPLEX,
    ):
        return SpnTransitionKind.BINDING
    if sk == SpnPlaceKind.METABOLITE or tk == SpnPlaceKind.METABOLITE:
        return SpnTransitionKind.CATALYSIS
    return SpnTransitionKind.BINDING


@dataclass
class SpnPetriNet:
    places: dict[str, SpnPlaceKind] = field(default_factory=dict)
    transitions: dict[str, SpnTransitionKind] = field(default_factory=dict)
    arcs: list[tuple[str, str]] = field(default_factory=list)
    # Signed GIN weights on place↔transition arcs (default 1.0 if missing).
    # Needed so histone⊣promoter (weight=-1) becomes SPN inhibition.
    arc_weights: dict[tuple[str, str], float] = field(default_factory=dict)
    place_tokens: dict[str, int] = field(default_factory=dict)
    synthetic_transitions: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        place_counts: dict[str, int] = {}
        for k in self.places.values():
            place_counts[k.value] = place_counts.get(k.value, 0) + 1
        trans_counts: dict[str, int] = {}
        for k in self.transitions.values():
            trans_counts[k.value] = trans_counts.get(k.value, 0) + 1
        return {
            "n_places": len(self.places),
            "n_transitions": len(self.transitions),
            "n_arcs": len(self.arcs),
            "place_kind_counts": place_counts,
            "transition_kind_counts": trans_counts,
            "n_synthetic_transitions": len(self.synthetic_transitions),
            "n_marked_places": sum(1 for v in self.place_tokens.values() if v > 0),
        }


def build_petri_net(
    rehydrated: RehydratedSubgraph,
    *,
    layer_map: Mapping[str, str] | None = None,
    seed_gin_names: set[str] | None = None,
) -> SpnPetriNet:
    """Convert rehydrated subgraph to a place/transition Petri net."""
    net = SpnPetriNet()
    seed_gin_names = seed_gin_names or set()

    for name in rehydrated.real_nodes:
        net.places[name] = classify_place(name, layer_map)
        if name in seed_gin_names and is_markable_seed_place(name):
            net.place_tokens[name] = 1
        else:
            net.place_tokens.setdefault(name, 0)

    synth_counter = 0
    arc_set: set[tuple[str, str]] = set()

    def _add_arc(u: str, v: str) -> None:
        if not allow_petri_arc(u, v):
            return
        if (u, v) not in arc_set:
            arc_set.add((u, v))
            net.arcs.append((u, v))

    def _new_synthetic(kind: SpnTransitionKind) -> str:
        nonlocal synth_counter
        synth_counter += 1
        tid = f"synth_{kind.value}_{synth_counter}"
        net.synthetic_transitions.append(tid)
        net.transitions[tid] = kind
        return tid

    for edge in rehydrated.logical_edges:
        src, tgt = edge.source, edge.target
        if src not in net.places or tgt not in net.places:
            continue
        if edge.via_virtual:
            v = edge.via_virtual
            net.transitions[v] = classify_virtual_transition(v)
            _add_arc(src, v)
            _add_arc(v, tgt)
        else:
            kind = infer_synthetic_transition(src, tgt, layer_map)
            tid = _new_synthetic(kind)
            _add_arc(src, tid)
            _add_arc(tid, tgt)

    for v in rehydrated.virtual_nodes:
        if v not in net.transitions:
            net.transitions[v] = classify_virtual_transition(v)

    return net


def _xml_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_mepn_graphml(path: Path | str, net: SpnPetriNet) -> None:
    """
    Write mEPN-styled GraphML for BioLayout Express3D.

    Schema follows GINtoSPN conventions: ``element_type`` (place|transition),
    ``mepn_kind``, display shape/color, optional ``tokens`` on places.
    """
    path = Path(path)
    root = ET.Element(f"{{{_GRAPHML_NS}}}graphml")

    keys = [
        ("k_element_type", "node", "element_type", "string"),
        ("k_mepn_kind", "node", "mepn_kind", "string"),
        ("k_label", "node", "label", "string"),
        ("k_shape", "node", "shape", "string"),
        ("k_color", "node", "color", "string"),
        ("k_tokens", "node", "tokens", "int"),
        ("k_gin_name", "node", "gin_name", "string"),
        ("k_edge_sign", "edge", "sign", "double"),
    ]
    key_ids: dict[str, str] = {}
    for kid, for_what, attr_name, attr_type in keys:
        key_el = ET.SubElement(root, f"{{{_GRAPHML_NS}}}key")
        key_el.set("id", kid)
        key_el.set("for", for_what)
        key_el.set("attr.name", attr_name)
        key_el.set("attr.type", attr_type)
        key_ids[attr_name] = kid

    graph = ET.SubElement(root, f"{{{_GRAPHML_NS}}}graph")
    graph.set("edgedefault", "directed")

    node_id_map: dict[str, str] = {}

    def _add_node(nid: str, element_type: str, mepn_kind: str, label: str, style: dict[str, str], tokens: int = 0, gin_name: str = "") -> None:
        node_id_map[nid] = nid
        node = ET.SubElement(graph, f"{{{_GRAPHML_NS}}}node")
        node.set("id", nid)
        for attr, val in (
            ("element_type", element_type),
            ("mepn_kind", mepn_kind),
            ("label", label),
            ("shape", style.get("shape", "ellipse")),
            ("color", style.get("color", "#CCCCCC")),
            ("tokens", str(tokens)),
            ("gin_name", gin_name or label),
        ):
            data = ET.SubElement(node, f"{{{_GRAPHML_NS}}}data")
            data.set("key", key_ids[attr])
            data.text = _xml_escape(str(val))

    for i, (name, kind) in enumerate(sorted(net.places.items())):
        pid = f"p{i}"
        node_id_map[name] = pid
        style = _PLACE_STYLE.get(kind, _PLACE_STYLE[SpnPlaceKind.OTHER])
        _add_node(
            pid,
            "place",
            kind.value,
            name,
            style,
            tokens=net.place_tokens.get(name, 0),
            gin_name=name,
        )

    for i, (name, kind) in enumerate(sorted(net.transitions.items())):
        tid = f"t{i}"
        node_id_map[name] = tid
        style = _TRANSITION_STYLE.get(kind, _TRANSITION_STYLE[SpnTransitionKind.GIN_VIRTUAL])
        display = name if len(name) <= 40 else name[:37] + "..."
        _add_node(tid, "transition", kind.value, display, style, gin_name=name)

    for ei, (src, tgt) in enumerate(net.arcs):
        if src not in node_id_map or tgt not in node_id_map:
            continue
        edge = ET.SubElement(graph, f"{{{_GRAPHML_NS}}}edge")
        edge.set("id", f"e{ei}")
        edge.set("source", node_id_map[src])
        edge.set("target", node_id_map[tgt])
        data = ET.SubElement(edge, f"{{{_GRAPHML_NS}}}data")
        data.set("key", key_ids["sign"])
        data.text = "1.0"

    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(path, encoding="UTF-8", xml_declaration=True)


@dataclass
class SpnExportResult:
    petri_net: SpnPetriNet
    rehydrated: RehydratedSubgraph
    graphml_path: Path | None = None

    def qc_dict(self) -> dict:
        amb_examples = [
            {
                "source": e.source,
                "target": e.target,
                "via": e.via_virtual,
                "alternatives": e.ambiguous_virtuals[:5],
            }
            for e in self.rehydrated.logical_edges
            if e.ambiguous_virtuals
        ][:20]
        unmapped = [
            n
            for n in self.rehydrated.real_nodes
            if classify_place(n) == SpnPlaceKind.OTHER
        ]
        return {
            "rehydrate": self.rehydrated.to_dict(),
            "petri_net": self.petri_net.to_dict(),
            "ambiguous_edge_examples": amb_examples,
            "unclassified_place_nodes": sorted(unmapped)[:50],
            "n_unclassified_places": len(unmapped),
        }


def export_signaling_to_graphml(
    compressed_edges: set[tuple[str, str]],
    real_nodes: set[str],
    ed,
    out_graphml: Path | str,
    *,
    layer_map: Mapping[str, str] | None = None,
    seed_gin_names: set[str] | None = None,
    qc_json_path: Path | str | None = None,
) -> SpnExportResult:
    """End-to-end: rehydrate compressed signaling subgraph and write GraphML + QC."""
    from gin_tool.spn_rehydrate import rehydrate_compressed_edges

    rehydrated = rehydrate_compressed_edges(compressed_edges, ed, real_nodes=real_nodes)
    net = build_petri_net(rehydrated, layer_map=layer_map, seed_gin_names=seed_gin_names)
    out_graphml = Path(out_graphml)
    write_mepn_graphml(out_graphml, net)
    result = SpnExportResult(petri_net=net, rehydrated=rehydrated, graphml_path=out_graphml)
    if qc_json_path is not None:
        Path(qc_json_path).write_text(json.dumps(result.qc_dict(), indent=2), encoding="utf-8")
    return result

"""Entity ID mapping and QC for aligning omics features to GIN node names."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence


def is_virtual_gin_node(name: str) -> bool:
    """GIN virtual / transition nodes contain ';' in the notebook convention."""
    return ";" in name


def gin_entity_layer(node: str) -> str:
    """Infer coarse layer from GIN node string suffix (best-effort heuristic)."""
    if is_virtual_gin_node(node):
        return "virtual"
    lower = node.lower()
    if lower.endswith(".rna"):
        return "rna"
    if lower.endswith(".protein"):
        return "protein"
    if ".promoter." in lower or lower.endswith(".promoter"):
        return "promoter_chromatin"
    if "metabol" in lower or lower.endswith(".met"):
        return "metabolite"
    return "other"


def load_layer_map(path: str | Path) -> dict[str, str]:
    """
    Load ``gin_node\\tlayer`` (header optional) for overrides to ``assign_entity_layer``.

    Lines starting with ``#`` are ignored. Duplicate gin_node: last wins.
    """
    p = Path(path)
    out: dict[str, str] = {}
    with p.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            node, layer = parts[0].strip(), parts[1].strip()
            if node and layer:
                out[node] = layer
    return out


def assign_entity_layer(node: str, layer_map: Mapping[str, str] | None = None) -> str:
    """Layer tag for a GIN node; ``layer_map`` overrides heuristic when present."""
    if layer_map and node in layer_map:
        return str(layer_map[node])
    return gin_entity_layer(node)


def real_entity_indices(index_to_name: Sequence[str]) -> list[int]:
    """Indices of non-virtual nodes (eligible as null random seeds)."""
    return [i for i, n in enumerate(index_to_name) if not is_virtual_gin_node(n)]


def indices_by_entity_layer(
    index_to_name: Sequence[str],
    *,
    layer_map: Mapping[str, str] | None = None,
    include_virtual: bool = False,
) -> dict[str, list[int]]:
    """Bucket node indices by ``assign_entity_layer`` (default skips virtual)."""
    buckets: dict[str, list[int]] = {}
    for i, name in enumerate(index_to_name):
        layer = assign_entity_layer(name, layer_map)
        if layer == "virtual" and not include_virtual:
            continue
        buckets.setdefault(layer, []).append(i)
    return buckets


def _looks_like_qualified_gin_name(uid: str) -> bool:
    """Heuristic: user supplied a GIN-style name that failed membership."""
    if is_virtual_gin_node(uid):
        return True
    lower = uid.lower()
    return lower.endswith(".rna") or lower.endswith(".protein") or ".promoter." in lower


@dataclass
class MappingQCReport:
    """Structured QC after mapping user IDs to GIN nodes."""

    gin_version: str
    input_count: int
    mapped_count: int
    unmapped_ids: list[str] = field(default_factory=list)
    unmapped_detail: list[dict[str, str]] = field(default_factory=list)
    ambiguous_mappings: list[dict[str, object]] = field(default_factory=list)
    layer_counts: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(
            {
                "gin_version": self.gin_version,
                "input_count": self.input_count,
                "mapped_count": self.mapped_count,
                "mapping_rate": self.mapped_count / max(self.input_count, 1),
                "unmapped_ids": self.unmapped_ids[:500],
                "unmapped_truncated": max(0, len(self.unmapped_ids) - 500),
                "unmapped_detail": self.unmapped_detail[:500],
                "unmapped_detail_truncated": max(0, len(self.unmapped_detail) - 500),
                "ambiguous_mappings": self.ambiguous_mappings[:100],
                "layer_counts": self.layer_counts,
                "warnings": self.warnings,
            },
            indent=2,
        )


def default_symbol_to_gin_candidates(
    symbol: str,
    *,
    prefer_layers: Sequence[str] = ("RNA", "protein", "rna"),
) -> list[str]:
    """
    Generate candidate GIN node names from a HGNC-like gene symbol.

    GINv3.11 uses ``GENE.RNA`` (uppercase) and ``GENE.protein``. Legacy graphs
    may use ``GENE.rna``; callers verify existence against the graph index.
    """
    sym = symbol.strip()
    if not sym:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for layer in prefer_layers:
        cand: str | None = None
        if layer in ("rna", "RNA"):
            cand = f"{sym}.RNA" if layer == "RNA" else f"{sym}.rna"
        elif layer == "protein":
            cand = f"{sym}.protein"
        if cand and cand not in seen:
            seen.add(cand)
            out.append(cand)
    return out


def load_symbol_blacklist(path: str | Path | None) -> set[str]:
    """Load one HGNC symbol per line; ``#`` comments and blanks are ignored."""
    if path is None:
        return set()
    p = Path(path)
    if not p.is_file():
        return set()
    return {
        ln.strip()
        for ln in p.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    }


def load_symbol_aliases(path: str | Path | None) -> dict[str, str]:
    """
    Load ``old_symbol\\tapproved_symbol`` (header optional).

    Keys and values are stripped; duplicate old symbols: last wins.
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.is_file():
        return {}
    out: dict[str, str] = {}
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        old_sym, new_sym = parts[0].strip(), parts[1].strip()
        if old_sym and new_sym:
            out[old_sym] = new_sym
    return out


def apply_symbol_aliases(
    user_ids: Iterable[str],
    aliases: Mapping[str, str],
) -> tuple[list[str], dict[str, str]]:
    """
    Expand user IDs through an alias table.

    Returns (resolved_id_list, alias_applied) where ``alias_applied`` maps each
    *input* id to the symbol used for mapping (unchanged when no alias).
    """
    uid_list = [u.strip() for u in user_ids if str(u).strip()]
    if not aliases:
        return uid_list, {u: u for u in uid_list}
    alias_applied: dict[str, str] = {}
    resolved: list[str] = []
    for uid in uid_list:
        mapped = aliases.get(uid, uid)
        alias_applied[uid] = mapped
        resolved.append(mapped)
    return resolved, alias_applied


def map_ids_to_gin_nodes(
    user_ids: Iterable[str],
    gin_node_is_member: Callable[[str], bool],
    *,
    gin_version: str = "GINv3.11",
    candidate_fn: Callable[[str], list[str]] | None = None,
    missing_strategy: str = "warn",
    blocked_ids: set[str] | None = None,
    symbol_aliases: Mapping[str, str] | None = None,
) -> tuple[dict[str, str], MappingQCReport]:
    """
    Map each user ID to one GIN node name if any candidate exists in the graph.

    Parameters
    ----------
    user_ids
        Gene symbols, UniProt accessions, or already-qualified GIN names.
    gin_node_is_member
        O(1) or cheap membership test (e.g. ``set_of_nodes.__contains__``).
    candidate_fn
        Maps one user ID to ordered GIN name candidates; first hit wins.
    blocked_ids
        HGNC symbols rejected before graph lookup (``reason=blacklisted``).
    symbol_aliases
        Optional old_symbol → approved_symbol map applied before candidate expansion.
    missing_strategy
        ``warn`` (append to report.warnings), ``silent``, or ``strict`` (raise
        KeyError on first unmapped id).
    """
    candidate_fn = candidate_fn or default_symbol_to_gin_candidates
    blocked = blocked_ids or set()
    aliases = symbol_aliases or {}
    uid_list = [u.strip() for u in user_ids if str(u).strip()]
    mapping: dict[str, str] = {}
    unmapped: list[str] = []
    unmapped_detail: list[dict[str, str]] = []
    ambiguous_mappings: list[dict[str, object]] = []
    for uid in uid_list:
        if uid in blocked:
            unmapped.append(uid)
            unmapped_detail.append({"id": uid, "reason": "blacklisted"})
            if missing_strategy == "strict":
                raise KeyError(f"Blacklisted user id: {uid!r}")
            continue
        lookup = aliases.get(uid, uid)
        if gin_node_is_member(lookup):
            mapping[uid] = lookup
            continue
        if lookup != uid and gin_node_is_member(uid):
            mapping[uid] = uid
            continue
        hits = [c for c in candidate_fn(lookup) if gin_node_is_member(c)]
        if len(hits) > 1:
            mapping[uid] = hits[0]
            ambiguous_mappings.append(
                {
                    "id": uid,
                    "reason": "ambiguous",
                    "candidates": hits,
                    "chosen": hits[0],
                }
            )
        elif len(hits) == 1:
            mapping[uid] = hits[0]
        else:
            unmapped.append(uid)
            reason = (
                "node_not_in_graph"
                if _looks_like_qualified_gin_name(uid)
                else "no_candidate"
            )
            unmapped_detail.append({"id": uid, "reason": reason})
            if missing_strategy == "strict":
                raise KeyError(f"No GIN node for user id: {uid!r}")

    layer_counts: dict[str, int] = {}
    for gnode in mapping.values():
        layer = gin_entity_layer(gnode)
        layer_counts[layer] = layer_counts.get(layer, 0) + 1

    warnings: list[str] = []
    if unmapped and missing_strategy == "warn":
        warnings.append(
            f"{len(unmapped)} ids could not be mapped to any GIN node "
            f"(showing up to 20): {unmapped[:20]!r}"
        )

    if ambiguous_mappings and missing_strategy == "warn":
        warnings.append(
            f"{len(ambiguous_mappings)} ids had multiple GIN candidates; first hit was used."
        )

    report = MappingQCReport(
        gin_version=gin_version,
        input_count=len(uid_list),
        mapped_count=len(mapping),
        unmapped_ids=unmapped,
        unmapped_detail=unmapped_detail,
        ambiguous_mappings=ambiguous_mappings,
        layer_counts=layer_counts,
        warnings=warnings,
    )
    return mapping, report


def build_metabolite_regex_mapper(
    pattern_to_gin: Mapping[str, str],
) -> Callable[[str], list[str]]:
    """Example hook: map HMDB / InChIKey-like strings via regex to GIN names."""

    compiled: list[tuple[re.Pattern[str], str]] = [
        (re.compile(p), g) for p, g in pattern_to_gin.items()
    ]

    def _fn(uid: str) -> list[str]:
        for rx, gnode in compiled:
            if rx.search(uid):
                return [gnode]
        return []

    return _fn


def merge_candidate_fns(
    *fns: Callable[[str], list[str]],
) -> Callable[[str], list[str]]:
    """Try several candidate generators in order."""

    def _merged(uid: str) -> list[str]:
        seen: set[str] = set()
        out: list[str] = []
        for fn in fns:
            for c in fn(uid):
                if c not in seen:
                    seen.add(c)
                    out.append(c)
        return out

    return _merged


def iter_unique_gin_nodes_from_sifi(path: str, encoding: str = "utf-8") -> Iterator[str]:
    """Stream unique node names from a .sifi / .sifi.gz (two columns, tab-separated)."""
    from gin_tool.io import open_sifi_text

    seen: set[str] = set()
    with open_sifi_text(path, encoding=encoding) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            for n in (a, b):
                if n not in seen:
                    seen.add(n)
                    yield n

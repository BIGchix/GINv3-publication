"""Load per-layer GIN node pools from ``nodeList/`` files for diagnostics and mapping."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Mapping, Sequence

# Filename in nodeList/ -> diagnostic layer id (seed pool + closure tag).
NODE_LIST_FILE_TO_LAYER: dict[str, str] = {
    "nodeList_mRNA.txt": "rna",
    "nodeList_miRNA.txt": "mirna",
    "nodeList_lncRNA.txt": "lncrna",
    "nodeList_protein.txt": "protein",
    "nodeList_TF.txt": "tf",
    "nodeList_promoter.txt": "promoter_chromatin",
    "nodeList_metabolite.txt": "metabolite",
    "nodeList_His.txt": "histone_mark",
    "nodeList_HME.txt": "hme",
    "nodeList_complex.txt": "complex",
}

# When a GIN node appears in multiple lists, the earliest layer here wins for closure tags.
LAYER_PRIORITY: tuple[str, ...] = (
    "complex",
    "metabolite",
    "promoter_chromatin",
    "histone_mark",
    "hme",
    "tf",
    "protein",
    "mirna",
    "lncrna",
    "rna",
)

DEFAULT_DIAGNOSTIC_LAYERS: tuple[str, ...] = (
    "rna",
    "protein",
    "tf",
    "promoter_chromatin",
    "metabolite",
)

# Scheme D: protein/TF use runtime null at k=1 only; never precompute null grids.
RUNTIME_ONLY_LAYERS: frozenset[str] = frozenset({"protein", "tf"})

# Max k-hop allowed per layer (diagnosis, analyze, runtime null).
MAX_K_BY_LAYER: dict[str, int] = {
    "protein": 1,
    "tf": 1,
    "hme": 1,  # HME seeds often yield very large k=2 balls on GINv3.10
}

# Layers included in ``precompute_null`` when ``--node-list-dir`` is set (default).
PRECOMPUTE_DEFAULT_LAYERS: tuple[str, ...] = (
    "rna",
    "mirna",
    "lncrna",
    "promoter_chromatin",
    "metabolite",
    "complex",
)

UNKNOWN_SEED_LAYER = "unlisted"

# All layer ids used for stratified seed sampling (nodeList + unlisted).
ALL_STRATIFIED_LAYERS: tuple[str, ...] = tuple(
    dict.fromkeys(list(NODE_LIST_FILE_TO_LAYER.values()) + [UNKNOWN_SEED_LAYER])
)

# Stratified test draws: each nodeList layer 10, tf 5, unlisted 10.
DEFAULT_STRATIFIED_LAYER_SEED_COUNTS: dict[str, int] = {
    layer: (5 if layer == "tf" else 10) for layer in ALL_STRATIFIED_LAYERS
}

ZERO_LAYER_SEED_COUNTS: dict[str, int] = {layer: 0 for layer in ALL_STRATIFIED_LAYERS}


def parse_layer_seed_counts(
    spec: str | None,
    *,
    defaults: Mapping[str, int] | None = None,
    layer_only: bool = False,
) -> dict[str, int]:
    """
    Parse ``layer:n,layer:n,...`` into per-layer sample sizes.

    Parameters
    ----------
    layer_only
        If True, start from all layers at 0 and **only** set counts listed in
        ``spec`` (non-listed layers stay 0). Requires a non-empty ``spec``.
        If False, merge ``spec`` onto ``defaults`` (default stratified table).
    """
    if layer_only:
        if not spec or not str(spec).strip():
            raise ValueError(
                "--layer-only requires --layer-seeds with at least one layer:count "
                "(e.g. --layer-seeds protein:10,tf:5)"
            )
        base = dict(ZERO_LAYER_SEED_COUNTS)
    else:
        base = dict(defaults or DEFAULT_STRATIFIED_LAYER_SEED_COUNTS)
        if not spec or not str(spec).strip():
            return base
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(
                f"Invalid layer:count segment {part!r}; use e.g. tf:5,protein:10"
            )
        layer, _, n_s = part.partition(":")
        layer = layer.strip()
        n = int(n_s.strip())
        if n < 0:
            raise ValueError(f"Layer count must be >= 0, got {n} for {layer!r}")
        base[layer] = n
    return base


def sample_seeds_by_layer_counts(
    rng: random.Random,
    index_to_name: Sequence[str],
    buckets: Mapping[str, Sequence[int]],
    unlisted_pool: Sequence[int],
    layer_counts: Mapping[str, int],
) -> tuple[list[str], dict[str, int], list[str]]:
    """
    Sample ``layer_counts[layer]`` node indices per layer; shuffle combined names.

    ``unlisted`` uses ``unlisted_pool``; other layers use ``buckets[layer]``.
    """
    names: list[str] = []
    actual: dict[str, int] = {}
    warnings: list[str] = []

    for layer, n_want in layer_counts.items():
        n_want = int(n_want)
        if n_want <= 0:
            actual[layer] = 0
            continue
        if layer == UNKNOWN_SEED_LAYER:
            pool = list(unlisted_pool)
        else:
            pool = list(buckets.get(layer, []))
        if not pool:
            warnings.append(f"layer {layer}: empty pool (wanted {n_want})")
            actual[layer] = 0
            continue
        n_pick = min(n_want, len(pool))
        if n_pick < n_want:
            warnings.append(
                f"layer {layer}: pool size {len(pool)} < requested {n_want}; picked {n_pick}"
            )
        picks = rng.sample(pool, n_pick)
        names.extend(index_to_name[i] for i in picks)
        actual[layer] = n_pick

    rng.shuffle(names)
    return names, actual, warnings


def max_k_for_layer(layer: str, *, default: int = 2) -> int:
    """Upper bound on k-hop for this entity layer."""
    return MAX_K_BY_LAYER.get(layer, default)


def default_k_for_layer(layer: str | None) -> int:
    """Default k-hop for analysis and runtime null: protein/tf/hme → 1, others → 2."""
    if layer is None:
        return 2
    return max_k_for_layer(layer)


def layer_for_seed(name: str, layer_map: Mapping[str, str]) -> str | None:
    """Layer tag from nodeList, or ``None`` if the node is not in any list."""
    return layer_map.get(name)


def k_for_seed(name: str, layer_map: Mapping[str, str], *, user_k: int | None = None) -> int:
    """k-hop depth for one seed (per-layer default, or global ``user_k`` override)."""
    if user_k is not None:
        return int(user_k)
    ly = layer_for_seed(name, layer_map)
    return default_k_for_layer(ly)


def summarize_seed_layers(
    gin_node_names: Sequence[str],
    layer_map: Mapping[str, str],
    *,
    user_k: int | None = None,
) -> tuple[dict[str, int], dict[str, str | None], dict[str, int]]:
    """
    Per-seed layer tags and k values.

    Returns
    -------
    layer_counts
        Count of seeds per layer tag (unlisted seeds grouped under ``unlisted``).
    layer_by_name
        ``gin_node -> layer`` (``None`` in map stored as ``unlisted`` in counts only).
    k_by_name
        ``gin_node -> k`` for subgraph expansion.
    """
    layer_counts: dict[str, int] = {}
    layer_by_name: dict[str, str | None] = {}
    k_by_name: dict[str, int] = {}

    for name in gin_node_names:
        ly = layer_for_seed(name, layer_map)
        layer_by_name[name] = ly
        bucket = ly if ly is not None else UNKNOWN_SEED_LAYER
        layer_counts[bucket] = layer_counts.get(bucket, 0) + 1
        k_by_name[name] = k_for_seed(name, layer_map, user_k=user_k)

    return layer_counts, layer_by_name, k_by_name


def resolve_k_and_seed_layer(
    gin_node_names: Sequence[str],
    layer_map: Mapping[str, str],
    *,
    user_k: int | None = None,
) -> tuple[int | None, str | None, list[str]]:
    """
    Legacy helper: single ``k`` when all seeds share one k; else ``k=None`` (use per-seed).

    Mixed protein + RNA no longer raises; callers should use ``summarize_seed_layers``.
    """
    layer_counts, _, k_by_name = summarize_seed_layers(
        gin_node_names, layer_map, user_k=user_k
    )
    notes: list[str] = []
    unique_ks = set(k_by_name.values())
    tagged_layers = {ly: c for ly, c in layer_counts.items() if ly != UNKNOWN_SEED_LAYER}

    if not tagged_layers and UNKNOWN_SEED_LAYER in layer_counts:
        k = int(user_k) if user_k is not None else 2
        notes.append("No seeds found in nodeList; using k=2.")
        return k, None, notes

    if len(unique_ks) > 1:
        detail = ", ".join(
            f"{ly}→k={default_k_for_layer(ly if ly != UNKNOWN_SEED_LAYER else None)}"
            for ly in sorted(tagged_layers)
        )
        notes.append(f"Mixed layers ({detail}); per-seed k-hop union subgraph.")
        return None, None, notes

    k_val = next(iter(unique_ks))
    seed_layer = max(tagged_layers.items(), key=lambda x: (x[1], x[0]))[0] if tagged_layers else None
    if user_k is None and seed_layer:
        notes.append(f"Uniform k={k_val} for layer {seed_layer!r}.")
    elif user_k is not None:
        notes.append(f"User override k={k_val} for all seeds.")
    return k_val, seed_layer, notes


def filter_k_values_for_layer(layer: str, k_values: Sequence[int]) -> list[int]:
    """Keep only k values allowed for ``layer``."""
    cap = max_k_for_layer(layer)
    return sorted({int(k) for k in k_values if 1 <= int(k) <= cap})


def infer_dominant_seed_layer(
    gin_node_names: Sequence[str],
    layer_map: Mapping[str, str],
) -> str | None:
    """
    Best-effort layer tag for a user seed list using nodeList-derived ``layer_map``.

    Returns the most frequent layer among mapped seeds, or None if none are listed.
    """
    counts: dict[str, int] = {}
    for name in gin_node_names:
        ly = layer_map.get(name)
        if ly:
            counts[ly] = counts.get(ly, 0) + 1
    if not counts:
        return None
    return max(counts.items(), key=lambda x: (x[1], x[0]))[0]


def validate_k_for_seed_layer(layer: str | None, k: int) -> str | None:
    """If ``k`` exceeds the cap for ``layer``, return an error message; else None."""
    if layer is None:
        return None
    cap = max_k_for_layer(layer)
    if int(k) > cap:
        return (
            f"layer {layer!r} allows k <= {cap} (got k={k}). "
            f"protein, tf, and hme are capped at k=1."
        )
    return None


def unlisted_entity_indices(
    index_to_name: Sequence[str],
    layer_map: Mapping[str, str],
) -> list[int]:
    """Real GIN nodes not present in any nodeList-derived ``layer_map`` entry."""
    from gin_tool.mapping import is_virtual_gin_node

    return [
        i
        for i, name in enumerate(index_to_name)
        if not is_virtual_gin_node(name) and name not in layer_map
    ]


def read_node_list_file(path: Path) -> list[str]:
    """One GIN node id or alias per line; blanks and ``#`` comments skipped."""
    ids: list[str] = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.append(line)
    return ids


def load_node_list_dir(node_list_dir: Path) -> dict[str, list[str]]:
    """
    Read all standard ``nodeList_*.txt`` files under ``node_list_dir``.

    Returns ``layer_id -> list of raw ids`` (file order preserved, duplicates kept).
    """
    out: dict[str, list[str]] = {}
    if not node_list_dir.is_dir():
        return out
    for fname, layer in NODE_LIST_FILE_TO_LAYER.items():
        p = node_list_dir / fname
        if p.is_file():
            out[layer] = read_node_list_file(p)
    return out


def resolve_id_to_gin_node(uid: str, name_to_index: Mapping[str, int]) -> str | None:
    """Map one list entry to a GIN node name present in the graph index."""
    uid = uid.strip()
    if not uid:
        return None
    if uid in name_to_index:
        return uid
    return None


def resolve_layer_pools(
    layer_to_ids: Mapping[str, Sequence[str]],
    name_to_index: Mapping[str, int],
) -> tuple[dict[str, list[int]], dict[str, list[str]], dict[str, int]]:
    """
    Resolve list entries to graph indices per layer.

    Returns
    -------
    buckets
        ``layer -> sorted unique node indices`` (eligible seed pools).
    unmapped_by_layer
        ``layer -> ids not found in ``name_to_index``.
    mapped_counts
        ``layer -> number of unique resolved nodes``.
    """
    buckets: dict[str, list[int]] = {}
    unmapped: dict[str, list[str]] = {}
    mapped_counts: dict[str, int] = {}

    for layer, ids in layer_to_ids.items():
        seen_idx: set[int] = set()
        missing: list[str] = []
        for uid in ids:
            gnode = resolve_id_to_gin_node(uid, name_to_index)
            if gnode is None:
                missing.append(uid)
                continue
            seen_idx.add(name_to_index[gnode])
        buckets[layer] = sorted(seen_idx)
        if missing:
            unmapped[layer] = missing
        mapped_counts[layer] = len(seen_idx)

    return buckets, unmapped, mapped_counts


def build_layer_map_from_pools(
    layer_to_gin_nodes: Mapping[str, Sequence[str]],
) -> dict[str, str]:
    """
    ``gin_node -> layer`` for closure composition.

    Nodes listed in multiple files get the layer with highest priority in
    ``LAYER_PRIORITY``.
    """
    layer_map: dict[str, str] = {}
    order = {ly: i for i, ly in enumerate(LAYER_PRIORITY)}
    default_rank = len(LAYER_PRIORITY)

    for layer, nodes in layer_to_gin_nodes.items():
        rank = order.get(layer, default_rank)
        for node in nodes:
            node = node.strip()
            if not node:
                continue
            prev = layer_map.get(node)
            if prev is None or rank < order.get(prev, default_rank):
                layer_map[node] = layer

    return layer_map


def build_pools_and_layer_map(
    node_list_dir: Path,
    index_to_name: Sequence[str],
    name_to_index: Mapping[str, int],
) -> tuple[
    dict[str, list[int]],
    dict[str, str],
    dict[str, list[str]],
    dict[str, int],
    dict[str, str],
]:
    """
    Full pipeline: load nodeList dir → resolve indices + closure layer map.

    Returns
    -------
    buckets, layer_map, unmapped_by_layer, mapped_counts, source_files
    """
    from gin_tool.mapping import is_virtual_gin_node

    layer_to_ids = load_node_list_dir(node_list_dir)
    buckets, unmapped, mapped_counts = resolve_layer_pools(layer_to_ids, name_to_index)

    # Exclude virtual nodes from seed pools.
    for layer, indices in list(buckets.items()):
        buckets[layer] = [
            i for i in indices if not is_virtual_gin_node(index_to_name[i])
        ]
        mapped_counts[layer] = len(buckets[layer])

    layer_to_gin: dict[str, list[str]] = {}
    for layer, indices in buckets.items():
        layer_to_gin[layer] = [index_to_name[i] for i in indices]

    layer_map = build_layer_map_from_pools(layer_to_gin)

    source_files = {
        layer: str((node_list_dir / fname).resolve())
        for fname, layer in NODE_LIST_FILE_TO_LAYER.items()
        if layer in layer_to_ids
    }

    return buckets, layer_map, unmapped, mapped_counts, source_files

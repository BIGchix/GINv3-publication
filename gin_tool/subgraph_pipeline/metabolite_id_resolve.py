"""Resolve metabolite synonyms to GIN node names (local index + optional PubChem).

GINv2.0 UC_ID convention (npj Syst Biol Appl 2024, BIGchix/GINv2.0):
external chemical ids are aggregated into unified types
``CID`` / ``SID`` / ``CAS`` / ``KEGG`` / ``HMDB`` / ``ChEBI``, stored in GIN as
``CID…``, ``SID…``, ``CAS…``, ``KEGGC…``, ``HMDB…``, ``CHEBI…`` (and composites
joined by ``_``; reaction intermediates may also embed UC_IDs separated by ``;``).
"""

from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_CAS_RE = re.compile(r"^CAS(?P<num>\d{2,7}-\d{2}-\d)$", re.I)
_CID_NODE_RE = re.compile(r"^CID(?P<cid>\d+)$", re.I)
_CHEBI_NODE_RE = re.compile(r"^CHEBI(?P<id>\d+)$", re.I)
_SID_NODE_RE = re.compile(r"^SID(?P<sid>\d+)$", re.I)
_KEGGC_NODE_RE = re.compile(r"^KEGGC(?P<id>\d+)$", re.I)
_HMDB_NODE_RE = re.compile(r"^HMDB(?P<id>\d+)$", re.I)
_KEGG_CPD_RE = re.compile(r"^(?:cpd:|KEGG:)?C(?P<id>\d{5})$", re.I)
_UC_TOKEN_RE = re.compile(
    r"^(CAS\d{2,7}-\d{2}-\d|CHEBI\d+|SID\d+|CID\d+|KEGGC\d+|HMDB\d+)$",
    re.I,
)
_UC_PREFIXES = ("CAS", "CHEBI", "SID", "CID", "KEGGC", "HMDB")

# GINv2 UC_ID source preference (paper Suppl. Fig. 1 / Methods):
# aggregate external IDs into one canonical type, preferring CID first.
UC_ID_PRIORITY: tuple[str, ...] = ("CID", "SID", "CAS", "KEGGC", "HMDB", "CHEBI")

_TABLE1_ALIASES_PATH = Path(__file__).resolve().parent / "data" / "ginv2_table1_uc_aliases.json"
_table1_aliases_cache: dict[str, dict] | None = None


def load_ginv2_table1_aliases(path: Path | None = None) -> dict[str, dict]:
    """Curated KEGG→UC_ID extras from GINv2 paper Table 1 (ATP/ADP hubs, etc.)."""
    global _table1_aliases_cache
    if path is None and _table1_aliases_cache is not None:
        return _table1_aliases_cache
    p = path or _TABLE1_ALIASES_PATH
    if not p.is_file():
        return {}
    raw = json.loads(p.read_text(encoding="utf-8"))
    aliases = {str(k): v for k, v in (raw.get("aliases") or {}).items()}
    if path is None:
        _table1_aliases_cache = aliases
    return aliases


def normalize_kegg_compound_id(raw: str) -> str | None:
    """Return canonical ``C#####`` from ``C00031`` / ``cpd:C00031`` / ``KEGGC00031``."""
    s = raw.strip()
    if not s:
        return None
    m = _KEGGC_NODE_RE.match(s)
    if m:
        return f"C{int(m.group('id')):05d}"
    m = _KEGG_CPD_RE.match(s)
    if m:
        return f"C{int(m.group('id')):05d}"
    return None


def _normalize_cas(raw: str) -> str | None:
    s = raw.strip().upper().replace(" ", "")
    m = _CAS_RE.match(s)
    if m:
        return m.group("num")
    if re.match(r"^\d{2,7}-\d{2}-\d$", s):
        return s
    return None


def iter_uc_id_tokens(node: str) -> list[str]:
    """Extract UC_ID tokens from a real or composite/virtual GIN node name."""
    out: list[str] = []
    for part in re.split(r"[;_]", node):
        p = part.strip()
        if _UC_TOKEN_RE.match(p):
            out.append(p)
    return out


def _alias_keys_for_gin_node(node: str) -> set[str]:
    keys: set[str] = {node, node.upper(), node.lower()}
    m = _CID_NODE_RE.match(node)
    if m:
        cid = m.group("cid")
        keys.update({f"CID{cid}", f"cid{cid}", cid, f"PUBCHEM.COMPOUND:{cid}"})
    m = _CHEBI_NODE_RE.match(node)
    if m:
        cid = m.group("id")
        keys.update({f"CHEBI:{cid}", f"chebi:{cid}", f"CHEBI{cid}", cid})
    m = _SID_NODE_RE.match(node)
    if m:
        sid = m.group("sid")
        keys.update({f"SID{sid}", f"pubchem.sid:{sid}", sid})
    m = _KEGGC_NODE_RE.match(node)
    if m:
        kid = f"{int(m.group('id')):05d}"
        keys.update(
            {
                f"KEGGC{kid}",
                f"C{kid}",
                f"cpd:C{kid}",
                f"CPD:C{kid}",
                f"KEGG:C{kid}",
                f"kegg:C{kid}",
            }
        )
    m = _HMDB_NODE_RE.match(node)
    if m:
        hid = m.group("id")
        keys.update({f"HMDB{hid}", f"hmdb{hid}", f"HMDB:{hid}"})
    cas = _normalize_cas(node)
    if cas:
        keys.update({cas, f"CAS{cas}", f"CAS-{cas}"})
    # Composite / virtual embeddings: CAS50-99-7_CHEBI4167 or …;SID3333;…
    for tok in iter_uc_id_tokens(node):
        if tok != node:
            keys.update(_alias_keys_for_gin_node(tok))
            keys.add(tok)
    return keys


def collect_chem_nodes_from_sifi(sifi_path: str | Path) -> set[str]:
    """Real (non-virtual) chemical UC_ID nodes present in a ``.sifi`` graph.

    Matches GINv2 UC_ID tokens (``SID3333``, ``CAS50-99-7``, composites
    ``CAS…_…``). Does **not** treat gene symbols like ``CASP3`` as chemicals.
    """
    path = Path(sifi_path)
    out: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            for n in parts[:2]:
                if ";" in n:
                    continue
                if _UC_TOKEN_RE.match(n) or (
                    "_" in n and any(_UC_TOKEN_RE.match(p) for p in n.split("_")[:1])
                ):
                    # full node if first segment is a UC_ID (composite CAS_x / SID_x)
                    if _UC_TOKEN_RE.match(n.split("_")[0]):
                        out.add(n)
    return out


@dataclass
class GinMetaboliteIndex:
    """Local alias index over GIN metabolite node names (GINv2 UC_IDs)."""

    nodes: set[str] = field(default_factory=set)
    _alias_to_node: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_node_list(cls, path: str | Path) -> GinMetaboliteIndex:
        idx = cls()
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            node = line.strip()
            if node:
                idx.register(node)
        return idx

    @classmethod
    def from_pool(cls, pool: set[str]) -> GinMetaboliteIndex:
        idx = cls()
        for node in pool:
            idx.register(node)
        return idx

    @classmethod
    def from_sifi(cls, sifi_path: str | Path) -> GinMetaboliteIndex:
        return cls.from_pool(collect_chem_nodes_from_sifi(sifi_path))

    def register(self, node: str, *, prefer: bool = False) -> None:
        """Register a node. Prefer real UC_ID nodes over virtual/composite aliases."""
        is_virtual = ";" in node
        if not is_virtual:
            self.nodes.add(node)
        for key in _alias_keys_for_gin_node(node):
            if prefer or key not in self._alias_to_node:
                # Do not let virtual nodes overwrite a real mapping.
                existing = self._alias_to_node.get(key)
                if existing is not None and ";" not in existing and is_virtual:
                    continue
                if prefer or existing is None or (";" in existing and not is_virtual):
                    self._alias_to_node[key] = node
                    self._alias_to_node[key.upper()] = node
                    self._alias_to_node[key.lower()] = node

    def lookup(self, query: str) -> str | None:
        q = query.strip()
        if not q:
            return None
        if q in self.nodes:
            return q
        for cand in (q, q.upper(), q.lower()):
            hit = self._alias_to_node.get(cand)
            if hit:
                return hit
        cas = _normalize_cas(q)
        if cas:
            for cand in (cas, f"CAS{cas}"):
                hit = self._alias_to_node.get(cand) or self._alias_to_node.get(cand.upper())
                if hit:
                    return hit
        kegg = normalize_kegg_compound_id(q)
        if kegg:
            for cand in (f"KEGGC{kegg[1:]}", kegg, f"cpd:{kegg}"):
                hit = self._alias_to_node.get(cand) or self._alias_to_node.get(cand.upper())
                if hit:
                    return hit
        return None


@dataclass
class MetaboliteResolveResult:
    query: str
    gin_node: str | None
    method: str
    pubchem_cid: int | None = None
    matched_alias: str | None = None
    warnings: list[str] = field(default_factory=list)


class PubChemMetaboliteResolver:
    """Online synonym/CAS/name → PubChem CID (NCBI PUG REST, no API key)."""

    def __init__(
        self,
        *,
        cache_path: Path | None = None,
        min_interval_s: float = 0.25,
        timeout_s: float = 20.0,
    ) -> None:
        self.cache_path = cache_path
        self.min_interval_s = min_interval_s
        self.timeout_s = timeout_s
        self._cache: dict[str, dict] = {}
        self._last_req = 0.0
        if cache_path and cache_path.is_file():
            try:
                self._cache = json.loads(cache_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                self._cache = {}

    def _save_cache(self) -> None:
        if not self.cache_path:
            return
        try:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            self.cache_path.write_text(json.dumps(self._cache, indent=2), encoding="utf-8")
        except OSError:
            return

    def _throttle(self) -> None:
        wait = self.min_interval_s - (time.monotonic() - self._last_req)
        if wait > 0:
            time.sleep(wait)
        self._last_req = time.monotonic()

    def _get_json(self, url: str) -> dict | None:
        if url in self._cache:
            return self._cache[url]
        self._throttle()
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout_s) as resp:
                data = json.loads(resp.read().decode("utf-8"))
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            self._cache[url] = None  # type: ignore[assignment]
            self._save_cache()
            return None
        self._cache[url] = data
        self._save_cache()
        return data

    def cids_for_name(self, name: str) -> list[int]:
        enc = urllib.parse.quote(name.strip())
        url = f"{_PUBCHEM_BASE}/compound/name/{enc}/cids/JSON"
        data = self._get_json(url)
        if not data:
            return []
        return list(data.get("IdentifierList", {}).get("CID", []) or [])

    def cids_for_sid(self, sid: int | str) -> list[int]:
        sid_s = str(sid).removeprefix("SID").removeprefix("sid")
        url = f"{_PUBCHEM_BASE}/substance/sid/{sid_s}/cids/JSON"
        data = self._get_json(url)
        if not data:
            return []
        info = data.get("InformationList", {}).get("Information", [])
        if not info:
            return []
        return list(info[0].get("CID", []) or [])

    def synonyms_for_cid(self, cid: int, *, limit: int = 50) -> list[str]:
        url = f"{_PUBCHEM_BASE}/compound/cid/{int(cid)}/synonyms/JSON"
        data = self._get_json(url)
        if not data:
            return []
        info = data.get("InformationList", {}).get("Information", [])
        if not info:
            return []
        syns = info[0].get("Synonym", []) or []
        return [str(s) for s in syns[:limit]]


@dataclass
class KeggCompoundXrefs:
    """External IDs for one KEGG compound, ready for GINv2-priority alignment."""

    cpd: str
    cids: list[str] = field(default_factory=list)  # CID#####
    sids: list[str] = field(default_factory=list)  # SID#####
    cases: list[str] = field(default_factory=list)  # CAS#-#-#
    keggc: list[str] = field(default_factory=list)  # KEGGC#####
    hmdbs: list[str] = field(default_factory=list)  # HMDB#####
    chebis: list[str] = field(default_factory=list)  # CHEBI#####

    def candidates_by_priority(self) -> list[tuple[str, str]]:
        """Return ``(uc_type, token)`` in GINv2 preference order."""
        buckets: dict[str, list[str]] = {
            "CID": self.cids,
            "SID": self.sids,
            "CAS": self.cases,
            "KEGGC": self.keggc,
            "HMDB": self.hmdbs,
            "CHEBI": self.chebis,
        }
        out: list[tuple[str, str]] = []
        seen: set[str] = set()
        for kind in UC_ID_PRIORITY:
            for tok in buckets.get(kind, []):
                if tok not in seen:
                    seen.add(tok)
                    out.append((kind, tok))
        return out


def align_kegg_xrefs_to_gin(
    xrefs: KeggCompoundXrefs,
    index: GinMetaboliteIndex,
    *,
    allow_virtual: bool = False,
) -> tuple[str | None, str | None]:
    """Map a KEGG compound to a GIN node using UC_ID priority (CID first).

    Returns ``(gin_node, matched_uc_type)`` or ``(None, None)``.
    """
    for kind, tok in xrefs.candidates_by_priority():
        hit = index.lookup(tok)
        if not hit:
            continue
        if ";" in hit and not allow_virtual:
            continue
        return hit, kind
    if allow_virtual:
        for kind, tok in xrefs.candidates_by_priority():
            hit = index.lookup(tok)
            if hit:
                return hit, kind
    return None, None


def build_kegg_compound_xrefs(
    cpd: str,
    *,
    chebi_map: dict[str, str] | None = None,
    sid_map: dict[str, str] | None = None,
    sid2cid: dict[str, list[int]] | None = None,
    cas_map: dict[str, list[str]] | None = None,
    hmdb_map: dict[str, list[str]] | None = None,
    curated_aliases: dict[str, dict] | None = None,
    include_table1_aliases: bool = True,
) -> KeggCompoundXrefs:
    """Assemble UC_ID candidates for ``C#####`` from KEGG + PubChem + curated hubs.

    Curated aliases (GINv2 paper Table 1) add charge-state SIDs such as
    ``SID8148096`` (ATP4-) / ``SID85646635`` (ADP3-), which exist in GINv2/v3
    even when KEGG's ``SID3304`` / ``SID3310`` do not.
    """
    canon = normalize_kegg_compound_id(cpd) or cpd.strip().upper()
    x = KeggCompoundXrefs(
        cpd=canon,
        keggc=[f"KEGGC{canon[1:]}"] if canon.startswith("C") else [],
    )
    if chebi_map and canon in chebi_map:
        x.chebis.append(chebi_map[canon])
    if sid_map and canon in sid_map:
        sid_tok = sid_map[canon]
        x.sids.append(sid_tok)
        sid_num = sid_tok.removeprefix("SID")
        for cid in (sid2cid or {}).get(sid_num, []) or (sid2cid or {}).get(sid_tok, []):
            x.cids.append(f"CID{int(cid)}")
    for cas in (cas_map or {}).get(canon, []):
        cas_n = _normalize_cas(cas) or _normalize_cas(f"CAS{cas}")
        if cas_n:
            x.cases.append(f"CAS{cas_n}")
    for hmdb in (hmdb_map or {}).get(canon, []):
        h = hmdb.strip().upper().replace("HMDB:", "")
        if h.startswith("HMDB"):
            x.hmdbs.append(h)
        elif h.isdigit():
            x.hmdbs.append(f"HMDB{h}")

    aliases = curated_aliases
    if aliases is None and include_table1_aliases:
        aliases = load_ginv2_table1_aliases()
    extra = (aliases or {}).get(canon) or {}
    for cid in extra.get("CID") or []:
        x.cids.append(cid if str(cid).startswith("CID") else f"CID{cid}")
    for sid in extra.get("SID") or []:
        x.sids.append(sid if str(sid).startswith("SID") else f"SID{sid}")
    for cas in extra.get("CAS") or []:
        cas_n = _normalize_cas(str(cas)) or _normalize_cas(f"CAS{cas}")
        if cas_n:
            x.cases.append(f"CAS{cas_n}")
    for hmdb in extra.get("HMDB") or []:
        h = str(hmdb).strip().upper().replace("HMDB:", "")
        x.hmdbs.append(h if h.startswith("HMDB") else f"HMDB{h}")
    for chebi in extra.get("CHEBI") or []:
        c = str(chebi).strip().upper().replace("CHEBI:", "")
        x.chebis.append(c if c.startswith("CHEBI") else f"CHEBI{c}")

    x.cids = list(dict.fromkeys(x.cids))
    x.sids = list(dict.fromkeys(x.sids))
    x.cases = list(dict.fromkeys(x.cases))
    x.hmdbs = list(dict.fromkeys(x.hmdbs))
    x.chebis = list(dict.fromkeys(x.chebis))
    return x


def resolve_metabolite_input(
    query: str,
    index: GinMetaboliteIndex,
    *,
    pubchem: PubChemMetaboliteResolver | None = None,
    allow_online: bool = True,
) -> MetaboliteResolveResult:
    """Resolve user metabolite id/name to a GIN node (GINv2 UC_ID order)."""
    q = query.strip()
    if not q:
        return MetaboliteResolveResult(query=query, gin_node=None, method="empty")

    hit = index.lookup(q)
    if hit:
        return MetaboliteResolveResult(query=q, gin_node=hit, method="gin_local")

    # Offline KEGG compound → prefer CID bridges are applied by callers via
    # ``align_kegg_xrefs_to_gin``; here keep a light local KEGGC fallback.
    kegg = normalize_kegg_compound_id(q)
    if kegg:
        for cand in (f"KEGGC{kegg[1:]}", kegg):
            hit = index.lookup(cand)
            if hit:
                return MetaboliteResolveResult(
                    query=q, gin_node=hit, method="kegg_local", matched_alias=cand
                )

    if not allow_online or pubchem is None:
        return MetaboliteResolveResult(
            query=q,
            gin_node=None,
            method="unresolved",
            warnings=["no local GIN match; online lookup disabled"],
        )

    # SID##### → CID via PubChem substance API (KEGG pubchem conv is SID)
    sid_m = _SID_NODE_RE.match(q) or (
        _SID_NODE_RE.match(f"SID{q}") if q.isdigit() else None
    )
    if sid_m:
        for cid in pubchem.cids_for_sid(sid_m.group("sid"))[:3]:
            hit = index.lookup(f"CID{cid}")
            if hit:
                return MetaboliteResolveResult(
                    query=q,
                    gin_node=hit,
                    method="pubchem_sid_cid",
                    pubchem_cid=cid,
                    matched_alias=f"CID{cid}",
                )

    cids: list[int] = []
    cas = _normalize_cas(q)
    if cas:
        cids = pubchem.cids_for_name(cas)
    if not cids:
        cids = pubchem.cids_for_name(q)

    for cid in cids[:3]:
        for alias in (f"CID{cid}", str(cid)):
            hit = index.lookup(alias)
            if hit:
                return MetaboliteResolveResult(
                    query=q,
                    gin_node=hit,
                    method="pubchem_cid",
                    pubchem_cid=cid,
                    matched_alias=alias,
                )
        for syn in pubchem.synonyms_for_cid(cid):
            hit = index.lookup(syn)
            if hit:
                return MetaboliteResolveResult(
                    query=q,
                    gin_node=hit,
                    method="pubchem_synonym",
                    pubchem_cid=cid,
                    matched_alias=syn,
                )

    return MetaboliteResolveResult(
        query=q,
        gin_node=None,
        method="pubchem_miss",
        pubchem_cid=cids[0] if cids else None,
        warnings=["PubChem hit but no matching GIN metabolite node"],
    )

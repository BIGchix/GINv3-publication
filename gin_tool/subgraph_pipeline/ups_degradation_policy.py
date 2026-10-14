"""UPS degradation policy: scale weak E3 edges; skip global deg for strong ESI."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Literal, Sequence

Strength = Literal["strong", "weak"]

_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_E3_MEMBERS = _ROOT / "seedLists/ups/e3_members.tsv"
DEFAULT_ESI_EVIDENCE = _ROOT / "seedLists/ups/esi_evidence.tsv"
DEFAULT_OVERRIDE = _ROOT / "seedLists/ups/e3_pair_override.tsv"
DEFAULT_NON_E3_ACTIVITY = _ROOT / "seedLists/ups/non_e3_activity_inhibitors.tsv"
DEFAULT_NON_E3_ACTIVITY_SCALE = 0.1

GENE_ALIASES: dict[str, tuple[str, ...]] = {
    "COP1": ("RFWD2", "COP1"),
    "RFWD2": ("RFWD2", "COP1"),
    "KMT5A": ("KMT5A", "SETD8"),
    "SETD8": ("KMT5A", "SETD8"),
}

# Deprecated pair patch; prefer non_e3_activity_inhibitors.tsv (gene-level rule).
NON_UPS_INHIBITOR_SCALES: dict[tuple[str, str], float] = {}

# Default pairs to drop for TP53 production-bottleneck ablation.
DROP_TP53_INHIBITOR_PAIRS: tuple[tuple[str, str], ...] = (
    ("MDM2", "TP53"),
    ("COP1", "TP53"),
    ("SMYD2", "TP53"),
)


def _norm(sym: str) -> str:
    return (sym or "").strip().upper()


def _expand_aliases(gene: str) -> set[str]:
    g = _norm(gene)
    out = {g}
    for canon, aliases in GENE_ALIASES.items():
        if g == canon or g in {_norm(a) for a in aliases}:
            out.update(_norm(a) for a in aliases)
            out.add(_norm(canon))
    return out


def _gene_from_protein_place(name: str) -> str | None:
    if name.endswith(".protein"):
        return name[: -len(".protein")]
    return None


def _is_protein_monomer_place(name: str) -> bool:
    return name.endswith(".protein") and ";" not in name


def _virtual_e3_substrate(name: str) -> tuple[str, str] | None:
    """Parse ``E3.protein;SUB.protein`` (order may vary)."""
    if ";" not in name:
        return None
    parts = name.split(";")
    prots = [p for p in parts if p.endswith(".protein") and p.count(".") == 1]
    if len(prots) < 2:
        # allow SID-prefixed extras: take two *.protein tokens
        prots = [p for p in parts if p.endswith(".protein")]
    if len(prots) != 2:
        return None
    g0 = _gene_from_protein_place(prots[0])
    g1 = _gene_from_protein_place(prots[1])
    if not g0 or not g1:
        return None
    return g0, g1


@dataclass
class EsiEvidenceTables:
    e3_members: frozenset[str]
    # (e3, substrate) -> strength
    pair_strength: dict[tuple[str, str], Strength]
    pair_flags: dict[tuple[str, str], str] = field(default_factory=dict)

    def is_e3(self, gene: str) -> bool:
        return bool(_expand_aliases(gene) & self.e3_members)

    def strength_for(self, e3: str, substrate: str) -> Strength | None:
        """Return strong/weak if e3 is a catalog E3; None if not an E3 pair."""
        if not self.is_e3(e3):
            return None
        for a in _expand_aliases(e3):
            for b in _expand_aliases(substrate):
                if (a, b) in self.pair_strength:
                    return self.pair_strength[(a, b)]
        # E3 in catalog but pair not in gold ESI → weak
        return "weak"

    def flags_for(self, e3: str, substrate: str) -> str:
        for a in _expand_aliases(e3):
            for b in _expand_aliases(substrate):
                if (a, b) in self.pair_flags:
                    return self.pair_flags[(a, b)]
        if self.is_e3(e3):
            return "kegg_e3_only"
        return "none"


def load_esi_evidence_tables(
    *,
    e3_members_tsv: str | Path | None = None,
    esi_evidence_tsv: str | Path | None = None,
    override_tsv: str | Path | None = None,
) -> EsiEvidenceTables:
    members_path = Path(e3_members_tsv) if e3_members_tsv else DEFAULT_E3_MEMBERS
    evidence_path = Path(esi_evidence_tsv) if esi_evidence_tsv else DEFAULT_ESI_EVIDENCE
    over_path = Path(override_tsv) if override_tsv is not None else DEFAULT_OVERRIDE

    members: set[str] = set()
    with members_path.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = _norm(row.get("gene") or "")
            if g:
                members |= _expand_aliases(g)

    pair_strength: dict[tuple[str, str], Strength] = {}
    pair_flags: dict[tuple[str, str], str] = {}
    with evidence_path.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            e3 = _norm(row.get("e3") or "")
            sub = _norm(row.get("substrate") or "")
            if not e3 or not sub:
                continue
            st = (row.get("strength") or "strong").strip().lower()
            strength: Strength = "strong" if st == "strong" else "weak"
            flags = (row.get("evidence_flags") or "").strip() or "none"
            for a in _expand_aliases(e3):
                for b in _expand_aliases(sub):
                    pair_strength[(a, b)] = strength
                    pair_flags[(a, b)] = flags

    if over_path.is_file():
        with over_path.open(newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                e3 = _norm(row.get("e3") or "")
                sub = _norm(row.get("substrate") or "")
                st = (row.get("strength") or "").strip().lower()
                if not e3 or not sub or st not in ("strong", "weak"):
                    continue
                for a in _expand_aliases(e3):
                    for b in _expand_aliases(sub):
                        pair_strength[(a, b)] = st  # type: ignore[assignment]
                        pair_flags[(a, b)] = "override"

    return EsiEvidenceTables(
        e3_members=frozenset(members),
        pair_strength=pair_strength,
        pair_flags=pair_flags,
    )


@dataclass
class UpsInboundHit:
    e3: str
    substrate: str
    strength: Strength
    evidence_flags: str
    arc_e3_to_v: tuple[str, str, float] | None = None
    arc_v_to_sub: tuple[str, str, float] | None = None


@dataclass
class UpsPolicyStats:
    n_in: int = 0
    n_out: int = 0
    n_ups_pairs: int = 0
    n_strong_pairs: int = 0
    n_weak_pairs: int = 0
    n_arcs_scaled: int = 0
    n_non_ups_arcs_scaled: int = 0
    non_ups_scaled_pairs: list[str] = field(default_factory=list)
    weak_scale: float = 0.3
    strong_e3_inhibition_weight: float = 1.0
    """Target effective |w_in×w_out| for strong ESI inhibitory legs (via √w on each arc)."""
    n_strong_arcs_scaled: int = 0
    skip_protein_places: list[str] = field(default_factory=list)
    pair_records: list[dict] = field(default_factory=list)
    n_strong_both_gold: int = 0
    n_strong_ubinet_only: int = 0
    n_strong_ubibrowser_only: int = 0
    n_weak_kegg_only: int = 0

    def to_dict(self) -> dict:
        return asdict(self)


def _classify_ups_hits(
    gin_edges: Sequence[tuple[str, str, float]],
    tables: EsiEvidenceTables,
) -> list[UpsInboundHit]:
    """Find E3.protein → V → SUB.protein inhibitory paths."""
    # index negative arcs
    neg = [(s, t, float(w)) for s, t, w in gin_edges if float(w) < 0]
    hits: list[UpsInboundHit] = []
    seen: set[tuple[str, str]] = set()

    for s, t, w in neg:
        if not _is_protein_monomer_place(t):
            continue
        sub_gene = _gene_from_protein_place(t)
        if not sub_gene:
            continue
        parsed = _virtual_e3_substrate(s)
        if parsed is None:
            continue
        g0, g1 = parsed
        # substrate must be one side; E3 the other
        if _norm(g0) == _norm(sub_gene):
            e3_gene, other = g1, g0
        elif _norm(g1) == _norm(sub_gene):
            e3_gene, other = g0, g1
        else:
            continue
        if _norm(other) != _norm(sub_gene):
            continue
        strength = tables.strength_for(e3_gene, sub_gene)
        if strength is None:
            continue
        key = (_norm(e3_gene), _norm(sub_gene))
        if key in seen:
            continue
        seen.add(key)
        # find E3 → V arc
        e3_place = f"{e3_gene}.protein"
        # gene symbol case: use actual place from edges
        arc_e3_v = None
        for s2, t2, w2 in neg:
            if t2 != s:
                continue
            if _gene_from_protein_place(s2) and _norm(
                _gene_from_protein_place(s2) or ""
            ) in _expand_aliases(e3_gene):
                arc_e3_v = (s2, t2, float(w2))
                e3_place = s2
                break
        flags = tables.flags_for(e3_gene, sub_gene)
        hits.append(
            UpsInboundHit(
                e3=_gene_from_protein_place(e3_place) or e3_gene,
                substrate=sub_gene,
                strength=strength,
                evidence_flags=flags,
                arc_e3_to_v=arc_e3_v,
                arc_v_to_sub=(s, t, float(w)),
            )
        )
    return hits


def classify_ups_inbound_edges(
    gin_edges: Sequence[tuple[str, str, float]],
    tables: EsiEvidenceTables | None = None,
) -> list[UpsInboundHit]:
    return _classify_ups_hits(gin_edges, tables or load_esi_evidence_tables())


def substrates_skip_global_protein_degradation(
    hits: Iterable[UpsInboundHit],
) -> set[str]:
    return {f"{h.substrate}.protein" for h in hits if h.strength == "strong"}


def load_non_e3_activity_inhibitors(
    path: str | Path | None = None,
    *,
    default_scale: float = DEFAULT_NON_E3_ACTIVITY_SCALE,
) -> dict[str, float]:
    """Load gene → scale for non-E3 activity-repressing enzymes (not UPS E3s)."""
    p = Path(path) if path else DEFAULT_NON_E3_ACTIVITY
    out: dict[str, float] = {}
    if not p.is_file():
        return out
    with p.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = _norm(row.get("gene") or "")
            if not g or g.startswith("#"):
                continue
            try:
                sc = float(row.get("scale") or default_scale)
            except ValueError:
                sc = float(default_scale)
            for a in _expand_aliases(g):
                out[a] = sc
    return out


def _collect_non_e3_activity_inhibitor_arc_scales(
    gin_edges: Sequence[tuple[str, str, float]],
    *,
    gene_scales: dict[str, float] | None = None,
    e3_members: frozenset[str] | None = None,
    pair_scales: dict[tuple[str, str], float] | None = None,
) -> dict[tuple[str, str], float]:
    """Map (src,tgt) → scale for inhibitory drains catalyzed by non-E3 activity enzymes.

    Applies to negative arcs onto protein substrates when the catalyst gene is in
    ``gene_scales`` and is **not** an E3 catalog member (UPS edges stay under UPS rules).
    Optional ``pair_scales`` keeps legacy (enz, sub) overrides.
    """
    genes = gene_scales if gene_scales is not None else load_non_e3_activity_inhibitors()
    e3s = e3_members  # None → do not filter by E3; caller should pass tables.e3_members
    want_pairs: dict[tuple[str, str], float] = {}
    cfg_pairs = pair_scales if pair_scales is not None else NON_UPS_INHIBITOR_SCALES
    for (enz, sub), sc in (cfg_pairs or {}).items():
        for a in _expand_aliases(enz):
            for b in _expand_aliases(sub):
                want_pairs[(a, b)] = float(sc)

    def _enz_scale(enz: str, sub: str) -> float | None:
        en = _norm(enz)
        su = _norm(sub)
        if (en, su) in want_pairs:
            return want_pairs[(en, su)]
        if e3s is not None and (_expand_aliases(enz) & set(e3s)):
            return None  # UPS / E3 catalog — do not treat as activity inhibitor
        for a in _expand_aliases(enz):
            if a in genes:
                return float(genes[a])
        return None

    arc_scale: dict[tuple[str, str], float] = {}
    for s, t, w in gin_edges:
        if float(w) >= 0:
            continue
        eg = _gene_from_protein_place(s) if _is_protein_monomer_place(s) else None
        sg = _gene_from_protein_place(t) if _is_protein_monomer_place(t) else None
        if eg and sg:
            sc = _enz_scale(eg, sg)
            if sc is not None:
                arc_scale[(s, t)] = sc
                continue
        if not _is_protein_monomer_place(t):
            continue
        sub_gene = _gene_from_protein_place(t)
        parsed = _virtual_e3_substrate(s)
        if not sub_gene or parsed is None:
            continue
        g0, g1 = parsed
        if _norm(g0) == _norm(sub_gene):
            enz_gene = g1
        elif _norm(g1) == _norm(sub_gene):
            enz_gene = g0
        else:
            continue
        sc = _enz_scale(enz_gene, sub_gene)
        if sc is not None:
            arc_scale[(s, t)] = sc
    return arc_scale


# Back-compat alias
def _collect_non_ups_inhibitor_arc_scales(
    gin_edges: Sequence[tuple[str, str, float]],
    scales: dict[tuple[str, str], float] | None = None,
) -> dict[tuple[str, str], float]:
    return _collect_non_e3_activity_inhibitor_arc_scales(
        gin_edges, pair_scales=scales, gene_scales=load_non_e3_activity_inhibitors()
    )


@dataclass
class NonE3ActivityScaleStats:
    n_in: int = 0
    n_out: int = 0
    n_arcs_scaled: int = 0
    default_scale: float = DEFAULT_NON_E3_ACTIVITY_SCALE
    scaled_arcs: list[str] = field(default_factory=list)
    genes_hit: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def apply_non_e3_activity_inhibitor_scale(
    gin_edges: Sequence[tuple[str, str, float]],
    *,
    gene_list_tsv: str | Path | None = None,
    default_scale: float = DEFAULT_NON_E3_ACTIVITY_SCALE,
    e3_members: frozenset[str] | None = None,
    pair_scales: dict[tuple[str, str], float] | None = None,
) -> tuple[list[tuple[str, str, float]], NonE3ActivityScaleStats]:
    """Scale inhibitory protein drains from non-E3 activity-repressor enzymes."""
    genes = load_non_e3_activity_inhibitors(gene_list_tsv, default_scale=default_scale)
    if e3_members is None:
        try:
            e3_members = load_esi_evidence_tables().e3_members
        except FileNotFoundError:
            e3_members = frozenset()
    arc_sc = _collect_non_e3_activity_inhibitor_arc_scales(
        gin_edges,
        gene_scales=genes,
        e3_members=e3_members,
        pair_scales=pair_scales,
    )
    stats = NonE3ActivityScaleStats(
        n_in=len(gin_edges), default_scale=float(default_scale)
    )
    hit_genes: set[str] = set()
    out: list[tuple[str, str, float]] = []
    for s, t, w in gin_edges:
        wf = float(w)
        if (s, t) in arc_sc:
            sc = float(arc_sc[(s, t)])
            wf = wf * sc
            stats.n_arcs_scaled += 1
            stats.scaled_arcs.append(f"{s}->{t}*{sc}")
            # best-effort gene tag
            eg = _gene_from_protein_place(s) if _is_protein_monomer_place(s) else None
            if eg:
                hit_genes.add(_norm(eg))
            else:
                parsed = _virtual_e3_substrate(s)
                if parsed:
                    g0, g1 = parsed
                    sub = _gene_from_protein_place(t)
                    if sub and _norm(g0) == _norm(sub):
                        hit_genes.add(_norm(g1))
                    elif sub and _norm(g1) == _norm(sub):
                        hit_genes.add(_norm(g0))
        out.append((s, t, wf))
    stats.n_out = len(out)
    stats.genes_hit = sorted(hit_genes)
    return out, stats


def apply_ups_edge_scale(
    gin_edges: Sequence[tuple[str, str, float]],
    *,
    tables: EsiEvidenceTables | None = None,
    weak_scale: float = 0.3,
    strong_e3_inhibition_weight: float = 1.0,
    e3_members_tsv: str | Path | None = None,
    esi_evidence_tsv: str | Path | None = None,
    override_tsv: str | Path | None = None,
    non_ups_inhibitor_scales: dict[tuple[str, str], float] | None = None,
) -> tuple[list[tuple[str, str, float]], UpsPolicyStats, set[str]]:
    """Scale weak UPS arcs; return edges, stats, skip_protein_places.

    Non-E3 activity-inhibitor scaling is handled separately by
    ``apply_non_e3_activity_inhibitor_scale`` (gene list rule).
    ``non_ups_inhibitor_scales`` is accepted for API compat but ignored here.
    """
    del non_ups_inhibitor_scales  # moved to apply_non_e3_activity_inhibitor_scale
    tab = tables or load_esi_evidence_tables(
        e3_members_tsv=e3_members_tsv,
        esi_evidence_tsv=esi_evidence_tsv,
        override_tsv=override_tsv,
    )
    hits = _classify_ups_hits(gin_edges, tab)
    scale = float(weak_scale)
    strong_w = float(strong_e3_inhibition_weight)
    stats = UpsPolicyStats(
        n_in=len(gin_edges),
        weak_scale=scale,
        strong_e3_inhibition_weight=strong_w,
    )
    stats.n_ups_pairs = len(hits)

    scale_arcs: set[tuple[str, str]] = set()
    strong_scale_arcs: set[tuple[str, str]] = set()
    strong_arc_factor = strong_w**0.5 if strong_w > 0.0 and strong_w != 1.0 else 1.0
    for h in hits:
        rec = {
            "e3": h.e3,
            "substrate": h.substrate,
            "strength": h.strength,
            "evidence_flags": h.evidence_flags,
        }
        stats.pair_records.append(rec)
        if h.strength == "strong":
            stats.n_strong_pairs += 1
            fl = h.evidence_flags
            if fl == "both_gold":
                stats.n_strong_both_gold += 1
            elif fl == "ubinet":
                stats.n_strong_ubinet_only += 1
            elif fl in ("ubibrowser_known", "override"):
                stats.n_strong_ubibrowser_only += 1
            if strong_arc_factor != 1.0:
                if h.arc_e3_to_v is not None:
                    strong_scale_arcs.add((h.arc_e3_to_v[0], h.arc_e3_to_v[1]))
                if h.arc_v_to_sub is not None:
                    strong_scale_arcs.add((h.arc_v_to_sub[0], h.arc_v_to_sub[1]))
        else:
            stats.n_weak_pairs += 1
            if h.evidence_flags == "kegg_e3_only":
                stats.n_weak_kegg_only += 1
            if h.arc_e3_to_v is not None:
                scale_arcs.add((h.arc_e3_to_v[0], h.arc_e3_to_v[1]))
            if h.arc_v_to_sub is not None:
                scale_arcs.add((h.arc_v_to_sub[0], h.arc_v_to_sub[1]))

    out: list[tuple[str, str, float]] = []
    for s, t, w in gin_edges:
        wf = float(w)
        if (s, t) in scale_arcs and scale != 1.0:
            wf = wf * scale
            stats.n_arcs_scaled += 1
        elif (s, t) in strong_scale_arcs and strong_arc_factor != 1.0:
            wf = -abs(wf) * strong_arc_factor
            stats.n_strong_arcs_scaled += 1
        out.append((s, t, wf))

    skip = substrates_skip_global_protein_degradation(hits)
    stats.skip_protein_places = sorted(skip)
    stats.n_out = len(out)
    return out, stats, skip


@dataclass
class DropInhibitionStats:
    n_in: int = 0
    n_out: int = 0
    n_dropped: int = 0
    dropped_arcs: list[str] = field(default_factory=list)
    pairs: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def drop_protein_inhibition_pairs(
    gin_edges: Sequence[tuple[str, str, float]],
    pairs: Sequence[tuple[str, str]] | None = None,
) -> tuple[list[tuple[str, str, float]], DropInhibitionStats]:
    """Remove inhibitory arcs for configured (enzyme, substrate) protein pairs.

    Drops both ``ENZ.protein → ENZ.protein;SUB.protein`` (if w<0) and
    ``…;ENZ…;SUB… → SUB.protein`` (if w<0), matching MDM2/COP1/SMYD2⊣TP53 encoding.
    """
    want = {
        (_norm(a), _norm(b))
        for enz, sub in (pairs or DROP_TP53_INHIBITOR_PAIRS)
        for a in _expand_aliases(enz)
        for b in _expand_aliases(sub)
    }
    stats = DropInhibitionStats(
        n_in=len(gin_edges),
        pairs=[f"{e}->{s}" for e, s in (pairs or DROP_TP53_INHIBITOR_PAIRS)],
    )

    def _matches_pair(enz: str, sub: str) -> bool:
        return (_norm(enz), _norm(sub)) in want

    out: list[tuple[str, str, float]] = []
    for s, t, w in gin_edges:
        wf = float(w)
        drop = False
        if wf < 0:
            # monomer → monomer
            eg = _gene_from_protein_place(s) if _is_protein_monomer_place(s) else None
            sg = _gene_from_protein_place(t) if _is_protein_monomer_place(t) else None
            if eg and sg and _matches_pair(eg, sg):
                drop = True
            # complex → substrate
            if not drop and _is_protein_monomer_place(t):
                sub_gene = _gene_from_protein_place(t)
                parsed = _virtual_e3_substrate(s)
                if sub_gene and parsed is not None:
                    g0, g1 = parsed
                    if _norm(g0) == _norm(sub_gene) and _matches_pair(g1, sub_gene):
                        drop = True
                    elif _norm(g1) == _norm(sub_gene) and _matches_pair(g0, sub_gene):
                        drop = True
            # enzyme → complex (negative only; SMYD2 often +1 here — leave positive)
            if not drop and _is_protein_monomer_place(s):
                eg = _gene_from_protein_place(s)
                parsed = _virtual_e3_substrate(t)
                if eg and parsed is not None:
                    g0, g1 = parsed
                    if _matches_pair(eg, g0) or _matches_pair(eg, g1):
                        # only if complex involves both enzyme and a configured substrate
                        for sub_c in (g0, g1):
                            if _matches_pair(eg, sub_c):
                                drop = True
                                break
        if drop:
            stats.n_dropped += 1
            stats.dropped_arcs.append(f"{s}->{t}@{wf}")
            continue
        out.append((s, t, wf))
    stats.n_out = len(out)
    return out, stats

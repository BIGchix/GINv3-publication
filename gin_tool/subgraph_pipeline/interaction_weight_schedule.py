"""Time-varying interaction weight schedules (e.g. MDM2⊣TP53 inhibition pulse)."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Sequence


@dataclass(frozen=True)
class InteractionWeightPulse:
    """Scale matching interactions' weights over blocks (Nutlin-like pulse).

    Timeline (0-based block index, applied *before* firing that block):
      [0, onset)                 → base_factor (usually 1)
      [onset, onset+hold)        → low_factor
      restore=step:
        [onset+hold, ∞)          → base_factor
      restore=ramp:
        [onset+hold, onset+hold+ramp_blocks) → linear low→base
        [onset+hold+ramp_blocks, ∞)          → base_factor
    """

    source: str = "MDM2.protein"
    target: str = "TP53.protein"
    kind: str = "inhibition"
    onset: int = 50
    hold: int = 20
    low_factor: float = 0.1
    restore: str = "step"  # step | ramp
    ramp_blocks: int = 20
    base_factor: float = 1.0

    def factor_at_block(self, block: int) -> float:
        b = int(block)
        onset = max(0, int(self.onset))
        hold = max(0, int(self.hold))
        low = float(self.low_factor)
        base = float(self.base_factor)
        if b < onset:
            return base
        if b < onset + hold:
            return low
        mode = str(self.restore or "step").strip().lower()
        if mode == "step":
            return base
        if mode == "ramp":
            ramp = max(1, int(self.ramp_blocks))
            t0 = onset + hold
            if b >= t0 + ramp:
                return base
            frac = (b - t0 + 1) / float(ramp)
            frac = min(1.0, max(0.0, frac))
            return low + (base - low) * frac
        raise ValueError(f"Unknown restore mode {self.restore!r}")

    def matches(self, source: str, target: str, kind: str) -> bool:
        if self.source and source != self.source:
            return False
        if self.target and target != self.target:
            return False
        if self.kind and kind != self.kind:
            return False
        return True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any] | None) -> InteractionWeightPulse | None:
        if not d:
            return None
        return cls(
            source=str(d.get("source", "MDM2.protein")),
            target=str(d.get("target", "TP53.protein")),
            kind=str(d.get("kind", "inhibition")),
            onset=int(d.get("onset", 50)),
            hold=int(d.get("hold", 20)),
            low_factor=float(d.get("low_factor", 0.1)),
            restore=str(d.get("restore", "step")),
            ramp_blocks=int(d.get("ramp_blocks", 20)),
            base_factor=float(d.get("base_factor", 1.0)),
        )


def scale_interaction_tuples_for_block(
    interactions: Sequence[tuple],
    pulse: InteractionWeightPulse | None,
    block: int,
) -> list[tuple]:
    """Return interaction tuples with matched weights × factor(block).

    Tuple layout: (src, tgt, kind, weight, mode, tid) — see InteractionTuple.
    """
    if pulse is None:
        return list(interactions)
    fac = float(pulse.factor_at_block(block))
    if abs(fac - 1.0) < 1e-15:
        return list(interactions)
    out: list[tuple] = []
    for row in interactions:
        src, tgt, kind, weight, mode, tid = row
        if pulse.matches(str(src), str(tgt), str(kind)):
            out.append((src, tgt, kind, float(weight) * fac, mode, tid))
        else:
            out.append(row)
    return out

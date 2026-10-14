"""SPN model structures (Ruths et al. 2008 / BioLayout mEPN)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

InteractionKind = Literal["activation", "inhibition"]


@dataclass
class SPNInteraction:
    """One signaling edge u→v or u⊣v."""

    source: str
    target: str
    kind: InteractionKind
    transition_id: str = ""
    weight: float = 1.0


@dataclass
class SPNModel:
    """Cell-specific signaling Petri net."""

    places: list[str]
    interactions: list[SPNInteraction]
    place_labels: dict[str, str] = field(default_factory=dict)
    fixed_high: set[str] = field(default_factory=set)
    fixed_low: set[str] = field(default_factory=set)

    def with_perturbation(
        self,
        *,
        knockdown: set[str] | None = None,
        constitutive_high: set[str] | None = None,
    ) -> SPNModel:
        """Return perturbed copy (knockdown → fixed low, removes activations)."""
        knockdown = knockdown or set()
        constitutive_high = constitutive_high or set()
        low = set(self.fixed_low) | knockdown
        high = (set(self.fixed_high) | constitutive_high) - low
        kept: list[SPNInteraction] = []
        for edge in self.interactions:
            if edge.target in low and edge.kind == "activation":
                continue
            if edge.source in low and edge.kind == "activation":
                continue
            if edge.target in high and edge.kind == "inhibition":
                continue
            kept.append(edge)
        return SPNModel(
            places=list(self.places),
            interactions=kept,
            place_labels=dict(self.place_labels),
            fixed_high=high,
            fixed_low=low,
        )

    def add_self_inhibition_loops(self, places: list[str] | None = None) -> None:
        """Add u⊣u loops for places without explicit inhibitory incoming edges."""
        targets_with_inhibition = {
            e.target for e in self.interactions if e.kind == "inhibition"
        }
        for place in places or self.places:
            if place in targets_with_inhibition:
                continue
            if any(e.source == place and e.target == place for e in self.interactions):
                continue
            self.interactions.append(
                SPNInteraction(
                    source=place,
                    target=place,
                    kind="inhibition",
                    transition_id=f"self_{place}",
                )
            )

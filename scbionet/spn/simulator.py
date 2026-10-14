"""Signaling Petri Net simulator (Ruths et al. 2008; BioLayout SPN).

Improvements over BioLayout GUI-only workflow:
- CLI / batch runs with reproducible seeds
- uniform or truncated-normal token sampling
- optional edge weights; differential perturbation + t-test
- exports trajectories (CSV/JSON), not only end-of-block snapshots
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

from scbionet.spn.model import SPNInteraction, SPNModel

Distribution = Literal["uniform", "normal"]


@dataclass
class SimulationResult:
    place_names: list[str]
    mean_markings: np.ndarray  # (n_blocks+1, n_places)
    std_markings: np.ndarray
    n_runs: int
    n_blocks: int

    def to_dict(self) -> dict:
        return {
            "places": self.place_names,
            "n_runs": self.n_runs,
            "n_blocks": self.n_blocks,
            "mean_markings": self.mean_markings.tolist(),
            "std_markings": self.std_markings.tolist(),
        }

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    def save_csv(self, path: str | Path) -> None:
        import csv

        out = Path(path)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            header = ["block"] + self.place_names
            writer.writerow(header)
            for b in range(self.mean_markings.shape[0]):
                row = [b] + [f"{v:.4f}" for v in self.mean_markings[b]]
                writer.writerow(row)


@dataclass
class SPNSimulator:
    """Non-parametric SPN executor."""

    n_blocks: int = 20
    n_runs: int = 400
    high_tokens: int = 10
    low_tokens: int = 0
    distribution: Distribution = "uniform"
    seed: int | None = 42
    min_tokens: float = 0.0

    def _sample_amount(self, rng: np.random.Generator, upper: float) -> float:
        if upper <= 0:
            return 0.0
        hi = max(1, int(np.floor(upper)))
        if self.distribution == "normal":
            # BioLayout optional Standard Normal → truncate to [1, hi]
            val = int(round(rng.normal(loc=(1 + hi) / 2, scale=max(1.0, hi / 3))))
            val = int(np.clip(val, 1, hi))
            return float(val)
        return float(rng.integers(1, hi + 1))

    def _build_index(self, model: SPNModel) -> tuple[list[str], dict[str, int], list[SPNInteraction]]:
        places = list(model.places)
        idx = {p: i for i, p in enumerate(places)}
        return places, idx, list(model.interactions)

    def _initial_marking(self, model: SPNModel, places: list[str], idx: dict[str, int]) -> np.ndarray:
        m = np.zeros(len(places), dtype=np.float64)
        for p in model.fixed_high:
            if p in idx:
                m[idx[p]] = self.high_tokens
        for p in model.fixed_low:
            if p in idx:
                m[idx[p]] = self.low_tokens
        return m

    def _apply_stimulus(
        self,
        marking: np.ndarray,
        idx: dict[str, int],
        stimulus: dict[str, float],
    ) -> None:
        for place, tokens in stimulus.items():
            if place in idx:
                marking[idx[place]] = max(marking[idx[place]], tokens)

    def _fire_transition(
        self,
        marking: np.ndarray,
        idx: dict[str, int],
        edge: SPNInteraction,
        rng: np.random.Generator,
    ) -> None:
        if edge.source not in idx or edge.target not in idx:
            return
        u = idx[edge.source]
        v = idx[edge.target]
        w = max(edge.weight, 1e-8)
        src_tokens = marking[u] * w
        if edge.kind == "activation":
            delta = self._sample_amount(rng, src_tokens)
            marking[v] += delta
        else:
            delta = self._sample_amount(rng, src_tokens)
            marking[v] = max(self.min_tokens, marking[v] - delta)

    def run(
        self,
        model: SPNModel,
        *,
        initial_marking: dict[str, float] | None = None,
        stimulus: dict[str, float] | None = None,
    ) -> SimulationResult:
        places, idx, interactions = self._build_index(model)
        if not interactions:
            raise ValueError("SPN model has no interactions")

        rng = np.random.default_rng(self.seed)
        n_places = len(places)
        n_trans = len(interactions)
        block_samples = np.zeros((self.n_runs, self.n_blocks + 1, n_places), dtype=np.float64)

        for run in range(self.n_runs):
            m = self._initial_marking(model, places, idx)
            if initial_marking:
                for p, val in initial_marking.items():
                    if p in idx:
                        m[idx[p]] = val
            if stimulus:
                self._apply_stimulus(m, idx, stimulus)

            block_samples[run, 0] = m.copy()
            for block in range(self.n_blocks):
                order = rng.permutation(n_trans)
                for t_idx in order:
                    self._fire_transition(m, idx, interactions[t_idx], rng)
                block_samples[run, block + 1] = m.copy()

        mean_m = block_samples.mean(axis=0)
        std_m = block_samples.std(axis=0, ddof=1) if self.n_runs > 1 else np.zeros_like(mean_m)
        return SimulationResult(
            place_names=places,
            mean_markings=mean_m,
            std_markings=std_m,
            n_runs=self.n_runs,
            n_blocks=self.n_blocks,
        )

    def differential_run(
        self,
        control: SPNModel,
        perturbed: SPNModel,
        *,
        stimulus: dict[str, float] | None = None,
    ) -> tuple[SimulationResult, SimulationResult, dict]:
        """Control vs perturbed with per-place t-test at final block."""
        res_c = self.run(control, stimulus=stimulus)
        res_p = self.run(perturbed, stimulus=stimulus)
        # re-run paired runs for stats (same seeds offset)
        rng = np.random.default_rng(self.seed)
        places_c, idx_c, inter_c = self._build_index(control)
        places_p, idx_p, inter_p = self._build_index(perturbed)
        common = sorted(set(places_c) & set(places_p))
        final_c = []
        final_p = []
        for run in range(self.n_runs):
            seed_c = int(rng.integers(0, 2**31 - 1))
            sim_c = SPNSimulator(
                n_blocks=self.n_blocks,
                n_runs=1,
                high_tokens=self.high_tokens,
                low_tokens=self.low_tokens,
                distribution=self.distribution,
                seed=seed_c,
            )
            sim_p = SPNSimulator(
                n_blocks=self.n_blocks,
                n_runs=1,
                high_tokens=self.high_tokens,
                low_tokens=self.low_tokens,
                distribution=self.distribution,
                seed=seed_c,
            )
            rc = sim_c.run(control, stimulus=stimulus)
            rp = sim_p.run(perturbed, stimulus=stimulus)
            ic = {p: rc.mean_markings[-1, rc.place_names.index(p)] for p in common}
            ip = {p: rp.mean_markings[-1, rp.place_names.index(p)] for p in common}
            final_c.append([ic[p] for p in common])
            final_p.append([ip[p] for p in common])

        fc = np.asarray(final_c)
        fp = np.asarray(final_p)
        from scipy import stats

        pvals = {}
        for j, p in enumerate(common):
            if fc[:, j].std() < 1e-12 and fp[:, j].std() < 1e-12:
                pvals[p] = 1.0
            else:
                _, pv = stats.ttest_ind(fc[:, j], fp[:, j], equal_var=False)
                pvals[p] = float(pv)
        return res_c, res_p, {"places": common, "p_values": pvals, "delta_mean": (fp - fc).mean(axis=0).tolist()}

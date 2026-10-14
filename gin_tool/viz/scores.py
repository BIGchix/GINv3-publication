"""RWR score distribution plots."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from gin_tool.mapping import gin_entity_layer


def render_score_distribution(
    viz_dir: Path,
    names: list[str],
    scores: np.ndarray,
    seed_indices: set[int],
) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    layers = [gin_entity_layer(n) for n in names]
    is_seed = [i in seed_indices for i in range(len(names))]

    import pandas as pd

    df = pd.DataFrame(
        {"layer": layers, "score": scores, "is_seed": is_seed, "node": names}
    )
    nonvirt = df[df["layer"] != "virtual"]
    if nonvirt.empty:
        nonvirt = df

    fig, ax = plt.subplots(figsize=(9, 5))
    order = sorted(nonvirt["layer"].unique())
    if len(nonvirt) >= 12:
        sns.violinplot(
            data=nonvirt,
            x="layer",
            y="score",
            order=order,
            inner="box",
            ax=ax,
            cut=0,
        )
    else:
        sns.boxplot(data=nonvirt, x="layer", y="score", order=order, ax=ax)
    seeds_df = nonvirt[nonvirt["is_seed"]]
    if not seeds_df.empty:
        sns.stripplot(
            data=seeds_df,
            x="layer",
            y="score",
            order=order,
            color="#c0392b",
            size=8,
            jitter=0.15,
            ax=ax,
            zorder=3,
        )
    ax.set_title("RWR score by molecular layer (red = seeds)")
    ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(viz_dir / "score_distribution.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def render_ecdf_null_comparison(
    viz_dir: Path,
    scores: np.ndarray,
    *,
    null_precompute_dir: Path | None = None,
    k: int | None = None,
    n_seeds: int | None = None,
    rwr_alpha: float | None = None,
) -> None:
    """Observed subgraph RWR ECDF vs null replicate subsamples (schema 1.2)."""
    import matplotlib.pyplot as plt

    x_obs = np.sort(scores)
    y_obs = np.arange(1, len(x_obs) + 1) / max(len(x_obs), 1)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x_obs, y_obs, label="observed subgraph", color="#2980b9", lw=2)

    used_null = False
    if null_precompute_dir is not None and k is not None and n_seeds is not None:
        try:
            from gin_tool.null_rwr import collect_null_rwr_ecdf_samples, null_ecdf_median_curve

            nearest_n, samples, max_scores = collect_null_rwr_ecdf_samples(
                null_precompute_dir,
                k,
                n_seeds,
                rwr_alpha=rwr_alpha,
            )
            if samples:
                used_null = True
                xmin = min(float(x_obs.min()), float(min(s.min() for s in samples)))
                xmax = max(float(x_obs.max()), float(max(s.max() for s in samples)))
                grid = np.linspace(xmin, xmax, 200)
                y_med = null_ecdf_median_curve(samples, grid)
                ax.plot(
                    grid,
                    y_med,
                    color="#e67e22",
                    lw=2,
                    ls="--",
                    label=f"null median ECDF (n_seeds={nearest_n}, n={len(samples)})",
                )
                for s in samples[: min(12, len(samples))]:
                    sg = np.sort(s)
                    yg = np.arange(1, len(sg) + 1) / max(len(sg), 1)
                    ax.plot(sg, yg, color="#e67e22", alpha=0.12, lw=0.8)
            if max_scores:
                ax.axvline(
                    float(np.median(max_scores)),
                    color="#95a5a6",
                    ls=":",
                    label="null median max score",
                )
        except (FileNotFoundError, ValueError) as exc:
            ax.text(
                0.02,
                0.02,
                f"null ECDF unavailable: {exc}",
                transform=ax.transAxes,
                fontsize=8,
            )

    ax.set_xlabel("RWR score")
    ax.set_ylabel("ECDF")
    title = "RWR score ECDF vs random-seed null"
    if not used_null:
        title += " (re-run precompute_null schema 1.2 for null curves)"
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(viz_dir / "ecdf_null_comparison.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

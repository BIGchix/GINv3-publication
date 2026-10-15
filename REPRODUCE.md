# Publication recipe (GINv3.11a) — real compute

This document rebuilds manuscript artifacts by **executing** `gin_tool` / **SigPN**
(signaling Petri net — not stochastic PN). It is not a copy of `expected/`.

For **your own** seeds and analyses, use [docs/USER_GUIDE.md](docs/USER_GUIDE.md)
and `examples/01_custom_context.sh`. For cell ChIP allowlists and custom
simulations (beyond the H1 B0 Pulse special case), see
[docs/CHIP_ALLOWLIST.md](docs/CHIP_ALLOWLIST.md) and
[docs/SIGPN_CUSTOM_SIMULATIONS.md](docs/SIGPN_CUSTOM_SIMULATIONS.md).

Determinism: `export PYTHONHASHSEED=0`. Publication freeze used `n_jobs=20`.

The GIN prior used throughout is **`data/GINv3.11a.sifi.gz`** (read directly).
It ships in the git tree; GitHub Releases do not re-upload it as an asset.
Standalone URL (tag `v0.2.0`):
https://github.com/BIGchix/GINv3-publication/raw/v0.2.0/data/GINv3.11a.sifi.gz
— see [data/README.md](data/README.md).

## 0. Prepare the prior

```bash
python scripts/reproduce_publication.py --steps prepare
```

Verifies bundled `data/GINv3.11a.sifi.gz` (sha256 `f9bb23c4…`). No unpack.

## 1. Stage A (recompute)

```bash
python scripts/reproduce_publication.py --steps stage_a --n-jobs 20
# or: bash examples/02_publication_stage_a.sh
```

Lock: **11,368** edges · sha256 `350b7cb9…`.  
Then: `python scripts/assert_p53_pub_gin11a_lock.py`.

> Bit-match to the lock is expected on the original Linux freeze host. Other
> platforms may produce a nearby topology — see
> `docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md`.

## 2. Fig. 3e hop series

```bash
python scripts/reproduce_publication.py --steps fig3e --n-jobs 20
```

## 3. Fig. 4a–d (six lines, noprune)

```bash
python scripts/reproduce_publication.py --steps fig4cd --n-jobs 20
python scripts/assert_p53_pub_fig4cd_noprune_lock.py
```

## 4. H1 prune + B0 Pulse（发表特例）

子图路径：Stage A → ChIP allowlist prune → H1细胞边表 → SigPN B0（MDM2⊣TP53）。  
细节与「可复用 / 需自写」划分见 [docs/SIGPN_CUSTOM_SIMULATIONS.md](docs/SIGPN_CUSTOM_SIMULATIONS.md)。

```bash
python scripts/reproduce_publication.py --steps prune_h1,pulse --n-jobs 20 --seed 43
```

## Full pipeline

```bash
python scripts/reproduce_publication.py --steps all --n-jobs 20
```

## About `expected/` and `--from-expected`

| Path / flag | Role |
|-------------|------|
| `expected/` | Gold fingerprints + reference PDFs for **post-run asserts** and packaging |
| `--from-expected` | **CI / offline integrity only** — installs references without computing. **Not** the reproduction path for reviewers |

## Hardware

| Step | Typical |
|------|---------|
| Stage A / Fig. 3e | 96–128 GB RAM, 20 cores, hours |
| Fig. 4a–d | 6 × 200×200 SigPN |
| Pulse | 350×350 |

`--smoke` uses tiny MC sizes and will **not** match lock numbers.

## Out of scope

- Fig. 1–2 artwork
- Rebuilding GINv3.11a from earlier GIN releases
- Re-downloading raw ENCODE/ChIP peaks (precomputed runtime files are included)

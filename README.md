# GINv3-publication

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

Multi-omics **Global Integrative Network (GIN)** context extraction and **SigPN**
(**signaling** Petri net — not stochastic PN) simulation for mechanism-oriented
analyses. This repository is the manuscript **Code Availability** package:
an installable toolkit (`gin-multiomics`) plus a recipe that **recomputes** the
p53 / GINv3.11a publication numbers.

The GIN prior ships as **`data/GINv3.11a.sifi.gz`** (~75 MB). The toolkit and
publication recipe read it **directly** — no unpack step is required.

---

## Citation

Upon publication, please cite the accompanying paper (title / journal / DOI to
be updated here). Until then, cite this repository and the manuscript that
points to it under Code Availability.

---

## Overview

- Extract a **mechanistic context** subgraph from the GINv3 prior given gene
  (or GIN node) seeds.
- Optionally prune TF→promoter edges with a **cell-specific ChIP allowlist**
  to obtain a cell-line network.
- Run **SigPN** Monte-Carlo simulations for signaling interventions.
- **Recompute** the manuscript Stage A / prune / Pulse pipeline and assert
  against published fingerprints (locks under `docs/locks/`).

---

## Requirements

| Item | Notes |
|------|--------|
| OS | Linux recommended for bit-identical Stage A; macOS works for toolkit use |
| Python | ≥ 3.9 |
| Disk | Bundled `data/GINv3.11a.sifi.gz` (~75 MB) |
| RAM / CPU | Context with `--context-n-jobs 20` is cluster-scale; small demos run on a laptop |

---

## Installation

```bash
git clone https://github.com/BIGchix/GINv3-publication.git
cd GINv3-publication
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip && pip install -e ".[sim]"
export PYTHONHASHSEED=0
python scripts/reproduce_publication.py --steps prepare
# → verifies data/GINv3.11a.sifi.gz (sha256 f9bb23c4…)
```

After install, the CLI is available as:

```bash
python -m gin_tool.cli -h
# or: gin-multiomics -h
```

Default `--sifi` is `data/GINv3.11a.sifi.gz` (included in the clone).
Standalone download (Release does not attach the ~75 MB asset):
[raw v0.2.0](https://github.com/BIGchix/GINv3-publication/raw/v0.2.0/data/GINv3.11a.sifi.gz)
· [data/README.md](data/README.md). More detail:
[docs/USER_GUIDE.md](docs/USER_GUIDE.md).

---

## Quick start

**1. Context subgraph** (default prior = compressed GIN):

```bash
python -m gin_tool.cli \
  --seeds "TP53,MDM2,CDKN1A" \
  --out-dir gin_out/my_context \
  --context-mech-depth 2 \
  --context-n-jobs 8 \
  --node-list-dir ./nodeList \
  --no-debias
# → gin_out/my_context/subgraph_edges.tsv
```

```bash
bash examples/01_custom_context.sh
```

**2. SigPN demo** (bundled toy network; **no** GIN prior required):

```bash
bash examples/03_sigpn_demo.sh
# → gin_out/example_sigpn_demo/spn.graphml
# → gin_out/example_sigpn_demo/spn_sim_demo/…
```

See [examples/](examples/) for more scripts.

---

## Workflow

```text
seeds + data/GINv3.11a.sifi.gz
        │
        ▼
  context subgraph  ──►  (optional) ChIP prune  ──►  SigPN simulation
```

1. **Context** — seeds + GIN prior → mechanistic neighborhood (`python -m gin_tool.cli`).
2. **Cell prune** (optional) — ChIP allowlist → cell-specific edges
   ([docs/CHIP_ALLOWLIST.md](docs/CHIP_ALLOWLIST.md)).
3. **SigPN** — export the Petri net and run Monte-Carlo (demo below). The
   manuscript H1 B0 Pulse is a special-case intervention script; for other
   designs reuse the same engine
   ([docs/SIGPN_CUSTOM_SIMULATIONS.md](docs/SIGPN_CUSTOM_SIMULATIONS.md)).

### From subgraph to SigPN (executable demo)

Bundled toy data under `examples/demo_sigpn/` (TP53 / MDM2 / CDKN1A edges
excerpted from the publication H1-pruned network). Laptop-scale defaults:
`runs=2`, `blocks=10`.

```bash
# One-liner
bash examples/03_sigpn_demo.sh

# Equivalent explicit steps
mkdir -p gin_out/example_sigpn_demo
cp examples/demo_sigpn/subgraph_edges.tsv gin_out/example_sigpn_demo/
cp examples/demo_sigpn/place_levels.tsv   gin_out/example_sigpn_demo/

python examples/run_sigpn_demo.py \
  --subgraph-dir gin_out/example_sigpn_demo \
  --place-levels gin_out/example_sigpn_demo/place_levels.tsv \
  --chromhmm-emission-tsv runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv \
  --runs 2 \
  --blocks 10 \
  --seed 0
```

On your own context (after CLI or prune), point `--subgraph-dir` at a folder
that contains `subgraph_edges.tsv` (and optionally `place_levels.tsv`).

---

## Reproduce manuscript results

Recompute the publication pipeline (prior stays compressed throughout):

```bash
export PYTHONHASHSEED=0
python scripts/reproduce_publication.py --steps prepare
python scripts/reproduce_publication.py --steps stage_a --n-jobs 20
python scripts/assert_p53_pub_gin11a_lock.py

# Full pipeline (cluster-scale: prune + Pulse, etc.)
python scripts/reproduce_publication.py --steps all --n-jobs 20
```

Stage A lock target: **11,368** edges (details in `docs/locks/`).  
`expected/` is for **assert / CI**, not the primary reproduction path.  
Platform notes (e.g. Linux freeze host vs macOS): see
[docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md](docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md)
and the step-by-step guide [REPRODUCE.md](REPRODUCE.md).

---

## Documentation

| Document | Contents |
|----------|----------|
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | Getting started for custom analyses |
| [docs/CLI_REFERENCE.md](docs/CLI_REFERENCE.md) | Recommended CLI options |
| [docs/CHIP_ALLOWLIST.md](docs/CHIP_ALLOWLIST.md) | Cell ChIP allowlist format and use |
| [docs/SIGPN_CUSTOM_SIMULATIONS.md](docs/SIGPN_CUSTOM_SIMULATIONS.md) | SigPN reuse; Pulse as a special case |
| [REPRODUCE.md](REPRODUCE.md) | Full manuscript recompute recipe |
| [MANIFEST.md](MANIFEST.md) | Package contents checklist |

---

## Repository layout

| Path | Role |
|------|------|
| `gin_tool/` | Context / HME / dogma / SigPN CLI |
| `scbionet/spn/` | Vendored SigPN engine |
| `data/GINv3.11a.sifi.gz` | Frozen GIN prior (default `--sifi`) |
| `nodeList/` · `seedLists/` | Layer lists and seeds |
| `runtime/chip_gin_align/` | Example H1 allowlist |
| `runtime/omics_seed/` · `runtime/chromhmm/` | Example place levels / emissions |
| `examples/` | Context / Stage A / SigPN demos (`demo_sigpn/` toy data) |
| `scripts/reproduce_publication.py` | Manuscript driver (recomputes) |
| `expected/` | Gold fingerprints for assert |
| `docs/` | Guides, locks, experiment notes |

Historical path prefixes such as `runtime/spn_phys5/` still mean **SigPN**.

---

## License

MIT for code (see [LICENSE](LICENSE)). Network and omics artifacts: see
[NOTICE.md](NOTICE.md).

---

## Contact

- Email: [chix@big.ac.cn](mailto:chix@big.ac.cn)
- GitHub Issues on this repository (once published under Code Availability)

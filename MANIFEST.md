# Bundle contents

**GINv3-publication** — GIN multi-omics toolkit + GINv3.11a p53 publication recipe  
(lock `p53_pub_gin11a_hs0_20260929`). Simulation engine: **SigPN** (signaling Petri net).

## Product (general use)

| Area | Contents |
|------|----------|
| Code | `gin_tool/` CLI & pipeline; `scbionet/spn/` SigPN engine |
| Entry | `gin-multomics` / `python -m gin_tool.cli` |
| Docs | `USER_GUIDE.md` · `CLI_REFERENCE.md` · `CHIP_ALLOWLIST.md` · `SIGPN_CUSTOM_SIMULATIONS.md` |
| Prior | `data/GINv3.11a.sifi.gz` (read directly; optional `--materialize-sifi`) |
| Layers / seeds | `nodeList/`, `seedLists/` |
| Example H1 allowlist | `runtime/chip_gin_align/h1_wa01/` |
| Optional runtime | omics place levels, chromHMM emissions |

## Publication recipe

| Area | Contents |
|------|----------|
| Driver | `scripts/reproduce_publication.py` (**recomputes**) |
| Asserts | `scripts/assert_p53_pub_*.py` |
| Locks | `docs/locks/` |
| Stage A platform note | `docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md` |

## Reference only

| Area | Contents |
|------|----------|
| `expected/` | Gold fingerprints for assert after compute |
| `--from-expected` | CI / packaging only |

## Key fingerprints

| Artifact | sha256 prefix |
|----------|----------------|
| Bundled `GINv3.11a.sifi.gz` | `f9bb23c4…` |
| Unpacked `GINv3.11a.sifi` (optional) | `192b2cb5…` |
| Stage A edges | `350b7cb9…` (11,368) |
| H1 edges | `562d79bc…` (11,363) |

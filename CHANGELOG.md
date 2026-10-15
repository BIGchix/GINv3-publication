# Changelog

All notable changes to **GINv3-publication** / `gin-multiomics` are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] — 2026-10-15

First public Code Availability cut aligned with the GINv3.11a p53 publication
recipe. Installable as `gin-multiomics`; default GIN prior is the compressed
bundle under `data/`.

### Added

- Installable toolkit (`gin_tool`, CLI entry `gin-multiomics` /
  `python -m gin_tool.cli`) for context subgraph extraction, ChIP allowlist
  prune hooks, and SigPN (signaling Petri net) simulation.
- Vendored SigPN engine under `scbionet/spn/`.
- Bundled prior `data/GINv3.11a.sifi.gz` (~75 MB) with sidecar fingerprints
  (`f9bb23c4…` gzip; uncompressed `192b2cb5…`).
- Transparent `.sifi` / `.sifi.gz` I/O (`gin_tool.io.open_sifi_text`); pipeline
  and CLI **default to the compressed prior** (no unpack required).
- Publication driver `scripts/reproduce_publication.py` (prepare / stage_a /
  fig3e / fig4cd / prune_h1 / pulse / assert).
- Lock asserts under `scripts/assert_p53_pub_*.py` and fingerprints in
  `docs/locks/` (Stage A target **11,368** edges · `350b7cb9…` on the Linux
  freeze host).
- Examples: custom context, publication Stage A, laptop-scale SigPN demo
  (`examples/demo_sigpn/`).
- User docs: `README.md`, `docs/USER_GUIDE.md`, `CLI_REFERENCE.md`,
  `CHIP_ALLOWLIST.md`, `SIGPN_CUSTOM_SIMULATIONS.md`, `REPRODUCE.md`.
- Restored Active TF RWR null helper `scripts/active_tf_rwr_null_strength.py`
  (publication Stage A still uses the frozen ATF list by default).
- Smoke test `tests/test_sifi_gz_io.py` for plain vs gzip edge-stream parity.

### Changed

- `prepare` verifies `data/GINv3.11a.sifi.gz` only; optional
  `--materialize-sifi` remains an advanced escape hatch and does **not** change
  the default `--sifi` path.
- README / guides / examples omit unpack from the happy path; default CLI
  `--sifi` is `data/GINv3.11a.sifi.gz`.
- Clarified product naming: **GIN** = Global Integrative Network; **SigPN** =
  signaling Petri net (not stochastic PN).

### Notes for reproducers

- Bit-identical Stage A vs lock is expected on the original Linux freeze host
  (`PYTHONHASHSEED=0`, `n_jobs=20`). Other platforms may yield a nearby
  topology (see `docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md`).
- A `pip install` of the Python package ships **code** (`gin_tool` /
  `scbionet`); the large GIN prior lives in the git checkout under `data/`
  (not inside the wheel; not re-attached on the GitHub Release page).
  Download: https://github.com/BIGchix/GINv3-publication/raw/v0.2.0/data/GINv3.11a.sifi.gz

## [0.1.0] — 2026-10-13

Internal bootstrap of the Code Availability tree (pre-public packaging). Not
tagged on GitHub; listed here for lineage only.

[Unreleased]: https://github.com/BIGchix/GINv3-publication/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/BIGchix/GINv3-publication/releases/tag/v0.2.0
[0.1.0]: https://github.com/BIGchix/GINv3-publication/commits/21eee2e

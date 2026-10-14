# Third-party data notices

This package redistributes **derived** network and omics artifacts for manuscript
reproduction. Upstream licenses and attribution remain with the source databases:

- **GIN** bipartite prior (`data/GINv3.11a.sifi.gz`) — original GIN / lab build
  lineage; do not treat as a re-license of source pathway DBs.
- **ENCODE / Roadmap** chromHMM emissions and related weights under
  `runtime/chromhmm/` — follow ENCODE / Roadmap terms of use.
- **ChIP-Atlas / related TF evidence** summarized as H1 allowlists under
  `runtime/chip_gin_align/` — not raw peak beds; cite ChIP-Atlas / original
  experiments as appropriate.
- Cell-line RNA / metabolite place levels under `runtime/omics_seed/` — derived
  expression features used for SigPN seeding.

Code under `gin_tool/` and `scbionet/spn/` is MIT (see `LICENSE`) unless a file
header states otherwise.

# SigPN layer rules (v2)

> **Naming:** This document describes **SigPN** (signaling Petri net), not stochastic Petri nets.
> Historical filenames/flags may still say `spn` / `SPN`.

Reusable firing semantics for GIN → Petri → SigPN. Pathway-specific patches are out of scope; p53/H1 are acceptance cases only.

## Enable

```bash
python scripts/run_subgraph_spn_simulation.py \
  --subgraph-dir gin_out/<subgraph> \
  --layer-rules v2 \
  --stoich one_to_one \
  --tf-promoter-proteolysis 0
```

Default without `--layer-rules` remains legacy `gin_consume_v1` (`spn_v5`).

## Contract

| Layer | Source consumed? | How target grows | Degradation |
|---|---|---|---|
| Metabolic (token 1:1 / later KEGG ν) | Yes (≤ available; no mint) | Token transfer (≠ edge-weight balance) | Optional extras only — **no** default all-met degrade |
| Transcription | No (`promoter.active` read-only) | \(\Delta RNA = \alpha \min(m_{active}, M_{max})\) | RNA / active half-life **destroy** (no sink place) |
| Translation | No (RNA read-only) | \(\Delta protein = \alpha \min(m_{RNA}, M_{max})\) | Protein half-life destroy |
| Histone → chromatin | No | Regulatory modulate | Histone half-life destroy |
| TF → promoter | No (default recycle) | Regulatory modulate **via sigmoid** on activity places | Protein half-life; optional `--tf-promoter-proteolysis ε` |
| Enzyme → metabolite (e.g. PTEN) | **No** (catalyst; **token no-op**) | Chemistry via met↔met transfer only | Protein half-life destroy |
| Kinase → protein substrate | **No** (catalyst; **token no-op**) | No minting of protein copies; levels from translation / complex / degrade only | Protein half-life destroy |
| Complex assembly (`A_B`) | **Yes** | Monomer ↔ complex conserving transfer | Protein / complex half-life |
| Metabolite ↔ metabolite | Yes | Conserving transfer (± edge-balance) | Optional extras only |

**Anti-explosion contract:** REGULATORY never does unbounded `target += δ` (only logistic occupancy on activity places). CATALYST never does `target += f(enzyme)`. Copy-number growth is only transcription / translation / conserving transfer.

## Two meanings of「守恒」(do not conflate)

| 说法 | 实际在做什么 | 是否已在 v2 |
|---|---|---|
| **Token 守恒（`TRANSFER_CONSERVING`）** | 开火时从源扣 δ、给靶加 δ，且 δ≤源存量；禁止 w>1 造物 | 是 |
| **边权/度平衡（加权平均）** | 按节点 Σ\|w_in\|≈Σ\|w_out\| 降权重侧；缺计量时的拓扑先验 | 可选 `--edge-balance` |

### Out-route (first-order turnover + categorical outs)

先前用固定 `out_route_budget=1` 欠妥：边权只回答「走哪条」，总通量应是浓度的函数 \(δ_{\mathrm{tot}}=f(m)\)，不是常数。

现语义（**反应周转**，不是 sink 降解）：

- Group outs by `(source, mode, kind)` for `TRANSFER_CONSERVING` and `REGULATORY`
- Each block: sample \(δ_{\mathrm{tot}}\) with the same discrete half-life weight map as degradation (`first_order_turnover_weight`), then pick **one** out-edge with \(p_i \propto |w_i|\)
  - **TRANSFER**: move \(δ_{\mathrm{tot}}\) tokens (≤ available)
  - **REGULATORY**: one event with raw drive \(δ = δ_{\mathrm{tot}}\) (TF not consumed; logit step uses this drive)
- TF→chromatin included (weights may stay unscaled due to protection; still normalized when picking)
- Independent modes (transcription / translation / degrade / catalyst) still fire every edge
- CLI: `--out-route-half-life-blocks` (default **20**); `≤0` disables routing (falls back to independent edges)

Do **not** restore independent `δ∼m·w_i` on every out-edge (that reintroduces degree-driven explosion).

Token 守恒**不**修复 PIP2(1/30) 这类拓扑失衡。边权平衡默认建议仅 `metabolite`；若开 `protein`，**不得**降权 TF↔`*.promoter;state` / `promoter.active`（TF 回池，已在代码中保护）。

## miRNA → translation (repression)

**GIN master (`GINv3.11.sifi`)**: cross-gene `miRNA → Gene.*;translation` is
``negative`` (~1.8M). Cognate `MIR*.RNA → MIR*;translation` (+1 self-translation)
was **removed** (`scripts/drop_mir_self_translation_edges.py`; backup
`GINv3.11.sifi.pre_drop_mir_self_translation.bak`) — miRNA is not an mRNA for
protein synthesis. Both legs (`RNA→tid` and `tid→protein`) dropped (~1223 each).

**Subgraph bug (fixed)**: HME chain repair used to re-materialize R→V→R edges with
hardcoded ``weight=1.0``, wiping GIN negatives (e.g. `MIR3121→KDM6B.protein;translation`
became +1). Repair now copies signed weights from the GIN adjacency.

**SPN fire (current repression formula)** — when `kind=inhibition` or source is miRNA
on a `*;translation` transition:

\[
\begin{aligned}
\delta &\sim \alpha\cdot\mathrm{Unif}\{1,\ldots,\lfloor\min(m_{\mathrm{miR}},M_{\max})\rfloor\} \\
p &\leftarrow p - \min(p-p_{\min},\,\delta)
\end{aligned}
\]

i.e. same `production_amount` sampler as activating translation, but **subtracts**
from the protein place (ledger `translation_repress`) instead of minting. miRNA
itself is not consumed. This is a **protein-pool proxy**, not a true reduction of
cognate-mRNA translation rate (and in this histprune subgraph those miRNA
translation transitions have **no cognate mRNA input**).

Because repression does not consume miRNA, **degradation is the required outlet**:
miRNA places always get a `;degradation` arc (default half-life = mRNA's
`--rna-half-life-blocks`; override with `--mirna-half-life-blocks`).

## Metabolite upstream seeding

`met↔met` digraph is often **not** a DAG (cycles / large SCCs). Practical seeding: add tokens on **in-degree 0** sources that still have out-edges (`gin_tool/subgraph_pipeline/met_source_seed.py`, `scripts/build_met_source_place_levels.py`).

## Activity sigmoid (`promoter.active` / `promoter;state` / histone marks)

Regulatory (and conserving writes onto these places) use logit-space logistic updates so levels stay in `[0, M]`:

\[
m \leftarrow M\cdot\sigma\!\left(\mathrm{logit}(m/M) \pm k\cdot\delta/M\right)
\]

Applies to:
- `*.promoter.active`, `*.promoter;state`
- histone marks (`H3K*`, `H4K*`, … via `_is_histone_mark`)

Transcription reads a centered sigmoid of activity (`activity_readout`) before sampling RNA production.

CLI: `--activity-sigmoid-m-max` (default 20), `--activity-sigmoid-k` (default 1).

## ChromHMM signed mark weights (SPN v8, opt-in)

When `--chromhmm-emission-tsv` (or `--chromhmm-mark-weights-tsv`) is set with
`layer_rules=v2`, pipeline version becomes `spn_v8_chromhmm_emw`.

For **hist → `*.promoter.active` REGULATORY** only:

\[
w_j=\mathrm{mean}_{i\in A}e_{ij}-\mathrm{mean}_{i\in I}e_{ij},\quad
\delta_{\mathrm{signed}}=\delta\cdot w_j
\]

Sign of \(\delta_{\mathrm{signed}}\) selects activation vs inhibition. Weights are
**global per mark type** (not locus-specific). See `docs/SPN_V8_PIPELINE.md`.

## Ledger QC

Simulation summary / trajectory JSON include `ledger_totals`:

- `seed_inject`, `prod_transcription`, `prod_translation`
- `sink_remove` / `tokens_destroyed`, `transfer_moved`, `regulatory_delta`
- `catalyst_gate` (enzyme/kinase edges that fired as presence-only)
- `regulatory_skipped` (REGULATORY onto non-activity place — should stay ~0)
- `tf_proteolysis`, `illegal_create` (should stay 0 under v2 clamp)

## TF recycle vs UPS timer

Live-cell FRAP / single-molecule work shows most TFs exchange rapidly with a large free pool. Promoter-coupled proteasomal turnover exists for some activators but is not the default for every TF–DNA event. v2 therefore defaults to non-consuming regulatory edges; ε-proteolysis is opt-in.

## Stoichiometry

- `one_to_one`: metabolic edges use conserving transfer (phase B).
- `kegg`: BioPAX KEGGR coefficients (phase C; coverage reported in QC when wired).

## Code map

- `gin_tool/subgraph_pipeline/spn_layer_rules.py` — mode resolution
- `gin_tool/subgraph_pipeline/spn_simulator_gin.py` — firing + ledger
- `gin_tool/subgraph_pipeline/spn_degradation.py` — kind-wise sinks
- `gin_tool/subgraph_pipeline/spn_policy.py` — `SpnSimPolicy.layer_rules`

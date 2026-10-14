# p53 发表数据管线（冻结）· GINv3.11a

**版本 ID**：`p53_pub_gin11a_hs0_20260929`（**当前**；supersedes `p53_pub_gin11a_20260929`）  
**状态**：**FROZEN_PUBLICATION**  
**确定性**：`PYTHONHASHSEED=0` · `context_n_jobs=20` · `--no-debias`  
**锁文件（git）**：`docs/locks/p53_pub_gin11a_hs0_20260929.json`  
**仿真策略**：仍用 **`spn_std_h1_b0`**（见 `docs/SPN_P53_STANDARD_PIPELINE.md`）。

> **原则**：改 GIN / Stage A / H1 / Fig3E 边表前必须升版本 ID；**禁止**原地覆盖冻结路径。重跑须带相同 `PYTHONHASHSEED` 与 `n_jobs`。

---

## 1. 冻结产物（唯一发表路径）

| 阶段 | 路径 | 规模（锁定） |
|------|------|----------------|
| **GIN** | `GINv3.11a.sifi` | 由 `GINv3.11.sifi` 去掉 StageA−hist 多余边；**不**覆盖 `GINv3.11.sifi` |
| **A** 共享 context · rh1 · 无 GRN prune | `gin_out/p53_shared_context_v4_rh1_gin11a/` | **11372** 边 · **4630** 节点 · neg=**332** |
| **B** H1 GRN prune + SPN 导出 | `gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/` | **11367** 边 · SPN **1499** places |
| **C** 标准 B0 脉冲（展示种子） | `runtime/spn_phys5/p53_B0_std_stageA_gin11a_s48_h40_*` | **seed=48** · ΔTP53≈+31.8 · ΔMDM2_post≈+5.04 · follow≈99.7% |

**Context 版本串**：`context_v4_hme_soft_metenz3a_dogma_hme_seedhist_dropnsact_histprune`（`--context-pipeline-version v4`，`--context-mech-depth 2` ≡ real-hop 1）。

**种子**：`bash scripts/build_p53_pi3k_comma_seeds.sh` + `seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt`。

---

## 2. 与旧产物的关系（勿混用）

| 名称 | 路径 | 角色 |
|------|------|------|
| **本发表版** | `…_gin11a` / `GINv3.11a.sifi` | **唯一** Methods/Fig4 数据主线 |
| 旧 Stage A（2026-09-26） | `gin_out/p53_shared_context_v4_rh1`（12318） | 归档对照；**非**发表 |
| 历史复原子图 | `gin_out/p53_hist_repro_subgraph`（11374） | 20260815 lineage 锁；gin11a 与之 Jaccard≈**0.991** |
| 旧 Fig4 盘 | `gin_out/p53_cell_grn_prune_v4_seedhist/h1_wa01` | 历史标定盘；仿真语义仍 `spn_std_h1_b0` |

GINv3.11a 动机：去掉相对 hist_repro 的 Stage-A-only 多余边（主要为未测 TF→seed.promoter 的 HME soft-seed 边），使 Stage A 与历史复原子图对齐。

---

## 3. 复现入口（Slurm）

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 0 构建 GIN.a | `scripts/sbatch_build_ginv3_11a.sh` | 仅当需从 `GINv3.11.sifi` 重建；产物已冻结则跳过 |
| 1 Stage A | `scripts/sbatch_p53_stageA_gin11a.sh` | **勿**默认重跑覆盖；需新版本时改 `OUT=` |
| 2–3 H1 + 脉冲 | `scripts/sbatch_p53_stageA_gin11a_h1_std_pulse.sh` | 全量重建 H1 后默认 seed=42；发表展示用下行 |
| **发表脉冲（推荐）** | `scripts/sbatch_p53_pub_gin11a_h1_pulse.sh` | **只读**冻结 H1；`--seed 48`；跑前 `assert` 锁指纹 |

校验（不改盘）：

```bash
python3 scripts/assert_p53_pub_gin11a_lock.py
```

---

## 4. 发表展示脉冲（seed=48）读数

- **诱导臂**（主证据）：MDM2 脉冲后相对 control 跟随（ΔMDM2_post≈**+5.04**，follow≈**99.7%**）
- **抑制臂**：hold 窗 ΔTP53≈**+31.8**
- **峰时序**（step，t>50）：promoter **198** → RNA **199** → protein **204**（lag_pr=**-1**，lag_prot=**5**）
- **禁止**用 lag median / induction_strict / NF pass-rate 当诱导主结论

图：`runtime/spn_phys5/p53_B0_std_stageA_gin11a_s48_h40_track_means_peaks_gt50_wechat.{png,pdf}`

---

## 5. SHA256 指纹（摘要）

完整表见 `docs/locks/p53_pub_gin11a_20260929.json`。

| 文件 | sha256（前 16 hex） |
|------|---------------------|
| `GINv3.11a.sifi` | `192b2cb5a8147324…` |
| Stage A `subgraph_edges.tsv` | `37cfa1c926500ca9…` |
| H1 `subgraph_edges.tsv` | `be8a73329ed76629…` |

磁盘侧另写：

- `gin_out/p53_shared_context_v4_rh1_gin11a/PUBLICATION_LOCK.json`
- `gin_out/p53_stageA_gin11a_h1_grn_prune/PUBLICATION_LOCK.json`

（`gin_out/` 在 `.gitignore` 内；以 git 锁文件为准做 CI/断言。）

---

## 6. Figure 4C/D（六系相关热图 · **noprune** 冻结）

**版本 ID**：`p53_pub_fig4cd_noprune_gin11a_20261001`  
**状态**：**FROZEN_PUBLICATION** · 文章 **Fig. 4C（Simulated histones）/ 4D（Simulated promoters）** 标准管线  
**锁文件**：`docs/locks/p53_pub_fig4cd_noprune_gin11a_20261001.json`  
**备份**：`gin_out/locks/p53_pub_fig4cd_noprune_gin11a_20261001/`（含 `*_bundle.tar.gz`）

| 项 | 锁定内容 |
|----|----------|
| 输入 | Stage A `gin_out/p53_shared_context_v4_rh1_gin11a`（父锁 hs0；**无 GRN prune**） |
| Levels | 六系统一 `{cell}_rna_plus_met_sources_place_levels.tsv` |
| SPN | `spn_v5_degrade_anc_source` · 200 runs × 200 blocks · `maintain_seeds=True` |
| 工作目录 | `gin_out/p53_fig4_gin11a_noprune_cd/` |
| 图 | `manuscript/figures/fig4/fig4_{C,D}_gin11a_noprune_*.{pdf,png}` |
| 复现 | `sbatch scripts/sbatch_fig4cd_noprune_gin11a_stageA.sh`（**勿**默认覆盖冻结盘） |

校验：

```bash
python3 scripts/assert_p53_pub_fig4cd_noprune_lock.py
```

**明确排除**：prune+emw / `spn_std_h1_b0` 六系盘、旧 `corr_by_layer_phys6_*` **不是** Fig4 C/D 发表标准。风格参照 `fig4-2.pdf` C/D；色条按数据 offdiag nice-floor。

**Fig. 4A 标注数字**（与本冻结同源）：`docs/experiments/fig4a_stats_noprune_gin11a_frozen.md`  
→ places **1497** · transitions **3791** · 200×200 · layer n：protein **328** / histone **512** / metabolite **239** / promoter **81**。

---

## 7. 文档索引

- 本冻结说明：`docs/P53_PUBLICATION_GIN11A_PIPELINE.md`
- Fig4 C/D 锁：`docs/locks/p53_pub_fig4cd_noprune_gin11a_20261001.json`
- 管线审计手递：`docs/P53_PUBLICATION_PIPELINE_HANDOFF.md`（已指向本版本）
- SPN 数值标准：`docs/SPN_P53_STANDARD_PIPELINE.md`
- 实验笔记：`docs/experiments/p53_stageA_gin11a_*`

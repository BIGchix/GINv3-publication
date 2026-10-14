# p53 发表冻结管线 Handoff（笔记本 Cursor 交接 · 2026-10-04）

**用途**：在另一台机器 / 另一个 Cursor 会话中，**不依赖聊天历史**即可理解当前 **FROZEN_PUBLICATION** 主线、两臂仿真、读出规则、稿件草稿与禁止混用路径。  
**项目根**：`/home/chix/data/multiOmicsIntegration`（或笔记本上对应 clone；`gin_out/` 常在 `.gitignore`，大文件可能需单独同步）。

**权威锁文件（git 内，以 SHA 为准）**

| Lock ID | 文件 | 管什么 |
|---------|------|--------|
| `p53_pub_gin11a_hs0_20260929` | [`docs/locks/p53_pub_gin11a_hs0_20260929.json`](locks/p53_pub_gin11a_hs0_20260929.json) | GIN.a · Stage A · Fig3E · H1 prune · B0 脉冲 |
| `p53_pub_fig4cd_noprune_gin11a_20261001` | [`docs/locks/p53_pub_fig4cd_noprune_gin11a_20261001.json`](locks/p53_pub_fig4cd_noprune_gin11a_20261001.json) | Fig. 4A–D 六系 noprune |

**配套说明（勿当作唯一真相；冲突时以锁 JSON 为准）**

- [`docs/P53_PUBLICATION_GIN11A_PIPELINE.md`](P53_PUBLICATION_GIN11A_PIPELINE.md)
- 旧长审计手递（含历史路径）：[`docs/P53_PUBLICATION_PIPELINE_HANDOFF.md`](P53_PUBLICATION_PIPELINE_HANDOFF.md)
- SPN 标准旋钮：[`docs/SPN_P53_STANDARD_PIPELINE.md`](SPN_P53_STANDARD_PIPELINE.md)
- 层规则：[`docs/SPN_LAYER_RULES.md`](SPN_LAYER_RULES.md)
- chromHMM v8：[`docs/SPN_V8_PIPELINE.md`](SPN_V8_PIPELINE.md)

---

## 0. 30 秒心智模型

```text
GINv3.11a.sifi
    → context v4 · real-hop 1 · 无 GRN prune
         = Stage A  (共享边表 · Fig3E rh1)
         │
         ├─ Arm A (Fig 4A–D)  不 prune · 不 chromHMM · spn_v5 · maintain_seeds=True · 六系 RNA+metsrc · 200×200
         │
         └─ Arm B (Fig 4E–I)  H1-only allowlist prune (−5 边) · chromHMM EMW · spn_std_h1_b0 · 350×350 · B0 脉冲
```

**诱导臂判定（硬规则）**：只看脉冲后 **MDM2 相对 control 是否跟随**。禁止用 lag 符号 / induction_strict / NF pass-rate 当主结论。

**改冻结产物前**：升 lock `version_id`；**禁止**原地覆盖锁内路径。

---

## 1. 冻结产物与规模（以 hs0 + fig4cd 锁为准）

| 阶段 | 路径 | 锁定规模 / 指纹要点 |
|------|------|---------------------|
| **GIN** | `GINv3.11a.sifi` | sha256 `192b2cb5a8147324…`；由 `GINv3.11.sifi` 去掉 StageA−hist 多余边；**不覆盖** v3.11 |
| **Stage A** | `gin_out/p53_shared_context_v4_rh1_gin11a/` | 边 **11368** · 节点 **4626** · neg **334** · edges sha256 `350b7cb9…` |
| **Fig3E** | `gin_out/p53_fig3e_v4_gin11a_realhop{1,2,3}_gh{2,4,6}/` | rh1=Stage A；rh2: 17015 / 5097 / 50324；rh3: 24415 / 7461 / 83243 |
| **H1 prune** | `gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/` | 11368→**11363**；SPN **1497** places / **4384** transitions / 12192 arcs |
| **Arm A 工作盘** | `gin_out/p53_fig4_gin11a_noprune_cd/` | 边表 ≡ Stage A；SPN **1497** / **3791**；corr `corr_by_layer_noprune_gin11a_plateau100.json` |
| **脉冲报告** | `runtime/spn_phys5/p53_B0_std_stageA_gin11a_s43_h40_report.json`（锁内） | seed **43** · onset 150 · hold 40 · 350×350 |

> 注意：部分旧文档写 Stage A「11372 / 4630 / 1499」或脉冲「seed=48」。**当前 git 锁 JSON** 为 **11368 / 4626 / 1497** 与 **seed=43**。展示图可能另有 s48 跑次；写稿/复现以锁为准，或同时核对 `runtime/spn_phys5/` 下实际 tag。

**Context 版本串**：`context_v4_hme_soft_metenz3a_dogma_hme_seedhist_dropnsact_histprune`  
`--context-pipeline-version v4` · `--context-mech-depth 2` ≡ **real-hop 1**（graph-hop = 2 × real-hop）。

**种子**：`bash scripts/build_p53_pi3k_comma_seeds.sh` + `seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt`。

**确定性**：`PYTHONHASHSEED=0` · context / 六系 `n_jobs=20` · Stage A 构建 `--no-debias`。

---

## 2. 两条仿真臂（不可混用）

### Arm A — Fig. 4A–D（六系层相关 · **noprune**）

| 项 | 值 |
|----|-----|
| 拓扑 | Stage A **原样**（无 GRN prune） |
| chromHMM | **关** |
| SPN | `spn_v5_degrade_anc_source` · **1497** places / **3791** transitions |
| Levels | 六系统一 `{cell}_rna_plus_met_sources_place_levels.tsv` |
| maintain | **`maintain_seeds=True`**（RNA + met 源每 block 钳回） |
| 仿真 | 200 runs × 200 blocks · plateau = 末 100 blocks 均值 |
| 相关 | Pearson · linkage = average on `(1−r)` |
| 细胞 | H1, H9, HUVEC, IMR90, NHEK(`keratinocyte`), CD4(`tcell_cd4`) |
| 图 | `manuscript/figures/fig4/fig4_{C,D}_gin11a_noprune_*.{pdf,png}`（风格对齐 `fig4-2.pdf`；色条 data-driven） |
| Fig. 4A 数字 | [`docs/experiments/fig4a_stats_noprune_gin11a_frozen.md`](experiments/fig4a_stats_noprune_gin11a_frozen.md)：protein 328 / histone 512 / met 239 / promoter 81 |
| SI | protein / metabolite / TF.protein 热图：`manuscript/figures/fig4/si/fig4_SI_gin11a_noprune_*` |
| Runner | `scripts/run_fig4cd_noprune_gin11a_stageA.py` · `sbatch_fig4cd_noprune_gin11a_stageA.sh` |
| Assert | `python3 scripts/assert_p53_pub_fig4cd_noprune_lock.py` |

**代谢灌种（非代谢组）**：仅 met↔met **in-degree 0 且 out≥1** 的源；默认 level 多为 10；六系几乎共用，系间差异主要靠 **bulk ENCODE RNA TPM**。代谢层相关接近 1 → **不作**细胞系判别证据。

**明确排除（不是 Fig4 C/D）**：prune+emw 六系盘、`spn_std_h1_b0` 六系、旧 `corr_by_layer_phys6_*`、`p53_cell_grn_prune_v4_seedhist`。

### Arm B — Fig. 4E–I（H1 脉冲验收）

| 项 | 值 |
|----|-----|
| 拓扑 | Stage A → **仅 H1** ChIP allowlist prune |
| Allowlist | `runtime/chip_gin_align/h1_wa01/allowlist_union.tsv`（**不是**六系并集） |
| Prune 对象 | TF.protein → `*.promoter;state`；**−5** 边（11368→11363） |
| chromHMM | **开**：H1 E003 panel · \(w_j=\bar e_A-\bar e_I\) · 仅 hist→`promoter.active` REGULATORY · 6 marks（Roadmap core+K27ac）；me2 等无权则不乘 |
| SPN | `spn_std_h1_b0`（layer_rules v2 + 标准旋钮；EMW → v8 语义）· **1497** / **4384** |
| maintain | `maintain_seeds=False` · **`maintain_metabolite_sources=True`** |
| 脉冲 | MDM2.protein ⊣ TP53.protein × **0.1** · onset **150** · hold **40** · 350×350 |
| 读出 | 抑制：hold 窗 TP53↑ vs ctrl；诱导：**post-pulse MDM2 follow vs ctrl** |
| 锁内数值示例（s43） | ΔTP53_hold≈+32.0 · ΔMDM2_post≈+4.77 · follow≈99.4% |
| Assert | `python3 scripts/assert_p53_pub_gin11a_lock.py` |
| 发表脉冲 sbatch | `scripts/sbatch_p53_pub_gin11a_h1_pulse.sh`（只读冻结 H1） |

---

## 3. Context v4 步骤（共享 · 无 prune）

顺序意图：dogma 骨架 → HME soft-seed + histone 链修复（HME **不**直写 `promoter;state`）→ drop 非 seed `promoter.active/;state` → histprune → LCC + met-enzyme 回收。

Fig. 3E：**无** GRN prune、**无** chromHMM（规模普查）。

---

## 4. 复现入口（Slurm · 本机惯例 partition=`local`）

| 动作 | 脚本 | 注意 |
|------|------|------|
| 断言主锁 | `scripts/assert_p53_pub_gin11a_lock.py` | 只读 |
| 断言 Fig4 C/D | `scripts/assert_p53_pub_fig4cd_noprune_lock.py` | 只读 |
| 建 GINv3.11a | `scripts/sbatch_build_ginv3_11a.sh` | 已冻结则跳过 |
| Stage A | `scripts/sbatch_p53_stageA_gin11a.sh` | **禁止**无版本号覆盖冻结 `OUT=` |
| Fig3E hops | `scripts/sbatch_fig3e_v4_gin11a_context_hops.sh` | |
| Fig4 C/D | `scripts/sbatch_fig4cd_noprune_gin11a_stageA.sh` | 勿默认覆盖冻结盘 |
| 发表脉冲 | `scripts/sbatch_p53_pub_gin11a_h1_pulse.sh` | 只读 H1 |
| H1+脉冲全量 | `scripts/sbatch_p53_stageA_gin11a_h1_std_pulse.sh` | 新目录时改路径 |

长任务必须 **sbatch**；禁止 nohup 重仿真。内存峰值 ≤ 本机 RAM 50%。

---

## 5. 稿件草稿状态（2026-10-04）

| 文件 | 角色 |
|------|------|
| [`manuscript/draft/01_results.md`](../manuscript/draft/01_results.md) | Results（六系 noprune + 脉冲；bulk RNA） |
| [`manuscript/draft/03_discussion_v2_20261003.md`](../manuscript/draft/03_discussion_v2_20261003.md) | Discussion **新稿**（未替换旧 `03_discussion.md`） |
| [`manuscript/draft/04_methods_v2_20261004.md`](../manuscript/draft/04_methods_v2_20261004.md) | Methods **新稿**：共享底座 → Arm A → Arm B；公式/策略；无结果数字 |
| [`manuscript/draft/02_abstract_intro.md`](../manuscript/draft/02_abstract_intro.md) | Abstract/Intro（可能仍写旧 1501/4397 —— 需与冻结对齐） |
| [`manuscript/draft/06_si_cross_cell.md`](../manuscript/draft/06_si_cross_cell.md) | SI 跨细胞 |

图权威：Fig4 盘面曾用 `fig4-1001.pdf` / `fig4-2.pdf` 作风格参考；发表 C/D 面板以 `fig4_*_gin11a_noprune_*` 为准。

---

## 6. 禁止混用 / 常见陷阱

| 陷阱 | 正确做法 |
|------|----------|
| 用旧 `p53_shared_context_v4_rh1`（12318） | 用 **`…_gin11a`** Stage A |
| 用 `p53_cell_grn_prune_v4_seedhist` 当 Fig4 C/D | C/D = **noprune** Stage A |
| 六系也做 GRN prune / chromHMM | 仅 **Arm B H1** |
| `nogrnnprune` 并集当 Stage A | 那是 prune 后并集反推，**不是**发表 Stage A |
| Fig3E KEGG-only 旧 hop 表 | 用 **`fig3e_v4_gin11a_*`** |
| lag / induction_strict 当诱导 | **禁止**；诱导 = MDM2 post-pulse follow |
| Methods 写 Δ/r 结果 | Results 写效果；Methods 写 how |
| 本仓库有完整 GIN 从零构建脚本 | **无**；仅有 `build_ginv3_11a_from_stageA_extras.py` 等修补；全库构建在 scBionet 侧 |

---

## 7. 笔记本同步清单

**Git 内（应能 pull）**

- `docs/locks/p53_pub_gin11a_hs0_20260929.json`
- `docs/locks/p53_pub_fig4cd_noprune_gin11a_20261001.json`
- 本文件 + `docs/P53_PUBLICATION_GIN11A_PIPELINE.md`
- `manuscript/draft/0{1,3_discussion_v2_*,4_methods_v2_*}.md`
- `manuscript/figures/fig3/fig3_*gin11a*` · `manuscript/figures/fig4/*gin11a*`
- `scripts/assert_p53_pub_*.py` · `scripts/sbatch_*gin11a*` · `scripts/run_fig4cd_noprune_gin11a_stageA.py`
- `benchmarks/subgraph-speed/fig3e_v4_gin11a_context_hop_stats.json`（若已提交）

**常需单独同步（体积大 / 常被 ignore）**

- `GINv3.11a.sifi`（~825 MB）
- `gin_out/p53_shared_context_v4_rh1_gin11a/`
- `gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/`
- `gin_out/p53_fig4_gin11a_noprune_cd/`（轨迹 JSON 很大）
- `gin_out/p53_fig3e_v4_gin11a_realhop*/`
- `runtime/omics_seed/cellline_levels/{六系}/`
- `runtime/chip_gin_align/h1_wa01/`
- `runtime/chromhmm/panel_emissions/`
- `runtime/spn_phys5/p53_B0_std_stageA_gin11a_*`

到笔记本后先跑两个 `assert_*_lock.py`；缺大文件时 assert 会失败——属预期。

---

## 8. 给下一个 Agent 的操作约定

1. **先读本文件 + 两个 lock JSON**，再改稿或重跑。  
2. 长任务 **sbatch**；改冻结路径必须 **新 version_id**。  
3. 微信/结论句：诱导只报 follow；跨细胞标 descriptive。  
4. Discussion/Methods 新稿在 `*_v2_*.md`；合并进正式 `03_discussion.md` / `04_methods.md` 前需作者确认。  
5. Intro/旧 merged 稿里的 **1501/4397** 是历史数，应用 **1497 / 3791（A）** 与 **1497 / 4384（B）**。

---

## 9. 一句话版本戳

**发表主线 = GINv3.11a + Stage A gin11a（11368）+ Fig4 noprune 六系（spn_v5, 200×200）+ H1-only prune（11363）+ chromHMM EMW + B0 pulse（spn_std_h1_b0, 350×350, lock seed 43）；诱导 = post-pulse MDM2 follow vs control。**

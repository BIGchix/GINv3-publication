# p53标准管线经验总结（`spn_std_h1_b0`）· SigPN

> **命名更正**：本文所述引擎为**SigPN**（signaling Petri net），**不是**stochastic Petri net。  
> 文件名保留历史`SPN_`前缀；正文中的“SPN”均应理解为SigPN。  
> 自定义仿真与Pulse特例说明见 [`SIGPN_CUSTOM_SIMULATIONS.md`](SIGPN_CUSTOM_SIMULATIONS.md)；ChIP allowlist见 [`CHIP_ALLOWLIST.md`](CHIP_ALLOWLIST.md)。

锁定日期：2026-08-23  
来源：Job **382**（`v15cfea10e3ewe_o150b350_metsrc`）及v12→v18标定扫参；代码见`gin_tool/subgraph_pipeline/spn_standard_policy.py`。

本文档总结**子图构建 → 分层SigPN语义 → 标准参数 → B0脉冲协议 → 读数规则**，作为向WNT等通路泛化的基准。原则：**全局规则、可复用语义**，而非通路特化补丁。

---

## 1. 范围与目标

| 项 | 说明 |
|----|------|
| **通路** | p53 / MDM2 负反馈轴（主校准案例） |
| **子图** | `gin_out/p53_cell_grn_prune_v4_seedhist/{cell}` |
| **仿真版本** | `spn_std_h1_b0`（`SPN_PIPELINE_VERSION_STD`） |
| **细胞系** | 六系 phys panel：`h1_wa01`, `h9_wa09`, `huvec`, `imr90`, `keratinocyte`, `tcell_cd4` |
| **核心问题** | 脉冲干预后，**抑制臂**（TP53 抬高）与 **诱导臂**（MDM2 相对 control 跟随）是否可重复 |

**H1 标准结果（hold=40, runs=350）**：ΔTP53 hold **+31** · ΔMDM2 post **+4.42** · follow **99.1%**。

---

## 2. 子图构建（context v4 + seedhist prune）

### 2.1 管线版本

- **Context 提取**：`context_v4_hme_soft_metenz3a_dogma_hme_seedhist_dropnsact_histprune`（`gin_tool/pipeline_version.py`）
- **顺序（生物学意图）**：
  1. **Dogma 先行**：保证 seed 基因的 promoter/RNA/protein 骨架完整
  2. **HME soft-seed + 组蛋白链修复**：`hme_chain_repair` — HME 经 histone mark 间接调控 promoter，**禁止** HME 直接当 TF 写 `promoter;state`
  3. **非 seed 的 `promoter.active` / `promoter;state` 丢弃**，保留 RNA/protein 入口
  4. **Histprune**：去掉无 in-edge 的孤立 histone、无 promoter 连接的 mark 等
  5. **LCC + 代谢酶回收**（post-topology）：连通性修复，代谢超边酶不作 0-hop soft-seed（防 BFS 爆炸）

### 2.2 p53 GRN prune（`p53_cell_grn_prune_v4_seedhist`）

在 context 子图基础上，按 **ChIP/allowlist** 修剪 TF→`promoter;state`：

- 有测量证据的 TF→靶 promoter **保留**
- 无证据且非 spared 列表 → **删除**（减少虚高 GRN 边）
- 典型规模（H1）：`subgraph_edges.tsv` **~11.4k** 边、**~4.6k** GIN 节点；导出 SPN 后 **1501 places**、**10563 interactions**

### 2.3 子图构建时的优化原则（非仿真补丁）

| 优化 | 模块 | 生物学理由 |
|------|------|------------|
| HME 不直连 `;state` | `promoter_edge_policy` + `hme_chain_repair` | 染色质酶通过 mark 间接作用，不应与序列特异性 TF 混用 |
| 非 seed promoter 状态丢弃 | context v4 | 外围基因无 chromatin 支撑时不应驱动转录 |
| Histone 链修复 + eraser 边 | `hme_chain_repair` | 去甲基化/去乙酰化中间产物与 HME 催化链完整 |
| GRN allowlist prune | `p53_cell_grn_prune` | TF-DNA 边需表观/结合证据，抑制全连接虚边 |
| 最大连通子图 | `connectivity.py` | 仿真在单一连通组件上进行 |
| 代谢酶 3-hop 回收 | context v4 | 子图内代谢通量可闭合，又不无限扩张 |

### 2.4 Omics 种子与 chromHMM

- **Place levels**：`runtime/omics_seed/cellline_levels/{cell}/{cell}_rna_plus_met_sources_place_levels.tsv`  
  - RNA 丰度 + **代谢 in-degree-0 源**（`met_source_seed.py`）
- **chromHMM emission**：`runtime/chromhmm/panel_emissions/{cell}_{EID}_promoter_emission.tsv`  
  - promoter 区 active vs inactive 发射差 → **有符号** mark 权重（`spn_v8_chromhmm_emw`）
  - 仅作用于 **histone mark → `promoter.active`**，不消耗 histone

---

## 3. 标准仿真策略一览（`standard_h1_b0_spn_policy`）

### 3.1 数值参数表

| 类别 | 参数 | 标准值 | 说明 |
|------|------|--------|------|
| **版本** | `version` | `spn_std_h1_b0` | 指纹见 `is_standard_h1_b0_knobs()` |
| **时程** | `blocks` / `runs` | 350 / 350（panel 确认 1000） | 平台需 ≥150 block |
| **层规则** | `layer_rules` | `v2` | 见 §4 |
| **代谢** | `stoich` | `one_to_one`（v2 默认） | 代谢边守恒转移 |
| **边权** | `edge_balance` | `metabolite` | 代谢层归一化 |
| **Dogma** | `prod_alpha_transcription` | **1.0** | 转录产出系数 |
| | `prod_alpha_translation` | 1.0 | |
| | `prod_m_max` | **0（无硬顶）** | ≤0 表示 uncapped，用满 promoter/RNA 水平 |
| | `translation_output_scale` | 1.0 | 与 `protein_production_output_scale` 乘积作用于翻译 |
| | `protein_production_output_scale` | 1.0 | 全局 PRODUCE_TRANSLATION 缩放 |
| **时标 fast×e** | `protein_half_life_blocks` | 10/e ≈ **3.68** | 相对 v15c 基准 ÷e |
| | `rna_half_life_blocks` | 5/e ≈ **1.84** | |
| | `promoter_active_half_life_blocks` | 5/e ≈ **1.84** | |
| | `out_route_half_life_blocks` | 20/e ≈ **7.36** | 出边路由稀释 |
| **Promoter** | `promoter_incoming_edge_scale` | **0.1** | 所有入边 → `*.promoter;state` / `*.promoter.active` ×0.1 |
| | `histone_to_promoter_state_scale` | 0.1（再叠加入边 scale → hist 路有效 ×0.01） | 组蛋白→;state |
| | `drop_hme_to_promoter_state` | **True** | 删除 HME→;state |
| | `sign_repressive_chromatin_to_state` | **True** | 甲基化轴抑制酶 → ;state 负权 |
| **TF→DNA** | `tf_dna_route_independent` | **True** | TF 调控独立于 out-route 扇出 |
| | `tf_dna_edge_weight` | **0.5**（× prom_in 后有效 **0.05**） | δ ∝ sample(m_TF × w) |
| **活性更新** | `activity_update_mode` | `relu_active` | promoter.active / ;state |
| | `activity_relu_drive_scale` | **0.2** | |
| | `activity_relu_scale_mode` | **`net`** | 每 block 后缩放净变化（含降解） |
| | `activity_sigmoid_m_max` | 40 | 容量上界 |
| **UPS / E3** | `ups_degradation_policy` | **True** | 见 §5 |
| | `strong_e3_inhibition_weight` | **e** | 强 ESI 抑制弧 ×√e |
| | `weak_e3_edge_scale` | **1/e** | 弱 ESI 抑制弧 ×1/e |
| | `scale_non_e3_activity_inhibitors` | **True** | 非 E3 活性抑制酶 ×0.1（基因表） |
| **维持** | `maintain_seeds` | **False** | 不每 block 重注入全 seed |
| | `maintain_metabolite_sources` | **True** | met↔met **in-degree-0** 源每 block 钳位 |
| **降解** | `add_degradation` | True | 分层 half-life 出口 |
| **种子** | `mark_mode` | `ancestor_all` | in-degree-0 祖先注入 |
| **脉冲** | 见 §6 | onset=150, hold=40 | MDM2⊣TP53 |

### 3.2 刻意不采用的“补丁”

| 废弃方向 | 原因 |
|----------|------|
| `prod_m_max=50` 硬顶（v12） | 卡死 RNA/蛋白，MDM2 蛋白无法跟随 |
| 仅去 cap 不压 promoter（v14） | runaway，QC 失败 |
| 同时 `α_tx×0.1` + `tl×0.1`（v15） | 蛋白轴塌缩，脉冲无效 |
| 用 lag / strict / NF pass-rate 作诱导主结论 | 双向反馈下易误导（见 §7） |

---

## 4. layer_rules=v2：反应类型与处理方式

定义：`gin_tool/subgraph_pipeline/spn_layer_rules.py` — **按 place/边类型自动分类**，非 p53 专用。

| FireMode | 典型边 | 行为 |
|----------|--------|------|
| **REGULATORY** | TF→`promoter;state` / histone→`promoter.active` | 饱和更新（sigmoid / relu）；**不消耗** TF/histone；chromHMM 符号权重 |
| **PRODUCE_TRANSCRIPTION** | `promoter.active`→RNA | ΔRNA = α·sample(level)；**不消耗** promoter.active |
| **PRODUCE_TRANSLATION** | RNA→protein | Δprot = α·sample(RNA)×global_prod_scale；**不消耗** RNA |
| **CATALYST** | 激酶/抑制酶⊣蛋白 | 激活：门控；抑制：只扣靶蛋白（MDM2⊣TP53 属此类） |
| **TRANSFER_CONSERVING** | 代谢、部分复合物 | 源→靶守恒转移 |
| **DEGRADE** | 自降解层 | 一级衰变到 sink；强 ESI 底物可跳过全局蛋白降解 |
| **LEGACY_CONSUME** | layer_rules=off | 旧 gin_consume_v1（已弃用于标准管线） |

**Out-route**（`out_route_half_life_blocks`）：REGULATORY/CATALYST 开火后，源蛋白的多个出边按 fan-out 竞争稀释 — 除 **TF→DNA 独立路由**外。

**TF→DNA 独立**（`tf_dna_route_independent=True`）：每条 TF.protein→promoter 边独立开火，δ ~ sample(m_TF × w)，不参与 out-route 分摊。

---

## 5. Promoter / 染色质边策略

模块：`promoter_edge_policy.py`（导出前作用于 `gin_edges`）

1. **删除** `HME.protein → *.promoter;state`（默认基因表 `nodeList_HME`）
2. **缩放** `histone_mark → *.promoter;state` × `histone_to_state_scale`（标准 0.1）
3. **缩放** 所有入边 → `*.promoter;state` / `*.promoter.active` × `promoter_incoming_edge_scale`（标准 0.1）  
   - H1 p53 子图约 **2628** 条入边被缩放  
   - TF→DNA 有效权重：0.5 × 0.1 = **0.05**
4. **符号校正**：`repressive_chromatin_enzymes.tsv` 中酶 → ;state 强制负权

---

## 6. UPS 降解与抑制强度

模块：`ups_degradation_policy.py` + 表 `seedLists/ups/*.tsv`

| 证据强度 | 处理 |
|----------|------|
| **强 ESI**（如 MDM2–TP53） | E3→virtual→底物 抑制弧 × **√strong_e3**（标准 √e）；底物 **跳过** 全局蛋白自降解 |
| **弱 ESI** | 抑制弧 × **weak_e3**（标准 1/e） |
| **非 E3 活性抑制酶** | 基因列表规则，抑制弧 ×0.1 |

**脉冲可作用边**：`InteractionWeightPulse` 在 block 窗内对匹配 `(source, target, kind)` 的弧权再乘 `low_factor`（标准 B0：`MDM2.protein` ⊣ `TP53.protein`，×0.1）。

---

## 7. B0 脉冲实验设计

### 7.1 时间轴

```
block:  0 -------- 150 -------- 190 -------- 350
        [  平台期 / 预平衡  ] [ hold=40 ] [ 恢复 + 观察 post ]
                              ↑ onset
```

- **onset=150**：控制组在 ~150 block 进入平台后再脉冲（此前 onset=50 时层级耦合差）
- **hold=40**：Nutlin 样抑制削弱窗
- **restore=step**：hold 结束后抑制权重恢复基线
- **ramp_blocks=20**：可选 ramp 恢复（标准用 step）

### 7.2 双臂与读数

| 臂 | 条件 | 主读数 |
|----|------|--------|
| **control** | 无 `interaction_weight_pulse` | 基线漂移、平台高度 |
| **step** | `standard_b0_pulse_step()` | 见下 |

| 指标 | 窗口 | 含义 |
|------|------|------|
| `mean_dTP53_hold_vs_ctrl` | onset : onset+hold | **抑制臂**：TP53 是否被释放抬高 |
| `mean_dMDM2_post_vs_ctrl` | onset+hold : onset+hold+40 | **诱导臂**：MDM2 是否跟随 |
| `rate_follow_vs_ctrl` | per-run (dTP53>0 ∧ dMDM2>0) | 跟随比例 |

**禁止**将 lagged Pearson、induction_strict、NF pass-rate 作为诱导臂主结论（双向负反馈下非因果）。

### 7.3 建议保存的轨迹

- **关键蛋白**：TP53, MDM2, CDKN1A, BAX, PTEN（`.protein`）
- **MDM2 层级**：`MDM2.promoter.active`, `MDM2.RNA`, `MDM2.protein`
- **跨细胞系**：control 轨迹用于 Figure 4 类 promoter.active / histone 相关性（plateau 窗 block 50–149）

---

## 8. 代谢维持

`maintain_metabolite_sources=True` 且 `maintain_seeds=False` 时：

- 每 block 将 **met↔met 有向图 in-degree-0** 的代谢物 place 钳回 seed level（`met_source_seed.py`）
- 防止代谢源在长跑中耗尽导致全图塌缩
- **不**重注入 RNA/蛋白 seed（允许转录翻译动态）

额外：`DEFAULT_EXTRA_DEGRADE_PLACES` 为 H₂O、磷酸等 hub 代谢物提供降解出口。

---

## 9. 标定历程摘要（为何是今日标准）

| 版本 | 关键旋钮 | 结果 | 教训 |
|------|----------|------|------|
| v12 | cap=50, α=0.15, 无 prom scale | promoter~950, 蛋白平 | cap 卡 dogma |
| v13 | 去 cap + 极低 α + tl×0.1 | 轴塌缩 | 多旋钮同时打压 |
| v14 | 仅去 cap | runaway | 需 promoter 尺度 + 时标 |
| v15 | prom×0.1 + α/tl 双降 | 蛋白~0.1 | 勿同时砍翻译 |
| v15b | prom×0.1 + v12 dogma | 可用但诱导弱 | 可行中间态 |
| v15c | v15b + 无 cap | 幅度过大 | 需 fast×e + E3 尺度 |
| **std** | fast×e + strong=e + weak=1/e + metsrc + o150 | follow 99% | **当前基准** |

**单旋钮原则**：优先改 **全局语义参数**（promoter 入边、E3 强弱、half-life 时标、代谢源维持），避免单基因边手工改权。

---

## 10. 代码入口

| 用途 | 路径 |
|------|------|
| 策略定义 | `gin_tool/subgraph_pipeline/spn_standard_policy.py` |
| 策略类型 | `gin_tool/subgraph_pipeline/spn_policy.py` |
| 仿真导出 | `gin_tool/subgraph_pipeline/spn_simulate.py` |
| H1 脉冲 | `scripts/run_p53_b0_std_pulse.py` |
| 六系 panel | `scripts/run_p53_b0_std_panel_n1000.py` |
| Slurm | `scripts/sbatch_p53_B0_std_pulse.sh`, `sbatch_p53_B0_std_panel_n1000.sh` |
| 脉冲调度 | `gin_tool/subgraph_pipeline/interaction_weight_schedule.py` |

---

## 11. 向 WNT 泛化的检查清单（预告）

在 **不改全局语义** 的前提下换子图：

1. **子图**：WNT KEGG / context 提取是否含破坏复合体 hub 与 CTNNB1 层级  
2. **参数**：默认复用 `spn_std_h1_b0` 数值；仅当 plateau/QC 失败再扫 onset / prom_in / fast×e  
3. **脉冲边**：优先 **降解抑制**（GSK3 hub ⊣ CTNNB1），格式同 `InteractionWeightPulse`  
4. **读数**：激活臂 = ΔCTNNB1_hold；诱导臂 = ΔCCND1/MYC_post vs control  
5. **问题回流模型**：若 hub 过密或旁路过多 → 回流 **promoter 边策略 / UPS 分类 / 子图 prune 规则**，而非单条边补丁

---

## 12. 参考文献（仓库内）

- `docs/SPN_V7_PIPELINE.md` — dogma、miRNA、代谢源 v7  
- `docs/SPN_V8_PIPELINE.md` — chromHMM signed emission  
- `.cursor/cli-worker-context.md` — 当前 handoff 与 Job 状态  
- `docs/WECHAT_BRIEF.md` — 标定对话流水  

---

*本文档随 `spn_std_h1_b0` 更新；参数变更请同步修改 `spn_standard_policy.py` 与本文件。*

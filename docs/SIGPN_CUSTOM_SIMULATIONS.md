# SigPN仿真：发表Pulse特例与自定义实验指南

> **SigPN** = signaling Petri net（本项目的信号Petrinet仿真）。  
> 请勿与 **stochastic Petri net（随机PN）** 混淆。  
> 仓库中部分路径/文件名仍含历史前缀`spn_`（如`spn_export_qc.json`、`runtime/spn_phys5/`），语义均为SigPN。

发表稿中的**B0 Pulse（Fig.4e–i）是一个高度特化的示例实验**（H1细胞 · MDM2⊣TP53 · 固定onset/hold/读数）。  
本包提供可复用的**导出与仿真引擎**，但**不**提供通用“任意扰动实验”的完整GUI；其它生物学问题需要用户在引擎之上自行编写实验脚本。

---

## 1. 发表Pulse所用子图是如何做成的

整条链路（与锁`p53_pub_gin11a_hs0_20260929`一致）：

```text
GINv3.11a.sifi
    │  gin_tool.cli  context v4 · mech-depth 2 · no-debias
    ▼
共享 Stage A 边表
  gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv
    │  细胞特异 ChIP allowlist prune（仅滤 TF→promoter;state）
    │  scripts/build_realhop1_cell_spn.py
    ▼
H1 细胞子图
  gin_out/p53_stageA_gin11a_h1_grn_prune/h1_wa01/subgraph_edges.tsv
  （发表：11368 → 11363 边）
    │  导出 bipartite → SigPN 模型 + 标准政策 spn_std_h1_b0
    │  加载 H1 place levels + chromHMM emission
    │  scripts/run_p53_b0_std_pulse.py
    ▼
B0 Pulse 轨迹 / 报告
  runtime/spn_phys5/p53_B0_std_stageA_gin11a_s43_* 
```

### 1.1 Stage A（共享context）

- 代码：`python -m gin_tool.cli`（或`scripts/reproduce_publication.py --steps stage_a`）  
- 要点：默认context/`v4`即可（不必改mode/version）· `--context-mech-depth 2` · `PYTHONHASHSEED=0` · 发表用`n_jobs=20`  
- 产物：机制邻域边表，**尚未**按细胞ChIP修剪  

详见`docs/USER_GUIDE.md`、`REPRODUCE.md`。

### 1.2 H1 GRN allowlist prune → 细胞子图

- 代码：`scripts/build_realhop1_cell_spn.py`  
- 核心库：`gin_tool/subgraph_pipeline/grn_allowlist_prune.py`  
- 输入：  
  - Stage A（或任意context）的`subgraph_edges.tsv`  
  - `runtime/chip_gin_align/h1_wa01/allowlist_union.tsv` + `tf_gate_tfs.txt`  
- 输出：`…/h1_wa01/subgraph_edges.tsv`、`prune_qc.json`；默认还可做SigPN**导出**（`--runs 0`时不做长MC）  

Allowlist含义与自制方法：`docs/CHIP_ALLOWLIST.md`。

### 1.3 B0 Pulse（发表特例）

- 代码：`scripts/run_p53_b0_std_pulse.py`  
- 政策：`gin_tool/subgraph_pipeline/spn_standard_policy.py`（`spn_std_h1_b0`）  
- 仿真入口：`gin_tool/subgraph_pipeline/spn_simulate.py` → `export_and_simulate_subgraph_dir`  
- 发表设定（摘要）：  
  - 细胞子图：H1 prune后目录  
  - 干预：MDM2 ⊣ TP53（抑制边权重×0.1量级的标准B0步骤）  
  - onset=150 · hold=40 · blocks=350 · runs=350 · RNG seed=**43**  
  - place levels：`runtime/omics_seed/cellline_levels/h1_wa01/…tsv`  
  - chromHMM：`runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv`  
  - 读数：ΔTP53_hold、ΔMDM2_post、follow 等（见锁JSON）  

**这不是通用脉冲框架**：换基因、换干预类型、换读数窗口，都需要新脚本或改该脚本的参数/逻辑。

标准政策与分层规则的长文（历史文件名仍含SPN，正文请按SigPN理解）：  
`docs/SPN_P53_STANDARD_PIPELINE.md`、`docs/SPN_LAYER_RULES.md`。

---

## 2. 一键复现发表Pulse（计算，非复制）

```bash
export PYTHONHASHSEED=0
python scripts/reproduce_publication.py --steps prepare
python scripts/reproduce_publication.py --steps stage_a --n-jobs 20
python scripts/reproduce_publication.py --steps prune_h1,pulse --n-jobs 20 --seed 43
```

或分步调用`build_realhop1_cell_spn.py`与`run_p53_b0_std_pulse.py`（见脚本`--help`）。

---

## 3. 可复用 vs 需用户自写

### 3.1 建议直接复用（本包已提供）

|模块 / 脚本|作用|
|---|---|
|`python -m gin_tool.cli`|任意种子的context子图（默认）|
|`grn_allowlist_prune.py` + `build_realhop1_cell_spn.py`|按细胞allowlist修剪TF→`;state`，导出细胞边表 / SigPN骨架|
|`spn_simulate.export_and_simulate_subgraph_dir`|边表→SigPN模型→Monte-Carlo轨迹|
|`spn_standard_policy` / `spn_policy` / `spn_layer_rules`|层语义、半衰期、边权平衡等**可调政策**|
|`runtime/omics_seed/…place_levels.tsv`|六系示例初始marking（可换成你的细胞TSV，格式对齐）|
|`runtime/chromhmm/panel_emissions/…`|示例emission权重（可换成你的细胞）|
|`load_place_levels_tsv`等辅助|读入place水平|

### 3.2 需要用户自行编写的部分

|需求|说明|
|---|---|
|**新的扰动协议**|例如敲低其它节点、边权重时间表、多脉冲、药物组合——发表的`run_p53_b0_std_pulse.py`只实现MDM2⊣TP53 B0|
|**新的读数 / 统计**|锁文件中的ΔTP53/ΔMDM2/follow是该特例的报告；其它表型需自己从`*_trajectory.json` / `pair_runs.npz`汇总|
|**新细胞allowlist**|按`docs/CHIP_ALLOWLIST.md`制作；本包只带H1成品|
|**新通路Stage A以外的实验设计**|种子选择、深度、是否HME/dogma，由你用CLI组合；发表脚本只锁死p53配方|
|**大规模可视化 / 图稿**|可用`plot_fig4_*.py`作参考，但不是通用画图库|

### 3.3 自定义仿真的最小模板（概念）

```text
1. CLI 得到 context 边表
2. （可选）allowlist prune → 细胞边表目录 cell_dir/
3. 准备 place_levels.tsv（及可选 chromHMM tsv）
4. 调用 export_and_simulate_subgraph_dir(
      cell_dir,
      place_levels=...,
      # 以及你的 runs/blocks/seed、政策 kwargs
      # 在政策或后处理中编码「你的」干预
    )
5. 解析 trajectory，写你的指标
```

具体kwargs见`standard_h1_b0_simulate_kwargs()`与`export_and_simulate_subgraph_dir`签名；发表pulse脚本是最完整的**特例参考实现**。

---

## 4. Fig.4a–d（六系noprune）与Pulse的区别

|项|Fig.4a–d noprune|Fig.4e–i Pulse|
|---|---|---|
|子图|共享Stage A，**不做**ChIP prune|Stage A → **H1 allowlist prune**|
|脚本|`scripts/run_fig4cd_noprune_gin11a_stageA.py`|`run_p53_b0_std_pulse.py`|
|目的|六系相关结构 / 层相关|H1上MDM2–TP53脉冲动力学|
|是否特例|是（固定六系与政策）|是（固定干预与读数）|

二者都建立在同一SigPN引擎上；换科学问题仍需自写实验层。

---

## 5. 相关文档索引

|文档|内容|
|---|---|
|`docs/USER_GUIDE.md`|工具包CLI总览|
|`docs/CHIP_ALLOWLIST.md`|allowlist格式与制作|
|`docs/SPN_P53_STANDARD_PIPELINE.md`|标准政策与层语义（按SigPN阅读）|
|`docs/SPN_LAYER_RULES.md`|分层更新规则|
|`REPRODUCE.md`|发表计算复现步骤|

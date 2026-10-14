# 细胞特异ChIP allowlist：用途、格式与制作方法

本文说明下游**细胞特异GRN修剪**所用的ChIP allowlist：它是什么、文件长什么样、如何制作、如何接到本工具包。

> 命名：本项目的仿真引擎称为**SigPN**（signaling Petri net / 信号Petrinet），**不是**stochastic Petri net（随机PN）。文档一律写SigPN；代码路径里可能仍出现历史目录名`spn_*`（如`runtime/spn_phys5/`、`spn_export_qc.json`），语义上均指SigPN。

---

## 1. 用途（为什么需要）

Stage A等**共享context子图**在GIN上是通路/机制邻域，**不带细胞特异ChIP过滤**。  
进入细胞系仿真前，通常要对边表做一次**GRN allowlist prune**，只动一类边：

`TF.protein → GENE.promoter;state`（TF ∈ `nodeList`中的TF集合）

规则（实现：`gin_tool/subgraph_pipeline/grn_allowlist_prune.py`）：

| 情况 | 行为 |
|------|------|
| 该TF在本细胞**未测到**（不在`tf_gate_tfs.txt`） | **保留**其在子图中的全部TF→`;state`边（图谱不全≠该TF沉默） |
| 该TF在本细胞**已测到** | **仅当**`(TF, gene)`出现在`allowlist_union.tsv`时保留 |
| 其它边（dogma、代谢、miRNA、histone链等） | **不动** |
| 修剪后度数为0的节点 | 删除 |

因此allowlist回答的是：在该细胞里，**哪些“已测TF → 靶基因”有ChIP（或同类）结合证据**。

发表H1例子目录：

`runtime/chip_gin_align/h1_wa01/`  
（约64.5万对；`tf_gate_tfs.txt`约99个已测TF）

---

## 2. 目录与文件结构

每个细胞一个子目录（`<cell_id>`可自定，如`h1_wa01`、`my_hepg2`）：

```text
runtime/chip_gin_align/<cell_id>/
  allowlist_union.tsv   # 必需：TF–基因允许对（可多库合并）
  tf_gate_tfs.txt       # 必需：本细胞“已测TF”符号表（一行一个）
  qc.json               # 可选：构建过程摘要（本包H1带示例）
  peaks/                # 可选：原始/中间peak BED（本发布包通常不附带高峰文件）
```

接入修剪脚本时通过：

```bash
--chip-root runtime/chip_gin_align --cell <cell_id>
```

解析为`runtime/chip_gin_align/<cell_id>/allowlist_union.tsv`与`tf_gate_tfs.txt`。

---

## 3. 文件格式示例

### 3.1 `allowlist_union.tsv`（制表符分隔，含表头）

|列|含义|
|---|---|
|`tf`|TF基因符号（须能映射到`nodeList` TF / `TF.protein`）|
|`gene`|靶基因符号（对应子图中`GENE.promoter;state`）|
|`dbs`|可选：来源库，逗号分隔（如`chip_atlas,gtrd,remap`）|
|`n_dbs`|可选：库个数|

示例（节选，真实H1文件极大）：

```tsv
tf	gene	dbs	n_dbs
ASH2L	AAAS	chip_atlas,gtrd,remap	3
ASH2L	AACS	chip_atlas,remap	2
CTCF	CDKN1A	chip_atlas,gtrd	2
TP53	MDM2	chip_atlas	1
```

最小合法文件：只要表头含`tf`与`gene`两列即可；加载时只使用这两列做集合`(tf, gene)`。

### 3.2 `tf_gate_tfs.txt`

一行一个HGNC符号，`#`开头为注释：

```text
# Measured TFs in this cell (peak experiments mapped to nodeList_TF)
ASH2L
CTCF
TP53
MDM2
```

**语义**：出现在此文件中的TF视为“本细胞已测”；prune时对它们启用allowlist硬过滤。未出现的TF走“未测→全保留”分支。

---

## 4. 推荐制作流程（用户侧）

本发布包**不捆绑**完整peak下载与对齐流水线（体积与许可证原因）。用户可按下列步骤自建，或使用实验室内`gin_tool.chip_gin_align`构建器（若你持有完整源码树）。

### 步骤A — 准备启动子窗口

- 用人基因组TSS±窗口（发表构建常用约TSS −1000 / +100 bp量级）得到`promoters.bed`：`chrom start end gene_symbol …`
- 基因符号尽量与GIN/`nodeList`一致（HGNC）。

### 步骤B — 获取该细胞的TF ChIP peaks

任选或合并：

- **ChIP-Atlas**（按细胞/组织的AllAg或实验BED）
- **GTRD** / **ReMap** 等细胞匹配的peak集

将TF实验标签**映射为基因符号**（去掉抗体别名、epitope tag等无法映射的条目）。

### 步骤C — Peak ∩ 启动子 → TF–基因对

对每个数据库：

1. `bedtools intersect`（或等价）将peaks与`promoters.bed`求交  
2. 得到原始`(tf_label, gene)`对  
3. 过滤无法映射的TF/基因  
4. 写出该库的pair表

### 步骤D — 合并为`allowlist_union.tsv`

- 对`(tf, gene)`取并集  
- 可选填写`dbs` / `n_dbs`  
- 写出带表头的TSV

### 步骤E — 写`tf_gate_tfs.txt`

- 取“本细胞至少一个peak实验成功映射”的TF符号集合  
- 建议与`nodeList`中TF池求交，避免孤儿符号

### 步骤F — 接入本包修剪

```bash
# 1) 已有共享context边表（Stage A或自定义context）
# 2) 已准备 runtime/chip_gin_align/my_cell/{allowlist_union.tsv,tf_gate_tfs.txt}

python scripts/build_realhop1_cell_spn.py \
  --context-edges gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv \
  --context-dir   gin_out/p53_shared_context_v4_rh1_gin11a \
  --out-root      gin_out/my_cell_grn_prune \
  --cell          my_cell \
  --chip-root     runtime/chip_gin_align \
  --node-list-dir ./nodeList \
  --runs 0
```

输出：`gin_out/my_cell_grn_prune/my_cell/subgraph_edges.tsv` + `prune_qc.json` 等，可作为后续SigPN导出/仿真的细胞子图。

---

## 5. 质量检查建议

- `n_allowlist_pairs`、已测TF数是否合理（过少→几乎不剪；过多且噪声大→过剪）  
- 对照`prune_qc.json`中`n_tf_to_state_removed` / `n_tf_to_state_spared_unmeasured`  
- 抽查通路关键边（如TP53→靶）是否仍在修剪后子图中  
- TF/基因符号与`nodeList`不一致是最常见失败原因

---

## 6. 与发表H1产物的关系

|项|发表H1|
|---|---|
|目录|`runtime/chip_gin_align/h1_wa01/`|
|用途|Fig.4e–i H1 prune + B0 pulse的细胞边过滤|
|注意|这是**H1 allowlist**，不是六系并集；其它细胞需各自目录|

手写/自建allowlist时，**不必**复现H1的64万对；按你的细胞与数据库重新生成即可，只要格式与prune规则一致。

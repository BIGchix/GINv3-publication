# 用户指南：用GINv3-publication做自己的分析

本包是**通用GIN多组学工具包**；p53文稿管线只是其上的一套配方。

> **命名**：仿真引擎写作**SigPN**（signaling Petri net）。不要与stochastic Petri net（随机PN）混淆。代码/目录里历史名`spn_*`均指SigPN。

---

## 1. 安装与解压先验

```bash
pip install -e ".[sim]"
export PYTHONHASHSEED=0
python scripts/reproduce_publication.py --steps prepare
# → ./GINv3.11a.sifi（约787MB），校验sha256 192b2cb5…
```

也可将`--sifi`指向其它兼容的`GINv*.sifi`。

---

## 2. 核心CLI（任意种子）

```bash
python -m gin_tool.cli ...
# 等价：gin-multomics ...
python -m gin_tool.cli -h
```

**参数说明**：[CLI_REFERENCE.md](CLI_REFERENCE.md)（子图模式与pipeline版本用默认即可，不必改）。

**必填**：`--seeds`（逗号分隔基因符号或完整GIN节点名）。  
**常用**：`--sifi`、`--out-dir`、`--node-list-dir ./nodeList`。

### 2.1 Context子图

```bash
python -m gin_tool.cli \
  --sifi ./GINv3.11a.sifi \
  --seeds "TP53,MDM2,ATM,ATR" \
  --out-dir gin_out/demo_tp53 \
  --context-mech-depth 2 \
  --context-n-jobs 8 \
  --node-list-dir ./nodeList \
  --no-debias
```

|标志|含义|
|---|---|
|`--context-mech-depth`|图跳预算（发表约定：2 ≡ real-hop 1）|
|`--context-n-jobs`|并行度（发表冻结用**20**）|
|`--active-tfs-file`|跳过TF排序，使用给定列表|
|`--no-dogma-chain-repair` / `--no-hme-chain-repair`|关闭相应修复|

产物含`subgraph_edges.tsv`、`run_meta.json`、各类repair QC等。

示例脚本：`examples/01_custom_context.sh`。

### 2.2 种子文件

```bash
SEEDS=$(bash scripts/build_p53_pi3k_comma_seeds.sh)   # 发表用p53+代谢物
# 或：SEEDS=$(grep -v '^#' my_genes.txt | paste -sd, -)
```

---

## 3. 下游：细胞特异修剪与SigPN

共享context之后，常见两步：

1. **ChIP allowlist prune** → 细胞特异边表（滤`TF→promoter;state`）  
2. **SigPN导出 / Monte-Carlo**（政策、初值、干预由你设定）

### 3.1 ChIP allowlist

完整说明（格式、制作步骤、与prune规则）：**[`docs/CHIP_ALLOWLIST.md`](CHIP_ALLOWLIST.md)**。

接入示例：

```bash
python scripts/build_realhop1_cell_spn.py \
  --context-edges gin_out/<your_context>/subgraph_edges.tsv \
  --context-dir   gin_out/<your_context> \
  --out-root      gin_out/<your>_cell_prune \
  --cell          <cell_id> \
  --chip-root     runtime/chip_gin_align \
  --node-list-dir ./nodeList \
  --runs 0
```

要求事先存在：

`runtime/chip_gin_align/<cell_id>/{allowlist_union.tsv,tf_gate_tfs.txt}`

本包预置发表用H1：`runtime/chip_gin_align/h1_wa01/`。其它细胞需自建。

### 3.2 SigPN仿真与发表Pulse特例

发表的H1 B0 Pulse（MDM2⊣TP53）是**特例**，不是通用实验框架。

请阅读：**[`docs/SIGPN_CUSTOM_SIMULATIONS.md`](SIGPN_CUSTOM_SIMULATIONS.md)**

其中说明：

- Pulse子图如何由Stage A → allowlist prune → 细胞边表得到  
- 对应脚本与政策代码路径  
- **哪些模块可复用**，**哪些干预/读数代码需用户自写**

标准政策长文（文件名历史含SPN，请按SigPN理解）：  
`docs/SPN_P53_STANDARD_PIPELINE.md`、`docs/SPN_LAYER_RULES.md`。

---

## 4. 发表配方 vs 自己的分析

|目标|入口|
|---|---|
|自己的种子 / 深度 / 输出|`python -m gin_tool.cli`（本文§2）|
|自建细胞allowlist并修剪|`docs/CHIP_ALLOWLIST.md` + `build_realhop1_cell_spn.py`|
|自定义SigPN干预实验|复用`spn_simulate`引擎；实验脚本自写（见`SIGPN_CUSTOM_SIMULATIONS.md`）|
|重算文稿Stage A / Fig3e / Fig4 / Pulse|`scripts/reproduce_publication.py --steps …`（**真实计算**）|
|与文稿锁比对|`scripts/assert_p53_pub_*.py`|
|仅检查发布包指纹（非复现）|`scripts/validate_bundle.py` / `--from-expected`|

---

## 5. 确定性提示

- 设置`PYTHONHASHSEED=0`  
- 比较边表哈希时固定`--context-n-jobs`  
- 与发表Stage A **bit一致**在原Linux冻结主机上已验证；其它平台可能HME桥接路径略有差异——见`docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md`

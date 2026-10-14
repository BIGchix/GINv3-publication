# `gin_tool.cli` 参数说明

入口（等价）：

```bash
python -m gin_tool.cli -h
gin-multomics -h
```

下文只说明**推荐使用的选项**。CLI里其它历史开关请勿依赖；本包文档不展开。发表配方见`REPRODUCE.md`。

子图提取保持程序默认即可（context · pipeline `v4`）；命令行不必写出`--subgraph-mode` / `--context-pipeline-version`。

---

## 1. 输入 / 输出

|选项|默认|说明|
|---|---|---|
|`--sifi`|`GINv3.11a.sifi`|GIN边表路径|
|`--seeds`|（必填）|逗号分隔基因符号或GIN节点名|
|`--out-dir`|`gin_out`|输出目录|
|`--node-list-dir`|仓库`nodeList/`|分层列表（自动k、debias、TF/HME池）|
|`--max-lines`|无|只读前N条边（测试用）|
|`--symbol-blacklist` / `--symbol-aliases`|仓库config下默认路径|映射前过滤/别名|

---

## 2. Context深度与并行

|选项|默认|说明|
|---|---|---|
|`--context-mech-depth`|`2`|有向BFS的**图跳**上限。严格R↔V二分图上：`2` ≡ **real-hop 1**（发表Stage A）|
|`--context-max-nodes`|`100000`|扩展节点上限|
|`--context-n-jobs`|`1`|并行度；发表锁用**20**。`0`/负数→自动取CPU|
|`--no-promoter-reverse-walk`|关（即默认开启反向）|从种子promoter沿入边回溯TF|
|`--no-lcc-bridge`|关（默认做LCC桥）|岛桥接|
|`--lcc-bridge-hops`|`2`|LCC桥的图跳预算|

---

## 3. Active TF

|选项|说明|
|---|---|
|`--active-tfs`|逗号分隔，跳过排序|
|`--active-tfs-file`|每行一个；发表用`seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt`|
|`--active-tf-method`|`rwr_cohesion`（默认）/ `rwr_rewire` / `fisher`|
|`--skip-active-tf`|完全跳过|
|`--context-include-active-tf-anchors`|遗留：把active TF也当BFS锚点（默认只要seeds）|

---

## 4. Dogma / HME repair

默认在context下**开启**；用`--no-*`关闭。

|选项|作用|
|---|---|
|`--no-dogma-chain-repair`|不补`*.RNA`/`*.promoter`链|
|`--no-hme-chain-repair`|不做组蛋白/HME桥接与修复|
|`--hme-max-real-hops`|PPI→HME桥real-hop上限（默认3）|
|`--hme-upstream-graph-hops`|从hist mark回溯HME的图跳（默认4）|
|`--no-eraser-chain-repair`|关闭eraser催化链补边|
|`--eraser-max-rounds`|去甲基/去乙酰阶梯轮数（默认3）|
|`--no-hme-active-hist-backfill`|关闭对seed `promoter.active`的同基因hist回填|
|`--hme-backfill-all-actives`|遗留：对所有active回填（默认仅seed）|
|`--no-drop-nonseed-promoter-active`|保留非seed的active/`;state`|
|`--no-hme-sourceless-hist-prune`|关闭无入边hist/孤岛清理|

---

## 5. RWR / debias / 分析

|选项|说明|
|---|---|
|`--alpha`|RWR重启概率（默认0.3）|
|`--debias` / `--no-debias`|运行时null校正（发表Stage A用**`--no-debias`**）|
|`--null-tail-alpha` / `--null-debias-mode` / `--null-runtime-replicates`|null细节|
|`--analyze`|Louvain模块 + 跨层流 + schematic|
|`--viz`及相关|可视化（需viz可选依赖）|

---

## 6. 发表Stage A常用命令对照

```bash
export PYTHONHASHSEED=0
SEEDS=$(bash scripts/build_p53_pi3k_comma_seeds.sh)

python -m gin_tool.cli \
  --sifi ./GINv3.11a.sifi \
  --seeds "$SEEDS" \
  --out-dir gin_out/p53_shared_context_v4_rh1_gin11a \
  --context-mech-depth 2 \
  --context-n-jobs 20 \
  --node-list-dir ./nodeList \
  --active-tfs-file seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt \
  --no-debias
```

---

## 7. 其它文档

|文档|内容|
|---|---|
|`docs/USER_GUIDE.md`|工具包使用流程|
|`docs/CHIP_ALLOWLIST.md`|细胞ChIP allowlist|
|`docs/SIGPN_CUSTOM_SIMULATIONS.md`|SigPN / Pulse特例|
|`REPRODUCE.md`|文稿逐步重算|

# Stage A复现对照记录（服务器 vs 本机Mac / Docker）

**文档用途**：留档说明发表锁Stage A拓扑的来源、服务器两次独立运行、本机（含Docker）重算结果，以及二者差异，供审稿或外部复现询问时引用。  
**记录日期**：2026-10-14  
**锁版本**：`p53_pub_gin11a_hs0_20260929`  
**权威锁文件**：`docs/locks/p53_pub_gin11a_hs0_20260929.json`  
**Stage A边表指纹（发表）**：sha256 `350b7cb9acc8eaa7b1454f3c912f5aa1d94742a5f05b7492ed01b447039e6a43` · **11368**边 · **4626**节点 · neg **334**

> **一句话结论**：同一套输入与CLI参数下，服务器Linux环境两次独立重算均bit对齐发表锁；本机macOS（arm64）与Docker Linux aarch64从零重算得到另一稳定指纹`5ac7f9d5…`（Jaccard≈0.9836）。发表数字与下游H1/pulse应使用**冻结边表**，而非本机重算边表。分叉点在HME chain repair的directed real-hop路径选择。

---

## 1. 共同输入与CLI契约（服务器与本机对齐部分）

以下在服务器备份树与本机复现目录中已逐项核对为一致（内容sha相同或字符串相同）。

|项目|取值|
|---|---|
|GIN先验|`GINv3.11a.sifi` · sha256 `192b2cb5a8147324fc41064a3f20c943a3d82c03972cff544e3d26219e7ace83` · 824892563字节|
|种子|`bash scripts/build_p53_pi3k_comma_seeds.sh`（85个逗号分隔symbol，两端字符串identical）|
|Active TFs|`seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt`|
|nodeList|仓库根`nodeList/`（33个文件，两端逐文件sha一致）|
|确定性|`PYTHONHASHSEED=0`|
|模式|`--subgraph-mode context` · `--context-pipeline-version v4` · `--context-mech-depth 2`（≡ real-hop 1）|
|并行|发表锁要求`context_n_jobs=20`（勿与`n_jobs=1`混用后声称同一锁）|
|Debias|`--no-debias`|
|Context版本串|`context_v4_hme_soft_metenz3a_dogma_hme_seedhist_dropnsact_histprune`|

标准CLI（与`scripts/sbatch_p53_pub_gin11a_hs0_refreeze.sh` / `…_repro.sh`一致）：

```bash
export PYTHONHASHSEED=0
SEEDS="$(bash scripts/build_p53_pi3k_comma_seeds.sh)"
python3 -u -m gin_tool.cli \
  --sifi GINv3.11a.sifi \
  --seeds "$SEEDS" \
  --out-dir <OUT> \
  --subgraph-mode context \
  --context-pipeline-version v4 \
  --context-mech-depth 2 \
  --context-n-jobs 20 \
  --active-tfs-file seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt \
  --no-debias
```

服务器脚本通常**不**显式传`--node-list-dir`（默认解析为工作目录下`nodeList/`）。`run_meta.json`记录为：

`/disk1/usrData/chixData/multiOmicsIntegration/nodeList`

本机显式传`--node-list-dir ./nodeList`时，文件内容与服务器nodeList一致。

### 1.1 直至HME之前的QC（服务器冻结 vs 本机stock重算）

下列计数在冻结`run_meta`与本机mac_rebuild`run_meta`中一致，说明分叉不在种子扩展/dogma/LCC等前序步骤：

|QC项|冻结 / 本机stock|
|---|---|
|context prune|88859 → 2335节点；anchors=331|
|role soft seeds|primary=231 · hme_soft=100|
|met-enzyme|+enz=196 · +edges=688|
|dogma-repair|+nodes=1362 · +edges=12245 · repaired=318|
|LCC bridge|cc 43→11 · bridge_edges=68|
|active TF|method=`user_supplied` · count=61|
|HME outside / bridged|18 / 18（两端相同）|

**唯一系统性分叉**：`hme_chain_repair`终态边/节点增量（见第4节）。

---

## 2. 服务器上获得Stage A的两次运行

工作根目录（freeze/repro脚本）：`/disk1/usrData/chixData/multiOmicsIntegration`  
（早期单跑脚本`sbatch_p53_stageA_gin11a.sh`写`cd /home/chix/data/multiOmicsIntegration`，与上为同项目树的路径别名。）

备份镜像：`/Volumes/Xu_work/data/multiomicsIntegration`（含冻结产物、repro产物、锁JSON与Slurm痕迹）。

### 2.1 运行A — 发表冻结（Slurm job **1190**）

|项|内容|
|---|---|
|脚本|`scripts/sbatch_p53_pub_gin11a_hs0_refreeze.sh`|
|作业名|`pub_gin11a_hs0` · job id **1190**|
|资源|#SBATCH `--cpus-per-task=20` · `--mem=128G` · `--partition=local` · `--time=24:00:00`|
|开始|锁记录`started`: **2026-09-29T15:04:05+08:00**；`locked_at`: 2026-09-29T07:21:47Z|
|环境变量|`PYTHONHASHSEED=0`（脚本`assert`强制）|
|`context_n_jobs`|20（=`SLURM_CPUS_PER_TASK`）|
|Stage A输出|`gin_out/p53_shared_context_v4_rh1_gin11a/`（**覆盖写入**规范发表路径）|
|同作业后续|Fig3E hops、H1 prune、B0 pulse，并写锁JSON|

**Stage A结果（权威）**

|指标|值|
|---|---|
|边数|11368|
|节点数|4626|
|neg|334|
|edges sha256|`350b7cb9acc8eaa7b1454f3c912f5aa1d94742a5f05b7492ed01b447039e6a43`|
|`hme_chain_repair.n_nodes_added`|1522|
|`hme_chain_repair.n_edges_added`|-7735|
|`hme_chain_repair.n_reprune_nodes_removed`|1219|

旁证文件：`gin_out/p53_shared_context_v4_rh1_gin11a/PUBLICATION_LOCK.json`、`docs/locks/p53_pub_gin11a_hs0_20260929.json`。  
Slurm err残留：`benchmarks/subgraph-speed/p53_pub_gin11a_hs0_1190.err`（含`/usr/bin/time` maxrss约3.4–3.5GB量级片段；完整stdout未随备份保留）。

### 2.2 运行B — 盲复现（Slurm job **1194**）

|项|内容|
|---|---|
|脚本|`scripts/sbatch_p53_pub_gin11a_hs0_repro.sh`|
|作业名|`pub_hs0_repro` · job id **1194**|
|资源|#SBATCH `--cpus-per-task=20` · `--mem=128G` · `--partition=local` · `--time=20:00:00`|
|时间窗|**2026-09-29T16:20:22+08:00** → **16:31:35+08:00**|
|环境变量|`PYTHONHASHSEED=0`|
|`context_n_jobs`|20|
|Stage A输出|`gin_out/p53_shared_context_v4_rh1_gin11a_repro/`（**不**覆盖冻结目录）|
|目的|相对锁JSON做拓扑/sha/脉冲数值盲比对|

**Stage A结果**

|指标|值|
|---|---|
|edges sha256|`350b7cb9…`（与锁**完全一致**）|
|边数|11368|
|与冻结边表Jaccard|**1.0**（only_frozen=0 · only_repro=0）|
|HME终态|与冻结相同（`n_nodes_added=1522`等）|
|整管线|`REPRO_OK: true`（Stage A + H1 + pulse ΔTP53/ΔMDM2/follow均exact）|

记录：`docs/experiments/p53_pub_gin11a_hs0_repro.md`、`docs/experiments/p53_pub_gin11a_hs0_repro_1194.json`。

**解释**：在同一服务器主机与同一代码/数据树上，Stage A发表指纹可在约11分钟级盲复现中bit对齐——说明锁在**该环境**下是可重复的，而非一次性手工编辑边表。

### 2.3 服务器环境要点（由备份推断）

|维度|证据|
|---|---|
|OS/调度|Linux + Slurm；工作路径`/disk1/usrData/chixData/...`|
|Python|备份中`gin_tool/__pycache__`以**cpython-312**为主（另有少量3.8痕迹）；作业调用`python3`|
|声明依赖|`numpy>=1.24` · `scipy>=1.11` · `networkx>=3.2`（无精确version pin留档）|
|架构|集群Linux x86_64类主机（与本机arm64不同；备份未另存`uname -a`）|

---

## 3. 本机Mac与Docker重算

复现工作目录：`github/ginv3-repro-stageA-mac/`（及trial下`…_mac_rebuild`）。  
主机：macOS 15.7.4 · **arm64** · 物理内存约24GB。  
本机venv：Python **3.9.6** · numpy 2.0.2 · scipy 1.13.1 · networkx 3.2.1。

### 3.1 主结果：stock代码从零重算（应与审稿对照的“本机重算”）

|运行|平台|`n_jobs`|edges sha256（前16）|边数|vs锁Jaccard|only_new / only_lock|
|---|---|---|---|---|---|---|
|Mac原生stock|macOS arm64|20（及trial中`n_jobs=6`的mac_rebuild）|`5ac7f9d5bc3a09bf…`|11390|0.983614|105 / 83|
|Docker Linux aarch64 stock|`python:3.9-slim` · numpy/scipy manylinux aarch64|6（内存所限；前序QC与冻结一致）|`5ac7f9d5bc3a09bf…`|11390|0.983614|105 / 83|

要点：

1. Mac与**aarch64 Linux Docker**得到**同一**边表指纹`5ac7f9d5…`，故差异**不能**简单归因于“macOS文件系统/大小写”，而是**本机所用CPU架构/软件栈**与服务器冻结主机不一致时，HME路径first-hit不同。  
2. `n_jobs=6`与`20`在本机stock路径上未改变该错误指纹（相对锁仍为`5ac7f9d5…`一类结果）；发表锁仍规定服务器侧使用20。  
3. HME终态相对锁：`n_nodes_added` 1524 vs 1522；`n_edges_added` -7713 vs -7735；`n_reprune_nodes_removed` 1215 vs 1219。

典型漂移机制（已在诊断中确认）：directed real-hop把含`;`的`GENE.promoter;state`当作virtual免费层，本机first-hit常走`*→HDAC4.promoter;state→…→miR→HME`捷径；冻结拓扑更多保留`ALDH4A1`/`TCF7L1`等桥接。详见`github/ginv3-repro-trial/logs/HOW_TO_MATCH_PUBLICATION_LOCK.md`。

### 3.2 诊断性实验（非发表路径；说明曾尝试算法对齐）

|实验|结果sha（前16）|vs锁Jaccard|备注|
|---|---|---|---|
|HME路径tie-break（排序sources/neighbors + `promoter;state`字典序）|`3dc05920…`|≈0.9967|PRDM9/ALDH4A1改善；仍残留ABL1/APP/BRCA1 vs冻结TCF7L1/ATR等|
|HME real-hop中`*.promoter;state`不计免费virtual|`e5b62c24…`|0.990856（+6/-98）|HDAC4边数可与锁同为9，但ALDH4A1回收不足（29 vs 72）；说明冻结本身仍按`;`virtual语义走桥|

上述修改**均未**达到`350b7cb9…`，也**未**写入发表锁；仅作机制记录。

### 3.3 本机保证发表锁一致的做法

```bash
rsync -a github/raw/gin_out/p53_shared_context_v4_rh1_gin11a/ \
        <repro>/gin_out/p53_shared_context_v4_rh1_gin11a/
# 或使用备份：/Volumes/Xu_work/data/multiomicsIntegration/gin_out/p53_shared_context_v4_rh1_gin11a/
```

安装后本机边表sha=`350b7cb9…`（MATCH）。下游H1/pulse/图应消费该冻结目录。

---

## 4. 服务器 vs 本机：对照总表

|维度|服务器 job1190 / 1194|本机Mac / Docker aarch64 stock|
|---|---|---|
|能否从零bit对齐锁`350b7cb9…`|**能**（1194盲复现Jaccard=1）|**否**（稳定得到`5ac7f9d5…`）|
|边数|11368|11390|
|Jaccard（对本机stock）|1.0（相对自身）|0.983614 vs锁|
|边差|—|多105 / 少83|
|HME `n_nodes_added`|1522|1524|
|前序QC（dogma/met/LCC/role等）|见§1.1|与服务器一致|
|OS / CPU|Linux Slurm · `/disk1/...`|macOS arm64；Docker亦为aarch64|
|Python|备份迹象为3.12为主|3.9.6（venv）/ 3.9.25（Docker）|
|`context_n_jobs`（发表规定）|20|本机主对照含20；Docker对照为6|
|内存申请|128G（freeze/repro）|主机~24G · Docker~16G|
|推荐对外陈述|锁由服务器环境生成并在同机复现|笔记本电脑从零重算不保证bit一致；发表用冻结边表|

### 4.1 若审稿人询问“能否复现Stage A？”

建议答复要点：

1. **发表拓扑**以锁文件`docs/locks/p53_pub_gin11a_hs0_20260929.json`中的edges sha为准；产物在`gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv`。  
2. 在**原计算服务器**上，用同一脚本与`PYTHONHASHSEED=0`、`n_jobs=20`，job1194已证明可**完全**复现（Jaccard=1，sha一致）。  
3. 在**消费级Mac（arm64）**或**aarch64 Linux容器**上，即使用相同sifi/nodeList/seeds与CLI，从零重算会稳定得到另一指纹（`5ac7f9d5…`，Jaccard≈0.984），差异限于HME bridge路径选择，不影响已冻结拓扑作为发表基准的使用。  
4. 跨机器bit级复现发表数字的推荐方式：分发/校验冻结Stage A边表（及锁JSON），再在其上跑H1与pulse；而非要求审稿人笔记本重算Stage A。

---

## 5. 关键路径索引

|角色|路径|
|---|---|
|锁JSON|`docs/locks/p53_pub_gin11a_hs0_20260929.json`|
|冻结Stage A|`gin_out/p53_shared_context_v4_rh1_gin11a/`|
|服务器盲复现Stage A|`gin_out/p53_shared_context_v4_rh1_gin11a_repro/`|
|1194比对JSON|`docs/experiments/p53_pub_gin11a_hs0_repro_1194.json`|
|freeze / repro脚本|`scripts/sbatch_p53_pub_gin11a_hs0_refreeze.sh` · `…_repro.sh`|
|本机Mac重算边表（stock）|`github/ginv3-repro-trial/gin_out/p53_shared_context_v4_rh1_gin11a_mac_rebuild/`|
|本机实验日志|`github/ginv3-repro-stageA-mac/logs/STAGE_A_REPRO_ATTEMPT.txt` · `stage_a_docker_host.log`|
|操作说明（用冻结边表）|`github/ginv3-repro-trial/logs/HOW_TO_MATCH_PUBLICATION_LOCK.md`|
|备份根|`/Volumes/Xu_work/data/multiomicsIntegration`|

---

## 6. 版本与维护

|字段|值|
|---|---|
|本文档ID|`p53_stageA_server_vs_mac_repro_20261014`|
|对应锁|`p53_pub_gin11a_hs0_20260929`|
|若重新冻结Stage A|须升锁版本ID，并更新本文对照表与sha|
|勿混淆|管道总览`docs/P53_PUBLICATION_GIN11A_PIPELINE.md`中部分历史边数（如11372）可能早于hs0终锁；**以锁JSON的11368 / `350b7cb9…`为准**|

# Release notes — v0.2.0

> **状态**：已验收（2026-10-15）。Release **不挂** `GINv3.11a.sifi.gz` 附件；先验从仓库 `data/` 下载（见下文链接）。

**tag**：`v0.2.0`（对齐 `pyproject.toml` / `CHANGELOG.md`）  
**日期**：2026-10-15

---

## 一句话

GINv3-publication首个公开Code Availability版本：可安装的`gin-multiomics`工具包 + 基于**压缩先验**`data/GINv3.11a.sifi.gz`的p53 / GINv3.11a发表复现配方（无需强制解压）。

---

## Highlights

- **默认先验**：`data/GINv3.11a.sifi.gz`（~75 MB）；CLI与`reproduce_publication.py`直接读取
- **工具包**：`pip install -e ".[sim]"` → `gin-multiomics` / `python -m gin_tool.cli`
- **发表配方**：Stage A →（可选）ChIP prune → SigPN；锁指纹见`docs/locks/`
- **命名**：GIN = Global Integrative Network；SigPN = signaling Petri net（非stochastic PN）

---

## GIN prior download（本Release不附资产）

| 用途 | URL |
|------|-----|
| 浏览 / 说明 | https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/data/GINv3.11a.sifi.gz |
| 直接下载（raw） | https://github.com/BIGchix/GINv3-publication/raw/v0.2.0/data/GINv3.11a.sifi.gz |
| 校验 sidecar | https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/data/GINv3.11a.sifi.gz.sha256 |

推荐：`git clone`整仓（先验已在`data/`）。若只取先验文件，下载后放到仓库的`data/GINv3.11a.sifi.gz`再跑`prepare`。

---

## Install

```bash
git clone https://github.com/BIGchix/GINv3-publication.git
cd GINv3-publication
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip && pip install -e ".[sim]"
export PYTHONHASHSEED=0
python scripts/reproduce_publication.py --steps prepare
# verifies data/GINv3.11a.sifi.gz
```

Quick context（默认`--sifi`已是压缩包）：

```bash
python -m gin_tool.cli \
  --seeds "TP53,MDM2,CDKN1A" \
  --out-dir gin_out/my_context \
  --context-mech-depth 2 \
  --context-n-jobs 8 \
  --node-list-dir ./nodeList \
  --no-debias
```

---

## Prior fingerprints

| Artifact | sha256 |
|----------|--------|
| `data/GINv3.11a.sifi.gz` | `f9bb23c44c30d0dd6eab73225600085b378d1912e59f1ff0d8a3a0b37cfd0a58` |
| Uncompressed content | `192b2cb5a8147324fc41064a3f20c943a3d82c03972cff544e3d26219e7ace83` |
| Stage A edges (Linux lock) | `350b7cb9…` · **11,368** edges |

---

## Reproduce manuscript (summary)

```bash
export PYTHONHASHSEED=0
python scripts/reproduce_publication.py --steps prepare
python scripts/reproduce_publication.py --steps stage_a --n-jobs 20
python scripts/assert_p53_pub_gin11a_lock.py
# full: --steps all --n-jobs 20
```

详情：`REPRODUCE.md`。平台差异：`docs/experiments/p53_stageA_server_vs_mac_repro_20261014.md`。

---

## Packaging note

- 本Release对应**源码树**；GIN先验在仓库路径`data/GINv3.11a.sifi.gz`（见上方下载链接），**不**作为Release附件上传。
- `pip`安装的wheel**不含**75MB先验（仅`gin_tool`/`scbionet`）；跑Stage A请clone本仓库或按链接下载先验到`data/`。

---

## GitHub Release正文

## GINv3-publication v0.2.0

First public **Code Availability** release of the GINv3 multi-omics toolkit
(`gin-multiomics`) and the GINv3.11a p53 publication recipe.

### What's included
- Context subgraph CLI + SigPN (signaling Petri net) simulation
- Default prior: `data/GINv3.11a.sifi.gz` (read directly; no unpack required)
- Publication driver: `scripts/reproduce_publication.py`
- Docs: README, USER_GUIDE, REPRODUCE, CHIP_ALLOWLIST, SIGPN_CUSTOM_SIMULATIONS

### Install

    git clone https://github.com/BIGchix/GINv3-publication.git
    cd GINv3-publication
    pip install -e ".[sim]"
    export PYTHONHASHSEED=0
    python scripts/reproduce_publication.py --steps prepare

### GIN prior download

This release does **not** attach the ~75 MB prior as a binary asset. Get it from the repository:

- File: https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/data/GINv3.11a.sifi.gz
- Raw download: https://github.com/BIGchix/GINv3-publication/raw/v0.2.0/data/GINv3.11a.sifi.gz
- Checksums: https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/data/GINv3.11a.sifi.gz.sha256

(`git clone` already includes `data/GINv3.11a.sifi.gz`.)

### Prior checksums
- gzip: `f9bb23c44c30d0dd6eab73225600085b378d1912e59f1ff0d8a3a0b37cfd0a58`
- uncompressed: `192b2cb5a8147324fc41064a3f20c943a3d82c03972cff544e3d26219e7ace83`

### Notes
- Stage A lock (**11,368** edges) is bit-matched on the Linux freeze host; other platforms may differ slightly — see `docs/experiments/`.
- Full changelog: https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/CHANGELOG.md
- Contact: chix@big.ac.cn

---

## 验收清单

- [x] 版本号`0.2.0` / tag名`v0.2.0`同意
- [x] Highlights与Install命令无误
- [x] 先验sha256与`data/GINv3.11a.sifi.gz.sha256`一致
- [x] Packaging note（wheel不含先验）表述可接受
- [x] GitHub Release正文可用（含先验下载链接）
- [x] Release **不挂**`GINv3.11a.sifi.gz`附件；说明文档给出仓库下载地址
- [x] 验收通过 · 2026-10-15

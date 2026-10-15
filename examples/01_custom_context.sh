#!/usr/bin/env bash
# Example: extract a small custom context with the toolkit CLI.
# Uses data/GINv3.11a.sifi.gz by default (no unpack required).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"

if [[ -f data/GINv3.11a.sifi.gz ]]; then
  SIFI=data/GINv3.11a.sifi.gz
elif [[ -f GINv3.11a.sifi ]]; then
  SIFI=./GINv3.11a.sifi
else
  echo "Missing GIN prior — run: python scripts/reproduce_publication.py --steps prepare"
  exit 1
fi

OUT="${OUT:-gin_out/example_custom_tp53}"
SEEDS="${SEEDS:-TP53,MDM2,CDKN1A,ATM,ATR}"
N_JOBS="${N_JOBS:-4}"

echo "[example] seeds=$SEEDS → $OUT (n_jobs=$N_JOBS) sifi=$SIFI"
python -m gin_tool.cli \
  --sifi "$SIFI" \
  --seeds "$SEEDS" \
  --out-dir "$OUT" \
  --context-mech-depth 2 \
  --context-n-jobs "$N_JOBS" \
  --node-list-dir ./nodeList \
  --no-debias

echo "[example] done: $OUT/subgraph_edges.tsv"
wc -l "$OUT/subgraph_edges.tsv"

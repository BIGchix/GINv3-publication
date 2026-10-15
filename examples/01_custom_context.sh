#!/usr/bin/env bash
# Example: extract a small custom context with the toolkit CLI.
# Default --sifi is data/GINv3.11a.sifi.gz (no unpack).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"

if [[ ! -f data/GINv3.11a.sifi.gz ]]; then
  echo "Missing data/GINv3.11a.sifi.gz — run: python scripts/reproduce_publication.py --steps prepare"
  exit 1
fi

OUT="${OUT:-gin_out/example_custom_tp53}"
SEEDS="${SEEDS:-TP53,MDM2,CDKN1A,ATM,ATR}"
N_JOBS="${N_JOBS:-4}"

echo "[example] seeds=$SEEDS → $OUT (n_jobs=$N_JOBS; default sifi=.gz)"
python -m gin_tool.cli \
  --seeds "$SEEDS" \
  --out-dir "$OUT" \
  --context-mech-depth 2 \
  --context-n-jobs "$N_JOBS" \
  --node-list-dir ./nodeList \
  --no-debias

echo "[example] done: $OUT/subgraph_edges.tsv"
wc -l "$OUT/subgraph_edges.tsv"

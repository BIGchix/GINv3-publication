#!/usr/bin/env bash
# Demo: toy subgraph → SigPN export + short Monte-Carlo (laptop-scale).
# Does not require the GIN prior (uses bundled examples/demo_sigpn/).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"
# Prefer this repo's installable package over a leftover PYTHONPATH.
unset PYTHONPATH
export PYTHONPATH="$ROOT"

DEMO_SRC="${DEMO_SRC:-examples/demo_sigpn}"
OUT="${OUT:-gin_out/example_sigpn_demo}"
RUNS="${RUNS:-2}"
BLOCKS="${BLOCKS:-10}"
SEED="${SEED:-0}"

CHROMHMM="${CHROMHMM:-runtime/chromhmm/panel_emissions/h1_wa01_E003_promoter_emission.tsv}"
if [[ ! -f "$CHROMHMM" ]]; then
  echo "Missing chromHMM TSV: $CHROMHMM"
  exit 1
fi
if [[ ! -f "$DEMO_SRC/subgraph_edges.tsv" ]]; then
  echo "Missing demo edges: $DEMO_SRC/subgraph_edges.tsv"
  exit 1
fi

mkdir -p "$OUT"
cp -f "$DEMO_SRC/subgraph_edges.tsv" "$OUT/subgraph_edges.tsv"
cp -f "$DEMO_SRC/place_levels.tsv" "$OUT/place_levels.tsv"

echo "[demo] SigPN export+sim → $OUT (runs=$RUNS blocks=$BLOCKS seed=$SEED)"
python examples/run_sigpn_demo.py \
  --subgraph-dir "$OUT" \
  --place-levels "$OUT/place_levels.tsv" \
  --chromhmm-emission-tsv "$CHROMHMM" \
  --runs "$RUNS" \
  --blocks "$BLOCKS" \
  --seed "$SEED"

echo "[demo] done"
ls -la "$OUT/spn.graphml" "$OUT/spn_export_qc.json" 2>/dev/null || true
ls -la "$OUT"/spn_sim_demo/* 2>/dev/null | head || true

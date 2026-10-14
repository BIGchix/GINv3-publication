#!/usr/bin/env bash
# Example: recompute publication Stage A (same flags as the freeze).
# This RUNS the pipeline — it does not copy expected/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PYTHONHASHSEED=0
N_JOBS="${N_JOBS:-20}"

if [[ ! -f GINv3.11a.sifi ]]; then
  python scripts/reproduce_publication.py --steps prepare
fi

echo "[pub Stage A] n_jobs=$N_JOBS HASHSEED=$PYTHONHASHSEED"
python scripts/reproduce_publication.py --steps stage_a --n-jobs "$N_JOBS"

echo "[pub Stage A] assert against lock (may fail on non-freeze hosts; see docs)"
python scripts/assert_p53_pub_gin11a_lock.py --skip-sifi || true

python - <<'PY'
import hashlib
from pathlib import Path
LOCK = "350b7cb9acc8eaa7b1454f3c912f5aa1d94742a5f05b7492ed01b447039e6a43"
p = Path("gin_out/p53_shared_context_v4_rh1_gin11a/subgraph_edges.tsv")
h = hashlib.sha256(p.read_bytes()).hexdigest()
print("RESULT sha", h[:16], "MATCH" if h == LOCK else "NOMATCH", "lines", sum(1 for _ in p.open()))
PY

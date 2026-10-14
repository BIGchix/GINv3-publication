#!/usr/bin/env bash
# Comma-separated --seeds for gin_tool.cli (KEGG p53 genes + curated metabolites).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
G="$ROOT/seedLists/p53_pi3k_curated/genes.txt"
M="$ROOT/seedLists/p53_pi3k_curated/metabolites.txt"
python3 - "$ROOT/seedLists/p53_pi3k_curated" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
parts: list[str] = []
for name in ("genes.txt", "metabolites.txt"):
    p = root / name
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            parts.append(line)
print(",".join(parts))
PY

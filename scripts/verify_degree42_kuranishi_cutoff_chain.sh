#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON:-$ROOT/.venv/bin/python}"

cd "$ROOT"

"$PYTHON_BIN" scripts/verify_degree42_conormal_rees_synchronization.py
"$PYTHON_BIN" scripts/verify_degree42_divisor_rees_reduction.py
"$PYTHON_BIN" scripts/verify_degree42_kuranishi_branches.py
"$PYTHON_BIN" scripts/verify_degree42_discriminant_quartics.py
"$PYTHON_BIN" scripts/verify_degree42_ab_residual_quartics.py
"$PYTHON_BIN" scripts/verify_degree42_higher_gcd_strata.py

echo "PASS: degree-42 Kuranishi nilpotence cutoff chain"

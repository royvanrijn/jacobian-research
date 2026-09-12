#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
DEFAULT_PYTHON=$(cd "$ROOT/../../.." && pwd)/venv/bin/python
PYTHON=${PYTHON:-$DEFAULT_PYTHON}

if [[ ! -x "$PYTHON" ]]; then
  echo "Python environment not found: $PYTHON" >&2
  exit 2
fi
if [[ ! -f "$ROOT/hard/h_certificate_exact.txt" ]]; then
  echo "Missing hard/h_certificate_exact.txt" >&2
  exit 2
fi

cd "$ROOT"

shasum -a 256 -c EXACT_SHA256SUMS.txt
shasum -a 256 -c RECONSTRUCTED_CERTIFICATES.sha256

"$PYTHON" verify_laurent_reduction.py
"$PYTHON" verify_case1_reduction.py
"$PYTHON" verify_firstblock_exact.py

# Regenerate the exact residual systems consumed by the certificate replay.
"$PYTHON" case2_exact_generate.py
"$PYTHON" case1_cascade_machine.py
"$PYTHON" derive_hne0.py
"$PYTHON" build_degree5.py
"$PYTHON" derive_hne0_branch2.py
"$PYTHON" build_degree5_branch2.py

"$PYTHON" verify_serialized_certificates.py
"$PYTHON" hard/verify_certificate_gmpy2.py
"$PYTHON" verify_hne0_branch_symmetry.py

echo "JC2_72_108_EXACT_REPLAY_PASS"

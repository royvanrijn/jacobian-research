#!/bin/bash
set -euo pipefail
ROOT="/Users/royvanrijn/Documents/jacobian-research"
cd "$ROOT"

PRIME="${1:-100003}"

echo "=== q24/orbit85 modular D12 signature ==="
sage -python elkies-k3/scripts/extract_h92_q24_d12_modp_signature.sage --prime "$PRIME"

echo
echo "=== q24 D12 -> A1 selected-path toolbox check ==="
sage -python elkies-k3/scripts/check_h92_q24_to_a1_toolbox.sage --prime "$PRIME"

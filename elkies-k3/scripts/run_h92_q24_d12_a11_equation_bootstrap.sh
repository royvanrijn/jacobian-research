#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
PRIME="${1:-100003}"

sage -python elkies-k3/scripts/anchor_h92_q24_d12_spinor_modp.sage --prime "$PRIME"
sage -python elkies-k3/scripts/scan_h92_q24_d12_a11_equation_friendly.sage --prime "$PRIME"

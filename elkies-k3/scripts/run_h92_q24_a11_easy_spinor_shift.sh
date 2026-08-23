#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/analyze_h92_q24_a11_easy_spinor_shift.sage --prime "${1:-100003}"

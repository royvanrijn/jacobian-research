#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/test_h92_q24_a11_reuse_f8_bisection.sage --prime "${1:-100003}"

#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/optimize_h92_q24_orbit42_direct_aj_combination.sage --prime "${1:-100003}"

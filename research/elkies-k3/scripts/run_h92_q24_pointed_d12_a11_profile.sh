#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/profile_h92_q24_pointed_d12_a11.sage --prime "${1:-100003}"

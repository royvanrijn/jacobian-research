#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/recover_h92_q24_pointed_zero_pole_sections.sage --prime "${1:-100003}"

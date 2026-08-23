#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/recover_h92_q24_a11_orbit23_p1_section.sage --prime "${1:-100003}"

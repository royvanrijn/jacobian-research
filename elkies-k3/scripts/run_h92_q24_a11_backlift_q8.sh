#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/backlift_h92_q24_a11_to_q8_equation_frame.sage --prime "${1:-100003}"

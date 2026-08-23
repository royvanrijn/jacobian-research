#!/bin/bash
# status: HISTORICAL_DIAGNOSTIC
# Superseded by the fixed q24 D42 resolved-RR construction route.
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/archive/recover_h92_q24_pointed_zero_pole_sections.sage --prime "${1:-100003}"

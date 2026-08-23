#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/audit_h92_q24_explicit_curves_for_a11.sage --prime "${1:-100003}"

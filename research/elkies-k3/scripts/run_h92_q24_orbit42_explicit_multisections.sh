#!/bin/bash
set -euo pipefail
cd /Users/royvanrijn/Documents/jacobian-research
sage -python elkies-k3/scripts/audit_h92_q24_orbit42_explicit_multisections.sage --prime "${1:-100003}"

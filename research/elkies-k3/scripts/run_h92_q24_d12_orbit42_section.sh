#!/bin/bash
set -euo pipefail
ROOT="/Users/royvanrijn/Documents/jacobian-research"
cd "$ROOT"
PRIME="${1:-100003}"

sage -python \
  elkies-k3/scripts/recover_h92_q24_d12_orbit42_section_modp.sage \
  --prime "$PRIME"

#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-/Users/royvanrijn/Documents/jacobian-research}"
P="${2:-100003}"
MAX="${3:-12}"
cd "$ROOT"
sage -python elkies-k3/scripts/scan_h92_q24_d12_local_module_splits_modp.sage --prime "$P" --max-order "$MAX"

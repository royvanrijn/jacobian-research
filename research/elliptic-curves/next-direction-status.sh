#!/usr/bin/env bash
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
research_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
exec python3 "$script_dir/cas/run_next_direction_benchmark.py" "${1:-status}" \
  --folder "$research_dir/artifacts/local/elliptic-curves/next-direction-benchmark-v1"

#!/usr/bin/env bash
set -euo pipefail

mkdir -p artifacts/local/elliptic-curves

echo "=== corrected profiler ==="
sage elliptic-curves/cas/analyze_anonymous_rank_candidate_v2.sage \
  2>&1 | tee artifacts/local/elliptic-curves/anonymous-profile-v2.log

echo
echo "=== initial point search ==="

sage -python elliptic-curves/cas/search_anonymous_initial_points.py \
  --heights 1000,10000,100000 \
  --timeout 30 \
  --e29-centers 12 \
  2>&1 | tee artifacts/local/elliptic-curves/anonymous-initial-point-search.log

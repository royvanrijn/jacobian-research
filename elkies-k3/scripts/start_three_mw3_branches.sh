#!/bin/bash
set -u
mkdir -p artifacts/local/elkies-k3/three-branches

echo "Starting three MW3 reconstruction branches..."

# A10 branch: actual local-Tate modular solve.
(
  PYTHONUNBUFFERED=1 python3 elkies-k3/scripts/run_mw3_local_probe.py     --stage p1 --p 101 --threads 4 --timeout 300     2>&1 | tee artifacts/local/elkies-k3/three-branches/a10.log
) &

# E6 branch: exact frame/MW preparation.
(
  PYTHONUNBUFFERED=1 sage elkies-k3/scripts/analyze_mw3_branch.sage     --frame elkies-k3/data/fibrations/mw3_e6_a3a3_a1a1_frame.txt     --name mw3-e6     2>&1 | tee artifacts/local/elkies-k3/three-branches/e6.log
  sage elkies-k3/scripts/describe_mw3_alt_fiber_model.sage --branch e6     2>&1 | tee -a artifacts/local/elkies-k3/three-branches/e6.log
) &

# A6/A4 branch: exact frame/MW preparation.
(
  PYTHONUNBUFFERED=1 sage elkies-k3/scripts/analyze_mw3_branch.sage     --frame elkies-k3/data/fibrations/mw3_a6_a4_a1x4_frame.txt     --name mw3-a6     2>&1 | tee artifacts/local/elkies-k3/three-branches/a6.log
  sage elkies-k3/scripts/describe_mw3_alt_fiber_model.sage --branch a6     2>&1 | tee -a artifacts/local/elkies-k3/three-branches/a6.log
) &

wait
echo "All three initial branch jobs finished."

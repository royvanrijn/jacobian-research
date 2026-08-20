#!/bin/bash
set -u
mkdir -p artifacts/local/elkies-k3/e6-attack

(
  sage elkies-k3/scripts/enumerate_e6_component_labels.sage     2>&1 | tee artifacts/local/elkies-k3/e6-attack/components.log
) &

(
  sage elkies-k3/scripts/build_e6_mw3_fiber_scaffold.sage     2>&1 | tee artifacts/local/elkies-k3/e6-attack/fibers.log
) &

wait
echo "E6 component and fiber scaffolds finished."

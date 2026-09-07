#!/bin/bash
set -u
mkdir -p artifacts/local/elkies-k3/e6-construction

echo "Inspecting E6 section stages..."
sage elkies-k3/scripts/build_e6_mw3_section_system.sage --stage p1 \
  | tee artifacts/local/elkies-k3/e6-construction/p1-meta.log
sage elkies-k3/scripts/build_e6_mw3_section_system.sage --stage p13 \
  | tee artifacts/local/elkies-k3/e6-construction/p13-meta.log
sage elkies-k3/scripts/build_e6_mw3_section_system.sage --stage all \
  | tee artifacts/local/elkies-k3/e6-construction/all-meta.log

echo "Starting only the P1 modular probe..."
python3 elkies-k3/scripts/run_e6_mw3_probe.py --stage p1 --p 101 --threads 8 --timeout 300 \
  | tee artifacts/local/elkies-k3/e6-construction/p1-solve.log

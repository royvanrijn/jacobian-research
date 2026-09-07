#!/usr/bin/env bash
# Operational supervisor only: does not change the frozen policy or budgets.
set -euo pipefail
cd /home/royvanrijn/src/jacobian-research
v3_run_dir=research/artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3
v3_cas=research/elliptic-curves/cas
v3_sage=/home/royvanrijn/.local/bin/sage
trap 'printf "SUPERVISOR FAILED at line %s; inspect the step logs.\n" "$LINENO" >&2' ERR
printf 'Resuming frozen V3 M17 replay.\n'
"$v3_sage" -python "$v3_cas/adaptive_visibility_cascade_v3.sage" run --start 17 >> "$v3_run_dir/replay-M17-run.log" 2>&1
printf 'Search terminal; independently replaying all stages.\n'
"$v3_sage" -python "$v3_cas/check_visibility_cascade_v3.sage" --start 17 --output "$v3_run_dir/replay-M17.json" --progress "$v3_run_dir/replay-M17-progress.json" >> "$v3_run_dir/replay-M17-check.log" 2>&1
printf 'Recomputing height decision forms.\n'
"$v3_sage" -python "$v3_cas/check_visibility_metric_v3.sage" --start 17 > "$v3_run_dir/metric-M17.log" 2>&1
printf 'Verifying V2 preservation.\n'
"$v3_sage" -python "$v3_cas/adaptive_visibility_cascade_v3.sage" preservation > "$v3_run_dir/preservation-final.log" 2>&1
"$v3_sage" -python "$v3_cas/diagnose_v2_terminal.sage" check-seal > "$v3_run_dir/seal-final.log" 2>&1
printf 'Binding certificates and packaging result.\n'
"$v3_sage" -python "$v3_cas/finalize_visibility_cascade_v3.py" > "$v3_run_dir/finalize.log" 2>&1
printf 'COMPLETE: research/artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v3.json\n'

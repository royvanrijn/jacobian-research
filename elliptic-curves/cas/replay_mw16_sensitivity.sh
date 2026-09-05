#!/usr/bin/env bash
# Explicitly bounded, checkpointed recovery and prospective replay from repo root.
set -euo pipefail
cd -- "$(dirname -- "$0")/../.."
work=artifacts/local/elliptic-curves/mw16-sensitivity
proof=artifacts/generated-results/elliptic-curves
mkdir -p "$work"
sage -python elliptic-curves/cas/calibrate_icarm_mw16_pointed_sieve.sage --height-bound 100000 --timeout-seconds 15 --output "$work/gauss-h100000-controls.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --seed "$work/gauss-h100000-controls.json" --specifications red --heights 10000,100000 --output "$work/red-initial.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --seed "$work/gauss-h100000-controls.json" --specifications 'red:1,2,0;red:2,1,0;red:2,2,1' --heights 100000 --output "$work/rational-initial.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --seed "$work/gauss-h100000-controls.json" --specifications red --heights 200000 --centres specialized,generic --seconds 60 --output "$work/red-centres-h200000.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --seed "$work/gauss-h100000-controls.json" --specifications 'red;metric:1/16;metric:16' --heights 10000,100000 --centres generic --output "$work/generic-metrics.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --seed "$work/gauss-h100000-controls.json" --specifications red --heights 100000 --centres generic --workers 8 --output "$work/red-generic-h100000.json"
python3 elliptic-curves/cas/freeze_mw16_sensitivity_setting.py --input "$work/generic-metrics.json" --output "$work/frozen-initial.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --mode adaptive --curves 398,400 --seed "$work/red-generic-h100000.json" --specifications red --heights 10000 --centres cvp,residue --workers 8 --output "$work/adaptive-red-active.json"
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --mode adaptive --curves 398,400 --seed "$work/frozen-initial.json" --specifications metric:16 --heights 10000 --centres cvp,residue --workers 8 --output "$work/adaptive-metric16-active.json"
python3 elliptic-curves/cas/freeze_mw16_adaptive_policy.py --initial "$work/red-generic-h100000.json" --trials "$work/adaptive-red-active.json" --output "$work/frozen-adaptive-red.json"
python3 elliptic-curves/cas/freeze_mw16_adaptive_policy.py --initial "$work/frozen-initial.json" --trials "$work/adaptive-metric16-active.json" --output "$work/frozen-adaptive.json"
python3 elliptic-curves/cas/verify_mw16_sensitivity.py --campaign "$work/red-initial.json" --campaign "$work/rational-initial.json" --campaign "$work/red-centres-h200000.json" --campaign "$work/frozen-adaptive-red.json" --campaign "$work/frozen-adaptive.json" --bundle "$proof/mw16_sensitivity_controls_v1.json.gz" --summary "$proof/mw16_sensitivity_controls_summary_v1.json"
python3 elliptic-curves/cas/build_mw16_sensitivity_gate.py --initial "$work/frozen-initial.json" --adaptive "$work/frozen-adaptive.json" --verified "$proof/mw16_sensitivity_controls_summary_v1.json" --bundle "$proof/mw16_sensitivity_controls_v1.json.gz" --output "$work/prospective-gate.json"
# Restore the frozen chart-selection ledger from retained evidence if needed.
python3 - <<'PY'
from pathlib import Path
import gzip
source=Path('artifacts/generated-results/elliptic-curves/icarm_mw16_pointed_sieve_h10000_v1.json.gz')
target=Path('artifacts/local/elliptic-curves/mw16-sensitivity/prospective-seed.json')
target.write_bytes(gzip.decompress(source.read_bytes()))
PY
sage -python elliptic-curves/cas/run_mw16_sensitivity.sage --mode prospective --seed "$work/prospective-seed.json" --calibration "$work/prospective-gate.json" --specifications metric:16 --heights 100000 --centres generic --workers 8 --output "$work/prospective.json"
python3 elliptic-curves/cas/verify_mw16_sensitivity.py --campaign "$work/prospective.json" --bundle "$proof/mw16_sensitivity_prospective_v1.json.gz" --summary "$proof/mw16_sensitivity_prospective_summary_v1.json"

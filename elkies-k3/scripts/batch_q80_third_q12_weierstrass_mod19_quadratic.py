#!/usr/bin/env python3
"""Run deterministic mapped Weierstrass samples for p=19 interpolation."""

import argparse
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKER = ROOT / "elkies-k3/scripts/sample_q80_third_q12_weierstrass_mod19_quadratic.sage"
LOCAL = ROOT / "artifacts/local/elkies-k3/q80-third-q12-p19-weierstrass-samples"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-weierstrass-sample-batch.json"
EXPECTED_STATUS = "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--attempts", type=int, default=72)
parser.add_argument("--workers", type=int, default=6)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not 1 <= args.attempts <= 361:
    raise ValueError("attempts must lie between 1 and 361")

LOCAL.mkdir(parents=True, exist_ok=True)
candidates = [(index % 19, index // 19) for index in range(args.attempts)]


def run_one(value):
    a, b = value
    output = LOCAL / f"T-{a:02d}-{b:02d}.json"
    if output.is_file():
        try:
            payload = json.loads(output.read_text())
            if (
                payload.get("status") == EXPECTED_STATUS
                and payload.get("worker", {}).get("sha256") == sha256(WORKER)
            ):
                return {"value": [a, b], "path": output, "cached": True, "error": None}
        except (OSError, json.JSONDecodeError):
            pass
    command = (
        "sage",
        "-python",
        str(WORKER),
        "--new-base",
        str(a),
        "--new-base-r",
        str(b),
        "--output",
        str(output),
    )
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
        check=False,
    )
    if completed.returncode:
        return {
            "value": [a, b],
            "path": output,
            "cached": False,
            "error": completed.stdout[-4000:],
        }
    payload = json.loads(output.read_text())
    if payload.get("status") != EXPECTED_STATUS:
        raise ArithmeticError(f"sample {a}+{b}r has an unexpected status")
    return {"value": [a, b], "path": output, "cached": False, "error": None}


with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    results = list(executor.map(run_one, candidates))

successes = []
failures = []
for record in results:
    if record["error"] is not None:
        failures.append(
            {
                "new_base_coefficients_1_r": record["value"],
                "diagnostic_tail": record["error"],
            }
        )
        continue
    path = record["path"]
    payload = json.loads(path.read_text())
    successes.append(
        {
            "new_base_coefficients_1_r": record["value"],
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "cached": record["cached"],
            "a1_a2_a3_a4_a6": payload["weierstrass"]["a1_a2_a3_a4_a6"],
            "discriminant": payload["weierstrass"]["discriminant"],
            "j": payload["weierstrass"]["j"],
            "inverse_weighted_bounds": [
                payload["birational_maps"]["inverse"]["W"]["weighted_bound"],
                payload["birational_maps"]["inverse"]["old_x"]["weighted_bound"],
            ],
        }
    )
if len(successes) < 56:
    raise ArithmeticError(
        f"only {len(successes)} successful samples from {len(candidates)} attempts"
    )

training_count = len(successes) - 8
output = {
    "schema": "elkies-k3.q80-third-q12-weierstrass-sample-batch-modp2.v1",
    "status": "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"},
    "attempt_count": len(candidates),
    "success_count": len(successes),
    "failure_count": len(failures),
    "training_count": training_count,
    "held_out_count": 8,
    "training_samples": successes[:training_count],
    "held_out_samples": successes[training_count:],
    "failures": failures,
    "worker": {"path": str(WORKER.relative_to(ROOT)), "sha256": sha256(WORKER)},
    "claim_boundary": (
        "Each retained fibre has exact forward/inverse maps to a coherently "
        "Laurent-normalized Weierstrass model. The batch itself does not yet "
        "interpolate or certify a generic child equation."
    ),
    "reproduce": (
        "python3 elkies-k3/scripts/batch_q80_third_q12_weierstrass_mod19_quadratic.py "
        f"--attempts {args.attempts} --workers {args.workers}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12WEIERSTRASSBATCH|attempts={len(candidates)}|"
    f"success={len(successes)}|failure={len(failures)}|"
    f"training={training_count}|heldout=8|"
    "status=PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_MOD19_QUADRATIC",
    flush=True,
)

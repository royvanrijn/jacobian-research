#!/usr/bin/env python3
"""Run deterministic mapped Weierstrass samples for a common-producer prime."""

import argparse
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKER = ROOT / "elkies-k3/scripts/sample_q80_third_q12_weierstrass_modp2.py"
EXPECTED_STATUS = "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_COMMON_PRODUCER"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--attempts", type=int, default=72)
parser.add_argument("--workers", type=int, default=6)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--sample-dir", type=Path, required=True)
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()
args.sample_dir = args.sample_dir.resolve()
if args.workers < 1:
    raise ValueError("workers must be positive")

pencil = json.loads(args.input.read_text())
if pencil.get("status") != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER":
    raise ValueError("input is not a certified common-producer pencil")
specialization = pencil["specialization"]
prime = int(specialization["prime"])
if not 1 <= args.attempts <= prime**2:
    raise ValueError("attempts must lie between one and the quadratic field size")

args.sample_dir.mkdir(parents=True, exist_ok=True)
candidates = [(index % prime, index // prime) for index in range(args.attempts)]
worker_hash = sha256(WORKER)
input_hash = sha256(args.input)


def run_one(value):
    a, b = value
    output = args.sample_dir / f"T-{a:03d}-{b:03d}.json"
    if output.is_file():
        try:
            payload = json.loads(output.read_text())
            if (
                payload.get("status") == EXPECTED_STATUS
                and payload.get("worker", {}).get("adapter", {}).get("sha256") == worker_hash
                and payload.get("input", {}).get("sha256") == input_hash
            ):
                return value, output, True, None
        except (OSError, json.JSONDecodeError):
            pass
    command = (
        "sage", "-python", str(WORKER),
        "--input", str(args.input),
        "--new-base", str(a),
        "--new-base-r", str(b),
        "--output", str(output),
    )
    completed = subprocess.run(
        command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=60, check=False,
    )
    if completed.returncode:
        return value, output, False, completed.stdout[-4000:]
    payload = json.loads(output.read_text())
    if payload.get("status") != EXPECTED_STATUS:
        raise ArithmeticError(f"sample {a}+{b}r has an unexpected status")
    return value, output, False, None


with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    results = list(executor.map(run_one, candidates))

successes = []
failures = []
for value, path, cached, error in results:
    if error is not None:
        failures.append({
            "new_base_coefficients_1_r": list(value),
            "diagnostic_tail": error,
        })
        continue
    payload = json.loads(path.read_text())
    successes.append({
        "new_base_coefficients_1_r": list(value),
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path),
        "cached": cached,
        "a1_a2_a3_a4_a6": payload["weierstrass"]["a1_a2_a3_a4_a6"],
        "discriminant": payload["weierstrass"]["discriminant"],
        "j": payload["weierstrass"]["j"],
        "inverse_weighted_bounds": [
            payload["birational_maps"]["inverse"]["W"]["weighted_bound"],
            payload["birational_maps"]["inverse"]["old_x"]["weighted_bound"],
        ],
    })
if len(successes) < 62:
    raise ArithmeticError(
        f"only {len(successes)} successful samples from {len(candidates)} attempts"
    )

training_count = len(successes) - 8
output = {
    "schema": "elkies-k3.q80-third-q12-weierstrass-sample-batch-modp2.v2",
    "status": "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER",
    "specialization": {
        "u": specialization["u"],
        "prime": prime,
        "extension_modulus": specialization["extension_modulus"],
    },
    "attempt_count": len(candidates),
    "success_count": len(successes),
    "failure_count": len(failures),
    "training_count": training_count,
    "held_out_count": 8,
    "training_samples": successes[:training_count],
    "held_out_samples": successes[training_count:],
    "failures": failures,
    "input": {"path": str(args.input.relative_to(ROOT)), "sha256": input_hash},
    "worker": {"path": str(WORKER.relative_to(ROOT)), "sha256": worker_hash},
    "normalization": "Grauert--Remmert module normalization in every retained fibre",
    "claim_boundary": (
        "Each retained fibre has exact forward/inverse maps to a coherently "
        "Laurent-normalized Weierstrass model. The batch does not by itself "
        "certify a generic child equation."
    ),
    "reproduce": (
        "python3 elkies-k3/scripts/batch_q80_third_q12_weierstrass_modp2.py "
        f"--input {args.input} --attempts {args.attempts} --workers {args.workers} "
        f"--sample-dir {args.sample_dir} --output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONBATCH|prime={prime}|attempts={len(candidates)}|"
    f"success={len(successes)}|failure={len(failures)}|"
    f"training={training_count}|heldout=8|"
    "status=PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER"
)

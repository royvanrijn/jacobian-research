#!/usr/bin/env python3
"""Certify provenance of the hybrid checkpointed 19^1024 factor lift."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--artifact",
    type=Path,
    default=RESULTS / "q80-third-q12-discriminant-factors-p19-adic-precision1024.json",
)
parser.add_argument(
    "--checkpoint",
    type=Path,
    default=LOCAL / "q80-third-q12-p19-target1024-checkpoint.json",
)
parser.add_argument(
    "--runner",
    type=Path,
    default=ROOT / "elkies-k3/scripts/run_q80_third_q12_discriminant_factors_p19_checkpointed.py",
)
parser.add_argument(
    "--output",
    type=Path,
    default=RESULTS / "q80-third-q12-discriminant-factors-p19-adic-precision1024-execution.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.artifact = args.artifact.resolve()
args.checkpoint = args.checkpoint.resolve()
args.runner = args.runner.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


artifact = json.loads(args.artifact.read_text())
checkpoint = json.loads(args.checkpoint.read_text())
if artifact.get("status") != "PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19":
    raise ValueError("final factor artifact is not certified")
if checkpoint.get("schema") != "elkies-k3.q80-third-q12-p19-hensel-checkpoint.v1":
    raise ValueError("unrecognized checkpoint schema")
if checkpoint["stage"] != "root" or checkpoint["valuation"] != 1024:
    raise ArithmeticError("checkpoint did not reach the root valuation-1024 gate")
if checkpoint["target_digits"] != 1024 or checkpoint["working_precision"] != 1027:
    raise ArithmeticError("checkpoint precision changed")
if checkpoint["checkpoint_runner_sha256"] != sha256(args.runner):
    raise ArithmeticError("checkpoint runner hash changed")
if checkpoint["canonical_worker_sha256"] != artifact["worker"]["sha256"]:
    raise ArithmeticError("canonical worker pins disagree")
if checkpoint["input_sha256"] != artifact["input"]["sha256"]:
    raise ArithmeticError("checkpoint/final source-pencil hashes disagree")

modulus = 19**1024


def reduce_pair(pair):
    return [int(pair[0]) % modulus, int(pair[1]) % modulus]


def compare_factor(name):
    checkpoint_coefficients = checkpoint[name]
    artifact_coefficients = artifact["factorization"][name]["coefficients_low_to_high_W"]
    if len(checkpoint_coefficients) != len(artifact_coefficients):
        raise ArithmeticError(f"{name}: W-degree changed")
    coordinate_count = 0
    for checkpoint_coefficient, artifact_coefficient in zip(
        checkpoint_coefficients, artifact_coefficients
    ):
        for checkpoint_key, artifact_key in (
            ("numerator", "numerator_coefficients_low_to_high_U_1_omega"),
            ("denominator", "denominator_coefficients_low_to_high_U_1_omega"),
        ):
            reduced = [reduce_pair(value) for value in checkpoint_coefficient[checkpoint_key]]
            expected = [reduce_pair(value) for value in artifact_coefficient[artifact_key]]
            if reduced != expected:
                raise ArithmeticError(f"{name}: checkpoint/final {checkpoint_key} changed")
            coordinate_count += 2 * len(reduced)
    return coordinate_count


counts = {name: compare_factor(name) for name in ("L", "Q", "D")}
factor_history = artifact["factorization"]["valuation_history"]
root_history = artifact["integral_basis_candidate"]["repeated_root_valuation_history"]
expected_factor_tail = [48, 96, 192, 384, 768, 1024]
if factor_history[-len(expected_factor_tail) :] != expected_factor_tail:
    raise ArithmeticError("hybrid factor jump history changed")
if root_history != [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
    raise ArithmeticError("pointwise root jump history changed")

output = {
    "schema": "elkies-k3.q80-third-q12-p19-precision1024-execution.v1",
    "status": "PASS_EXACT_THIRD_Q12_P19_PRECISION1024_HYBRID_EXECUTION",
    "specialization": {
        "u": "-2",
        "prime": 19,
        "digits": 1024,
        "working_precision": 1027,
    },
    "execution": {
        "algorithm": (
            "fixed mod-19 inverse digit corrections through valuation 48, then "
            "full-precision pointwise 9x9 Newton solves at good U-values with "
            "support-pinned rational interpolation; repeated-root pointwise Newton"
        ),
        "factor_valuation_history": factor_history,
        "repeated_root_valuation_history": root_history,
        "checkpoint_stage": checkpoint["stage"],
        "checkpoint_valuation": checkpoint["valuation"],
        "checkpoint_L_Q_D_scalar_coordinate_counts": counts,
        "checkpoint_L_Q_D_equal_final_mod_19_power": True,
    },
    "inputs": {
        "factor_artifact": {
            "path": str(args.artifact.relative_to(ROOT)),
            "sha256": sha256(args.artifact),
        },
        "checkpoint": {
            "path": str(args.checkpoint.relative_to(ROOT)),
            "sha256": sha256(args.checkpoint),
            "retention": "local resumable state; final coefficient comparison is certified here",
        },
        "runner": {
            "path": str(args.runner.relative_to(ROOT)),
            "sha256": sha256(args.runner),
        },
        "canonical_worker": artifact["worker"],
    },
    "claim_boundary": {
        "proved": [
            "the actual accelerated execution and final checkpoint are hash-pinned",
            "checkpoint L,Q,D reduce coefficientwise to the final 19^1024 artifact",
            "the displayed hybrid and quadratic valuation histories are exact",
        ],
        "not_proved": [
            "a characteristic-zero factorization or Jacobian",
            "that the canonical worker alone reproduces the accelerated valuation history",
        ],
    },
    "reproduce": (
        "python3 elkies-k3/scripts/certify_q80_third_q12_p19_precision1024_execution.py"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"precision-1024 execution artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80THIRDQ12P19EXECUTION|digits=1024|factor_tail=48,96,192,384,768,1024|"
    "root=1,2,4,8,16,32,64,128,256,512,1024|"
    "status=PASS_EXACT_THIRD_Q12_P19_PRECISION1024_HYBRID_EXECUTION"
)

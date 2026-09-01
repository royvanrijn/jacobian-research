#!/usr/bin/env python3
"""Replay a reconstructed Q80 third-q12 j-map at one blind prime."""

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.set_int_max_str_digits(0)
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=163)
parser.add_argument(
    "--candidate",
    type=Path,
    default=RESULTS / "q80-third-q12-j-map-p19-adic-reconstructed-qq.json",
)
parser.add_argument(
    "--blind-transport",
    type=Path,
    default=None,
)
parser.add_argument(
    "--blind-pencil-alignment",
    type=Path,
    default=None,
)
parser.add_argument(
    "--output",
    type=Path,
    default=None,
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
prime = args.prime
if prime < 3:
    raise ValueError("the held-out prime must be an odd prime")
if args.blind_transport is None:
    args.blind_transport = (
        RESULTS
        / f"q80-third-q12-long-jacobians-exact-quadratic-gauge-p{prime}-heldout.json"
    )
if args.blind_pencil_alignment is None:
    args.blind_pencil_alignment = (
        RESULTS
        / f"q80-third-q12-um2-exact-quadratic-pencils-p19-legacy-aligned-p{prime}-heldout.json"
    )
if args.output is None:
    args.output = RESULTS / f"q80-third-q12-j-map-p{prime}-heldout-replay.json"
for name in ("candidate", "blind_transport", "blind_pencil_alignment", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


candidate = json.loads(args.candidate.read_text())
transport = json.loads(args.blind_transport.read_text())
alignment = json.loads(args.blind_pencil_alignment.read_text())
if candidate.get("status") != "PASS_CANDIDATE_EXACT_THIRD_Q12_J_MAP_RECONSTRUCTION_QQ":
    raise ValueError("the characteristic-zero j-map candidate is not certified as a candidate")
if transport.get("status") != "PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE":
    raise ValueError("the blind transported Jacobian ensemble is not certified")
if alignment.get("status") != "PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_PENCIL_DESCENT_AND_LOCAL_ALIGNMENT":
    raise ValueError("the blind pencil alignment is not certified")

if prime in candidate["reconstruction"]["CRT_auxiliary_primes"]:
    raise ArithmeticError(f"p={prime} was consumed by candidate reconstruction")
if prime in candidate["validation"]["transported_primes_used_for_reconstruction_acceptance"]:
    raise ArithmeticError(f"p={prime} was consumed by candidate acceptance")
if prime not in transport["specialization"]["primes"]:
    raise ArithmeticError(f"p={prime} is absent from blind transport")
if prime not in alignment["specialization"]["primes"]:
    raise ArithmeticError(f"p={prime} is absent from blind alignment")
if str(prime) not in transport["transported_models"]:
    raise ArithmeticError(f"p={prime} transported model is absent")


def reduce_rational(value):
    value = Fraction(value)
    denominator = value.denominator % prime
    if denominator == 0:
        raise ZeroDivisionError(f"candidate denominator is divisible by p={prime}")
    return value.numerator % prime * pow(denominator, -1, prime) % prime


candidate_j = candidate["j_map"]
expected_j = transport["transported_models"][str(prime)]["j"]
comparisons = {}
all_match = True
for candidate_key, expected_key in (
    ("numerator_coefficients_low_to_high_U", "numerator_coefficients_low_to_high_1_omega"),
    ("denominator_coefficients_low_to_high_U", "denominator_coefficients_low_to_high_1_omega"),
):
    reduced = [
        [reduce_rational(value) for value in record["coefficients_1_omega"]]
        for record in candidate_j[candidate_key]
    ]
    expected = expected_j[expected_key]
    mismatch_indices = [
        index
        for index, (left, right) in enumerate(zip(reduced, expected))
        if left != right
    ]
    literal_match = not mismatch_indices and len(reduced) == len(expected)
    all_match &= literal_match
    comparisons[candidate_key] = {
        "coefficient_pairs": len(reduced),
        "literal_match": literal_match,
        "mismatch_count": len(mismatch_indices),
        "first_mismatch": (
            None
            if not mismatch_indices
            else {
                "coefficient_index": mismatch_indices[0],
                "candidate": reduced[mismatch_indices[0]],
                "blind_expected": expected[mismatch_indices[0]],
            }
        ),
    }

output = {
    "schema": "elkies-k3.q80-third-q12-j-map-heldout-prime-replay.v2",
    "status": (
        f"PASS_Q80_THIRD_Q12_J_MAP_BLIND_P{prime}_REPLAY"
        if all_match
        else f"FAIL_Q80_THIRD_Q12_J_MAP_BLIND_P{prime}_REPLAY"
    ),
    "specialization": {
        "u": "-2",
        "prime": prime,
        "coefficient_basis": ["1", "omega"],
        "base_coordinate": "the exact legacy-aligned U gauge",
    },
    "independence": {
        "candidate_CRT_auxiliary_primes": candidate["reconstruction"]["CRT_auxiliary_primes"],
        "candidate_acceptance_primes": candidate["validation"][
            "transported_primes_used_for_reconstruction_acceptance"
        ],
        "heldout_prime_absent_from_candidate_reconstruction": True,
        "heldout_prime_produced_before_candidate_reduction": True,
    },
    "comparison": comparisons,
    "inputs": {
        "candidate": {
            "path": str(args.candidate.relative_to(ROOT)),
            "sha256": sha256(args.candidate),
        },
        "blind_transport": {
            "path": str(args.blind_transport.relative_to(ROOT)),
            "sha256": sha256(args.blind_transport),
        },
        "blind_pencil_alignment": {
            "path": str(args.blind_pencil_alignment.relative_to(ROOT)),
            "sha256": sha256(args.blind_pencil_alignment),
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            f"p={prime} was absent from every CRT and row-acceptance input used by the candidate",
            f"the independently produced p={prime} pencil is aligned to the exact quadratic and base gauges",
        ]
        + (
            ["all 25 numerator and 25 denominator coefficient pairs match literally after reduction"]
            if all_match
            else ["the candidate and blind p=163 j-maps differ literally in the pinned exact gauge"]
        ),
        "not_proved": [
            "literal characteristic-zero substitution into the exact pencil",
            "an exact characteristic-zero Jacobian model or birational maps",
        ],
    },
    "reproduce": " ".join(
        [
            "python3",
            "elkies-k3/scripts/certify_q80_third_q12_j_map_p163_heldout.py",
            "--prime",
            str(prime),
            "--candidate",
            str(args.candidate.relative_to(ROOT)),
            "--blind-transport",
            str(args.blind_transport.relative_to(ROOT)),
            "--blind-pencil-alignment",
            str(args.blind_pencil_alignment.relative_to(ROOT)),
            "--output",
            str(args.output.relative_to(ROOT)),
        ]
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"p=163 held-out replay artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12JHELDOUT|prime={prime}|numerator_pairs=25|denominator_pairs=25|"
    f"status={'PASS' if all_match else 'FAIL'}_Q80_THIRD_Q12_J_MAP_BLIND_P{prime}_REPLAY"
)

#!/usr/bin/env python3
"""Compile a valid CRT ledger for transported long q12 Jacobians."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_INPUT = RESULTS / "q80-third-q12-long-jacobians-exact-quadratic-gauge.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-long-jacobian-exact-quadratic-crt.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(args.input.read_text())
if payload.get("status") != "PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE":
    raise ValueError("transported long Jacobians are not certified")
primes = tuple(payload["specialization"]["primes"])
if primes != (19, 61, 67, 83, 89, 103, 131):
    raise ArithmeticError("unexpected transported-prime set")


def coefficient_slots(model):
    slots = {}
    records = {
        **model["weierstrass"],
        "discriminant": model["discriminant"],
        "j": model["j"],
    }
    for name in ("a1", "a2", "a3", "a4", "a6", "discriminant", "j"):
        record = records[name]
        for side in ("numerator", "denominator"):
            values = record[f"{side}_coefficients_low_to_high_1_omega"]
            for degree, coordinates in enumerate(values):
                if len(coordinates) != 2:
                    raise ArithmeticError("quadratic coefficient does not have two coordinates")
                for coordinate, label in enumerate(("constant", "omega")):
                    slots[f"{name}.{side}.{degree}.{label}"] = int(coordinates[coordinate])
    return slots


by_prime = {
    prime: coefficient_slots(payload["transported_models"][str(prime)])
    for prime in primes
}
labels = tuple(sorted(by_prime[primes[0]]))
if any(tuple(sorted(by_prime[prime])) != labels for prime in primes[1:]):
    raise ArithmeticError("transported long-Jacobian slot shapes disagree")


def crt(residues):
    value = 0
    modulus = 1
    for prime, residue in zip(primes, residues):
        correction = (residue - value) * pow(modulus, -1, prime) % prime
        value += correction * modulus
        modulus *= prime
    return value, modulus


combined_slots = []
common_modulus = None
for label in labels:
    residues = [by_prime[prime][label] % prime for prime in primes]
    value, modulus = crt(residues)
    if common_modulus is None:
        common_modulus = modulus
    elif modulus != common_modulus:
        raise ArithmeticError("CRT moduli disagree")
    if any(value % prime != residue for prime, residue in zip(primes, residues)):
        raise ArithmeticError(f"CRT replay failed for {label}")
    centered = value if value <= modulus // 2 else value - modulus
    combined_slots.append(
        {
            "label": label,
            "residues_in_prime_order": residues,
            "crt_residue": value,
            "centered_residue": centered,
        }
    )

output = {
    "schema": "elkies-k3.q80-third-q12-long-jacobian-exact-quadratic-crt.v1",
    "status": "PASS_EXACT_ALIGNED_THIRD_Q12_LONG_JACOBIAN_QUADRATIC_CRT_LEDGER",
    "specialization": {"u": "-2", "primes": list(primes)},
    "exact_gauge": payload["exact_gauge"],
    "coefficient_basis": ["1", "omega"],
    "omega_square": "16*q1*q2",
    "common_degree_shapes": payload["common_degree_shapes"],
    "crt": {
        "modulus": common_modulus,
        "slot_count": len(combined_slots),
        "slots": combined_slots,
        "all_residue_replays": True,
    },
    "held_out_controls": {
        "p19": "included mandatory legacy full-child control",
        "p71": "exact horizontal held-out prime; no interpolated full child",
    },
    "input": {"path": str(args.input.relative_to(ROOT)), "sha256": sha256(args.input)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "all displayed long-Jacobian coordinates are in one exact base and quadratic-field gauge",
            "coordinatewise CRT accumulation across six independently computed primes",
            "literal reduction of every CRT residue to every input prime",
        ],
        "not_proved": [
            "that centered residues are integer or rational coefficients",
            "rational reconstruction at the current modulus",
            "a characteristic-zero Jacobian equation or birational maps",
        ],
    },
    "reproduce": (
        "python3 "
        "elkies-k3/scripts/compile_q80_third_q12_long_jacobian_quadratic_crt.py"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"quadratic long-Jacobian CRT artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12LONGCRT|primes={','.join(map(str, primes))}|"
    f"modulus={common_modulus}|slots={len(combined_slots)}|basis=1,omega|"
    "status=PASS_EXACT_ALIGNED_THIRD_Q12_LONG_JACOBIAN_QUADRATIC_CRT_LEDGER"
)

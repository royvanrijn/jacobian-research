#!/usr/bin/env python3
"""Compile a valid CRT ledger for transported two-way q12 maps."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_INPUT = RESULTS / "q80-third-q12-birational-maps-exact-quadratic-gauge.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-birational-maps-exact-quadratic-crt.json"
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
if payload.get("status") != "PASS_EXACT_TRANSPORTED_THIRD_Q12_BIRATIONAL_MAPS_COMMON_QUADRATIC_GAUGE":
    raise ValueError("transported birational maps are not certified")
primes = tuple(payload["specialization"]["primes"])
if primes != (19, 61, 67, 83, 89, 103, 131):
    raise ArithmeticError("unexpected transported-prime set")


def slots(value, path=()):
    result = {}
    if isinstance(value, dict):
        for key, child in sorted(value.items()):
            if key.endswith("_coefficients_low_to_high_1_omega"):
                for degree, coordinates in enumerate(child):
                    if len(coordinates) != 2:
                        raise ArithmeticError("map coefficient does not have two quadratic coordinates")
                    result[".".join(path + (key, str(degree), "constant"))] = int(coordinates[0])
                    result[".".join(path + (key, str(degree), "omega"))] = int(coordinates[1])
            else:
                result.update(slots(child, path + (key,)))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            result.update(slots(child, path + (str(index),)))
    return result


by_prime = {
    prime: slots(payload["transported_maps"][str(prime)]) for prime in primes
}
labels = tuple(sorted(by_prime[primes[0]]))
if any(tuple(sorted(by_prime[prime])) != labels for prime in primes[1:]):
    raise ArithmeticError("transported map slot shapes disagree")


def crt(residues):
    value = 0
    modulus = 1
    for prime, residue in zip(primes, residues):
        correction = (residue - value) * pow(modulus, -1, prime) % prime
        value += correction * modulus
        modulus *= prime
    return value, modulus


combined = []
common_modulus = None
for label in labels:
    residues = [by_prime[prime][label] % prime for prime in primes]
    value, modulus = crt(residues)
    common_modulus = modulus if common_modulus is None else common_modulus
    if modulus != common_modulus:
        raise ArithmeticError("CRT moduli disagree")
    if any(value % prime != residue for prime, residue in zip(primes, residues)):
        raise ArithmeticError(f"CRT replay failed for {label}")
    combined.append(
        {
            "label": label,
            "residues_in_prime_order": residues,
            "crt_residue": value,
            "centered_residue": value if value <= modulus // 2 else value - modulus,
        }
    )

output = {
    "schema": "elkies-k3.q80-third-q12-birational-maps-exact-quadratic-crt.v1",
    "status": "PASS_EXACT_ALIGNED_THIRD_Q12_BIRATIONAL_MAPS_QUADRATIC_CRT_LEDGER",
    "specialization": {"u": "-2", "primes": list(primes)},
    "exact_gauge": payload["exact_gauge"],
    "coefficient_basis": ["1", "omega"],
    "common_shape": payload["common_shape"],
    "crt": {
        "modulus": common_modulus,
        "slot_count": len(combined),
        "slots": combined,
        "all_residue_replays": True,
    },
    "held_out_controls": {
        "p19": "included mandatory legacy complete two-way-map control",
        "p71": "exact horizontal held-out prime; no interpolated full child",
    },
    "input": {"path": str(args.input.relative_to(ROOT)), "sha256": sha256(args.input)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "coordinatewise CRT for transported forward and inverse generic map coefficients",
            "one exact quadratic-field basis and base coordinate at all six primes",
            "literal reduction of every CRT slot to every input prime",
        ],
        "not_proved": [
            "rational reconstruction at the current modulus",
            "characteristic-zero forward or inverse map identities",
            "a characteristic-zero Jacobian equation",
        ],
    },
    "reproduce": (
        "python3 elkies-k3/scripts/compile_q80_third_q12_maps_quadratic_crt.py"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"quadratic map CRT artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12MAPCRT|primes={','.join(map(str, primes))}|modulus={common_modulus}|"
    f"slots={len(combined)}|maps=forward,inverse|"
    "status=PASS_EXACT_ALIGNED_THIRD_Q12_BIRATIONAL_MAPS_QUADRATIC_CRT_LEDGER"
)

#!/usr/bin/env sage-python
"""Reconstruct the compact q8/orbit376 horizontal over QQ by modular CRT.

The input files are exact good-prime outputs of
``identify_h92_q4o164_q8_horizontal_mod131.sage``.  Their fraction-field
coordinates have monic denominators, so every coefficient has a canonical
residue.  We combine those residues by CRT, try Sage rational reconstruction
coefficientwise, and recover each whole coordinate projectively by a small
LLL when the shared normalization scale is too large.  A QQ result is accepted
only after exact substitution in the compact Weierstrass equation and
reduction back to every input prime.

No Groebner basis is used.
"""

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from sage.all import CRT_list, PolynomialRing, QQ, ZZ, gcd, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o164-q8o376-horizontal-crt-qq.json",
)
args = parser.parse_args()
inputs = args.inputs or sorted(
    LOCAL.glob("q4o164-q8o376-horizontal-from-abel-trace-mod*.json")
)
inputs = [path if path.is_absolute() else ROOT / path for path in inputs]
output = args.output if args.output.is_absolute() else ROOT / args.output
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pad(values, size):
    assert len(values) <= size
    return list(values) + [0] * (size - len(values))


assert inputs, "no modular q8 horizontal artifacts found"
records = [json.loads(path.read_text()) for path in inputs]
primes = [ZZ(record["prime"]) for record in records]
assert len(set(primes)) == len(primes)
assert all(prime.is_prime() for prime in primes)
assert all(
    record["status"] == "PASS_EXACT_MODP_Q4O164_Q8O376_HORIZONTAL_FROM_ABEL_TRACE"
    for record in records
)
assert len({record["selected_identity"] for record in records}) == 1
assert all(
    record["selected"]["embedding_index"] == 15
    and record["selected"]["c8_opposite_coefficient"] == -1
    and record["selected"]["integral_B0_through_B7_word"] == [-1, 2, 1, -3, -1, -2, 0, 1]
    and record["selected"]["coordinate_degrees_x_num_x_den_y_num_y_den"] == [12, 8, 18, 12]
    for record in records
)

parts = (("x", "numerator", 13), ("x", "denominator", 9),
         ("y", "numerator", 19), ("y", "denominator", 13))
residues = {}
for coordinate, part, size in parts:
    key = f"{coordinate}_{part}"
    residues[key] = [
        pad(record["selected"][coordinate][f"{part}_coefficients_low_to_high"], size)
        for record in records
    ]
    if part == "denominator":
        assert all(values[-1] % prime == 1 for values, prime in zip(residues[key], primes))

modulus = math.prod(primes)
reconstructed = {}
failures = []
crt_coefficients = {}
for coordinate, part, size in parts:
    key = f"{coordinate}_{part}"
    coefficients = []
    for index in range(size):
        residue = ZZ(CRT_list([values[index] for values in residues[key]], primes))
        crt_coefficients.setdefault(key, []).append(residue)
        try:
            coefficient = residue.rational_reconstruction(modulus)
        except ArithmeticError as error:
            failures.append({"part": key, "index": index, "reason": str(error)})
            coefficient = None
        coefficients.append(coefficient)
    reconstructed[key] = coefficients

successful_coefficients = [
    coefficient
    for coefficients in reconstructed.values()
    for coefficient in coefficients
    if coefficient is not None
]
reconstruction_progress = {
    "total_coefficient_count": sum(size for _, _, size in parts),
    "successful_coefficient_count": len(successful_coefficients),
    "maximum_successful_numerator_bits": max(
        [abs(coefficient.numerator()).nbits() for coefficient in successful_coefficients],
        default=0,
    ),
    "maximum_successful_denominator_bits": max(
        [coefficient.denominator().nbits() for coefficient in successful_coefficients],
        default=0,
    ),
}


def simultaneous_candidates(coordinate, numerator_size, denominator_size):
    """Recover a primitive projective coefficient vector with a small LLL."""
    normalized = (
        crt_coefficients[f"{coordinate}_numerator"]
        + crt_coefficients[f"{coordinate}_denominator"][:-1]
    )
    dimension = len(normalized) + 1
    basis = matrix(ZZ, dimension, dimension)
    for index in range(len(normalized)):
        basis[index, index] = modulus
    basis[-1] = vector(ZZ, normalized + [1])
    reduced = basis.LLL(delta=0.99)
    candidates = []
    for row in sorted(reduced.rows(), key=lambda item: item * item):
        scale = ZZ(row[-1])
        if not scale or gcd(scale, modulus) != 1:
            continue
        if not all((row[index] - scale * normalized[index]) % modulus == 0
                   for index in range(len(normalized))):
            continue
        numerator = [QQ(row[index]) / scale for index in range(numerator_size)]
        denominator = [
            QQ(row[numerator_size + index]) / scale
            for index in range(denominator_size - 1)
        ] + [QQ.one()]
        candidates.append({
            "numerator": numerator,
            "denominator": denominator,
            "projective_scale": scale,
            "primitive_vector_max_bits": max(abs(ZZ(value)).nbits() for value in row),
            "euclidean_norm_squared": str(row * row),
        })
    return candidates


simultaneous = {
    "x": simultaneous_candidates("x", 13, 9),
    "y": simultaneous_candidates("y", 19, 13),
}

model = json.loads(MODEL.read_text())
R = PolynomialRing(QQ, "t")
t = R.gen()
equation_verified = False
reductions_verified = False
degree_verified = False
height_verified = False
fourfold_pole_degree = None
canonical_height = None
x = y = None
A = R([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])


def verify_reductions(x_value, y_value):
    verified = True
    for record, prime in zip(records, primes):
        for coordinate, value in (("x", x_value), ("y", y_value)):
            for part, polynomial in (("numerator", value.numerator()),
                                     ("denominator", value.denominator())):
                expected = record["selected"][coordinate][f"{part}_coefficients_low_to_high"]
                actual = [int((coefficient.numerator() * pow(int(coefficient.denominator()), -1, int(prime))) % prime)
                          for coefficient in polynomial.list()]
                if actual != expected:
                    verified = False
    return verified


candidate_pairs_tested = 0
if not failures:
    independent_x = R(reconstructed["x_numerator"]) / R(reconstructed["x_denominator"])
    independent_y = R(reconstructed["y_numerator"]) / R(reconstructed["y_denominator"])
    simultaneous["x"].insert(0, {
        "numerator": list(independent_x.numerator()),
        "denominator": list(independent_x.denominator()),
        "projective_scale": 1,
        "primitive_vector_max_bits": 0,
        "euclidean_norm_squared": "0",
    })
    simultaneous["y"].insert(0, {
        "numerator": list(independent_y.numerator()),
        "denominator": list(independent_y.denominator()),
        "projective_scale": 1,
        "primitive_vector_max_bits": 0,
        "euclidean_norm_squared": "0",
    })

for x_candidate in simultaneous["x"]:
    x_try = R(x_candidate["numerator"]) / R(x_candidate["denominator"])
    if [x_try.numerator().degree(), x_try.denominator().degree()] != [12, 8]:
        continue
    for y_candidate in simultaneous["y"]:
        y_try = R(y_candidate["numerator"]) / R(y_candidate["denominator"])
        candidate_pairs_tested += 1
        if [y_try.numerator().degree(), y_try.denominator().degree()] != [18, 12]:
            continue
        if y_try**2 != x_try**3 + A * x_try + B:
            continue
        if not verify_reductions(x_try, y_try):
            continue
        x, y = x_try, y_try
        equation_verified = reductions_verified = degree_verified = True
        break
    if equation_verified:
        break


def double_point(x_value, y_value):
    slope = (3 * x_value**2 + A) / (2 * y_value)
    x_result = slope**2 - 2 * x_value
    y_result = slope * (x_value - x_result) - y_value
    assert y_result**2 == x_result**3 + A * x_result + B
    return x_result, y_result


if equation_verified:
    x2, y2 = double_point(x, y)
    x4, unused_y4 = double_point(x2, y2)
    fourfold_pole_degree = max(
        x4.denominator().degree(), x4.numerator().degree() - 4
    )
    canonical_height = QQ(4 + fourfold_pole_degree) / 16
    height_verified = canonical_height == 11

passed = equation_verified and reductions_verified and degree_verified and height_verified


def qq_record(value):
    return {
        "numerator_coefficients_low_to_high": [str(coefficient) for coefficient in value.numerator().list()],
        "denominator_coefficients_low_to_high": [str(coefficient) for coefficient in value.denominator().list()],
        "degrees_numerator_denominator": [
            int(value.numerator().degree()), int(value.denominator().degree())
        ],
    }


payload = {
    "schema": "elkies-k3.q4o164-q8o376-horizontal-crt-qq.v1",
    "status": (
        "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
        if passed else "PARTIAL_Q4O164_Q8O376_HORIZONTAL_CRT"
    ),
    "primes": list(map(int, primes)),
    "crt_modulus": str(modulus),
    "crt_modulus_bits": int(ZZ(modulus).nbits()),
    "selected_identity": records[0]["selected_identity"],
    "rational_reconstruction_failures": failures,
    "rational_reconstruction_progress": reconstruction_progress,
    "simultaneous_projective_reconstruction": {
        "x_candidate_count": len(simultaneous["x"]),
        "y_candidate_count": len(simultaneous["y"]),
        "candidate_pairs_tested": candidate_pairs_tested,
        "selected_x_projective_scale": str(x_candidate["projective_scale"]) if passed else None,
        "selected_y_projective_scale": str(y_candidate["projective_scale"]) if passed else None,
        "selected_x_primitive_vector_max_bits": x_candidate["primitive_vector_max_bits"] if passed else None,
        "selected_y_primitive_vector_max_bits": y_candidate["primitive_vector_max_bits"] if passed else None,
    },
    "checks": {
        "exact_QQ_weierstrass_identity": equation_verified,
        "exact_degree_fingerprint": degree_verified,
        "exact_fourfold_pole_height": height_verified,
        "fourfold_x_pole_degree": int(fourfold_pole_degree) if fourfold_pole_degree is not None else None,
        "canonical_height": str(canonical_height) if canonical_height is not None else None,
        "reduction_to_every_input_artifact": reductions_verified,
        "large_Groebner_required": False,
    },
    "section": {"x": qq_record(x), "y": qq_record(y)} if passed else None,
    "proof_boundary": (
        "Exact CRT with simultaneous projective LLL reconstruction of the compact q8 horizontal, "
        "accepted only after exact QQ(t) substitution and replay at every input prime. "
        "This proves the horizontal section, not yet the resolved q8 Riemann--Roch pencil "
        "or the child 4A1/MW13 equation."
        if passed else
        "The available good-prime modulus does not yet produce a QQ section passing exact "
        "substitution and all modular replays. This is reconstruction progress, not a proof."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "model_path": str(MODEL.relative_to(ROOT)),
        "model_sha256": sha256(MODEL),
    },
    "runtime_seconds": time.monotonic() - started,
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164Q8CRT|primes={}|bits={}|failures={}|equation={}|reductions={}|status={}|output={}".format(
        list(map(int, primes)), ZZ(modulus).nbits(), len(failures), equation_verified,
        reductions_verified, payload["status"], output,
    ),
    flush=True,
)

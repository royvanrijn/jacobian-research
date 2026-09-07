#!/usr/bin/env sage
"""Projectively CRT-reconstruct the four q52 Abel-word sections over QQ.

Each modular input contains the old P.O=15 section used by the q52 RR divisor
and the three horizontal 5A1 roots (P.O=1,1,3).  Monic-denominator
normalization gives canonical residues prime by prime.  We reconstruct each
whole rational-function coordinate by a small projective LLL, accepting a
candidate only after literal QQ(t) substitution and replay at every input
prime.  No Groebner basis or nonlinear characteristic-zero solve is used.
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
MODEL = LOCAL / "fixed-reverse-4a1-rr-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "fixed-reverse-5a1-sections-crt-qq.json",
)
args = parser.parse_args()
inputs = args.inputs or sorted(
    path for path in LOCAL.glob("fixed-reverse-5a1-abel-word-seeds-mod*.json")
    if path.stem.rsplit("mod", 1)[1].isdigit()
    and int(path.stem.rsplit("mod", 1)[1]) > 1000000
)
inputs = [path if path.is_absolute() else ROOT / path for path in inputs]
output = args.output if args.output.is_absolute() else ROOT / args.output
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pad(values, size):
    assert len(values) <= size
    return list(values) + [0] * (size - len(values))


assert inputs, "no modular q52 Abel-word artifacts found"
records = [json.loads(path.read_text()) for path in inputs]
primes = [ZZ(record["prime"]) for record in records]
assert len(set(primes)) == len(primes)
assert all(prime.is_prime() for prime in primes)
assert all(
    record["status"] == "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_REGULAR_SEEDS"
    for record in records
)
assert len({json.dumps(record["words"]) for record in records}) == 1
assert all(len(record["sections"]) == 4 for record in records)
expected_poles = [15, 1, 1, 3]
assert all(
    [section["P_dot_O"] for section in record["sections"]] == expected_poles
    for record in records
)

modulus = math.prod(primes)
model = json.loads(MODEL.read_text())
R = PolynomialRing(QQ, "t")
t = R.gen()
A = R([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])


def normalized_modular_part(record, section_index, coordinate, part, size):
    section = record["sections"][section_index]
    numerator = pad(
        section[f"{coordinate}_numerator_coefficients_low_to_high"],
        2 * expected_poles[section_index] + 5
        if coordinate == "x" else 3 * expected_poles[section_index] + 7,
    )
    denominator = pad(
        section[f"{coordinate}_denominator_coefficients_low_to_high"],
        2 * expected_poles[section_index] + 1
        if coordinate == "x" else 3 * expected_poles[section_index] + 1,
    )
    prime = ZZ(record["prime"])
    scale = pow(int(denominator[-1]), -1, int(prime))
    values = numerator if part == "numerator" else denominator
    normalized = [(int(value) * scale) % int(prime) for value in values]
    assert normalized_modular_degree(normalized, prime) == size - 1
    return normalized


def normalized_modular_degree(values, prime):
    for index in range(len(values) - 1, -1, -1):
        if values[index] % prime:
            return index
    return -1


def coordinate_sizes(section_index, coordinate):
    pole = expected_poles[section_index]
    if coordinate == "x":
        return 2 * pole + 5, 2 * pole + 1
    return 3 * pole + 7, 3 * pole + 1


def crt_coordinate(section_index, coordinate):
    numerator_size, denominator_size = coordinate_sizes(section_index, coordinate)
    residues = {part: [
        normalized_modular_part(
            record, section_index, coordinate, part,
            numerator_size if part == "numerator" else denominator_size,
        )
        for record in records
    ] for part in ("numerator", "denominator")}
    crt = {}
    for part, size in (("numerator", numerator_size),
                       ("denominator", denominator_size)):
        crt[part] = [
            ZZ(CRT_list([values[index] for values in residues[part]], primes))
            for index in range(size)
        ]
    return crt, residues


def projective_candidates(crt, numerator_size, denominator_size):
    normalized = crt["numerator"] + crt["denominator"][:-1]
    dimension = len(normalized) + 1
    basis = matrix(ZZ, dimension, dimension)
    for index in range(len(normalized)):
        basis[index, index] = modulus
    basis[-1] = vector(ZZ, normalized + [1])
    reduced = basis.LLL(delta=0.99)
    candidates = []
    for row in sorted(reduced.rows(), key=lambda value: value * value):
        scale = ZZ(row[-1])
        if not scale or gcd(scale, modulus) != 1:
            continue
        if not all(
            (row[index] - scale * normalized[index]) % modulus == 0
            for index in range(len(normalized))
        ):
            continue
        numerator = [QQ(row[index]) / scale for index in range(numerator_size)]
        denominator = [
            QQ(row[numerator_size + index]) / scale
            for index in range(denominator_size - 1)
        ] + [QQ.one()]
        candidates.append({
            "value": R(numerator) / R(denominator),
            "scale": scale,
            "primitive_max_bits": max(abs(ZZ(entry)).nbits() for entry in row),
            "norm_squared": ZZ(row * row),
        })
    return candidates


def independent_candidate(crt):
    reconstructed = {}
    for part in ("numerator", "denominator"):
        values = []
        for residue in crt[part]:
            try:
                values.append(QQ(residue.rational_reconstruction(modulus)))
            except (ArithmeticError, ValueError):
                return None
        reconstructed[part] = values
    return {
        "value": R(reconstructed["numerator"]) / R(reconstructed["denominator"]),
        "scale": ZZ.one(),
        "primitive_max_bits": 0,
        "norm_squared": ZZ.zero(),
    }


def normalized_reduction(value, prime):
    numerator = list(value.numerator())
    denominator = list(value.denominator())
    leading = denominator[-1]
    scale = (
        int(leading.numerator())
        * pow(int(leading.denominator()), -1, int(prime))
    ) % int(prime)
    scale = pow(scale, -1, int(prime))

    def reduce(values):
        return [
            (
                int(coefficient.numerator())
                * pow(int(coefficient.denominator()), -1, int(prime))
                * scale
            ) % int(prime)
            for coefficient in values
        ]

    return reduce(numerator), reduce(denominator)


def verify_reductions(section_index, coordinate, value, residues):
    numerator_size, denominator_size = coordinate_sizes(section_index, coordinate)
    for input_index, prime in enumerate(primes):
        actual_num, actual_den = normalized_reduction(value, prime)
        if pad(actual_num, numerator_size) != residues["numerator"][input_index]:
            return False
        if pad(actual_den, denominator_size) != residues["denominator"][input_index]:
            return False
    return True


def qq_record(value):
    return {
        "numerator_coefficients_low_to_high": list(map(str, value.numerator().list())),
        "denominator_coefficients_low_to_high": list(map(str, value.denominator().list())),
        "degrees_numerator_denominator": [
            int(value.numerator().degree()), int(value.denominator().degree())
        ],
    }


section_results = []
all_passed = True
for section_index, pole in enumerate(expected_poles):
    coordinate_data = {}
    for coordinate in ("x", "y"):
        crt, residues = crt_coordinate(section_index, coordinate)
        numerator_size, denominator_size = coordinate_sizes(section_index, coordinate)
        candidates = projective_candidates(crt, numerator_size, denominator_size)
        independent = independent_candidate(crt)
        if independent is not None:
            candidates.insert(0, independent)
        coordinate_data[coordinate] = {
            "crt": crt,
            "residues": residues,
            "candidates": candidates,
        }
    selected = None
    pairs_tested = 0
    expected_degrees = {
        "x": [2 * pole + 4, 2 * pole],
        "y": [3 * pole + 6, 3 * pole],
    }
    for x_candidate in coordinate_data["x"]["candidates"]:
        x_value = x_candidate["value"]
        if [x_value.numerator().degree(), x_value.denominator().degree()] != expected_degrees["x"]:
            continue
        if not verify_reductions(
            section_index, "x", x_value, coordinate_data["x"]["residues"]
        ):
            continue
        for y_candidate in coordinate_data["y"]["candidates"]:
            y_value = y_candidate["value"]
            pairs_tested += 1
            if [y_value.numerator().degree(), y_value.denominator().degree()] != expected_degrees["y"]:
                continue
            if y_value ** 2 != x_value ** 3 + A * x_value + B:
                continue
            if not verify_reductions(
                section_index, "y", y_value, coordinate_data["y"]["residues"]
            ):
                continue
            selected = (x_candidate, y_candidate)
            break
        if selected is not None:
            break
    passed = selected is not None
    all_passed = all_passed and passed
    result = {
        "section_index": section_index,
        "P_dot_O": pole,
        "x_projective_candidate_count": len(coordinate_data["x"]["candidates"]),
        "y_projective_candidate_count": len(coordinate_data["y"]["candidates"]),
        "candidate_pairs_tested": pairs_tested,
        "exact_QQ_weierstrass_identity": passed,
        "reduction_to_every_input_artifact": passed,
    }
    if passed:
        x_candidate, y_candidate = selected
        result.update({
            "x": qq_record(x_candidate["value"]),
            "y": qq_record(y_candidate["value"]),
            "x_projective_scale": str(x_candidate["scale"]),
            "y_projective_scale": str(y_candidate["scale"]),
            "x_primitive_vector_max_bits": x_candidate["primitive_max_bits"],
            "y_primitive_vector_max_bits": y_candidate["primitive_max_bits"],
        })
    section_results.append(result)
    print(
        "FIXEDREVERSE5A1CRT|section={}|pole={}|x_candidates={}|y_candidates={}|passed={}|seconds={:.3f}".format(
            section_index, pole,
            result["x_projective_candidate_count"],
            result["y_projective_candidate_count"], int(passed),
            time.monotonic() - started,
        ),
        flush=True,
    )

payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-sections-crt-qq.v1",
    "status": (
        "PASS_EXACT_QQ_FIXED_REVERSE_5A1_SECTIONS_CRT"
        if all_passed else "PARTIAL_FIXED_REVERSE_5A1_SECTIONS_CRT"
    ),
    "primes": list(map(int, primes)),
    "crt_modulus_bits": int(ZZ(modulus).nbits()),
    "words": records[0]["words"],
    "sections": section_results,
    "method": {
        "canonical_monic_denominator_residues": True,
        "simultaneous_projective_LLL": True,
        "groebner_or_nonlinear_QQ_solve": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Exact projective CRT reconstruction, literal QQ(t) substitution, and replay at every input prime. "
        "This proves the old horizontal and three root sections; the q52 RR pencil, quartic, Jacobian and pointing remain separate gates."
        if all_passed else
        "The available CRT modulus has not reconstructed every section with exact QQ(t) substitution. This is progress only."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "model_path": str(MODEL.relative_to(ROOT)),
        "model_sha256": sha256(MODEL),
    },
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1CRT|primes={}|bits={}|sections={}/4|status={}|output={}".format(
        len(primes), payload["crt_modulus_bits"], sum(
            result["exact_QQ_weierstrass_identity"] for result in section_results
        ), payload["status"], output,
    ),
    flush=True,
)

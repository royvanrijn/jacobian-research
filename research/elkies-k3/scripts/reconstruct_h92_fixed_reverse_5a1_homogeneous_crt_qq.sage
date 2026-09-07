#!/usr/bin/env sage
"""Structured CRT reconstruction of the four q52 sections over QQ.

Reconstruct the homogeneous section coordinates separately:

    deg(X,Y,Z) = (2n+4, 3n+6, n),  Z monic.

The modular Abel compiler already supplies this normalization.  Removing the
redundant Z^2 and Z^3 denominator coefficients lowers the largest projective
LLL from dimension 98 to 53.  One small prime is withheld from CRT and filters
the LLL rows before the single exact QQ Weierstrass substitution.  No
Groebner basis or nonlinear characteristic-zero solve is used.
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
WITHHELD = LOCAL / "fixed-reverse-5a1-abel-word-seeds-mod167.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "fixed-reverse-5a1-homogeneous-crt-qq.json",
)
args = parser.parse_args()


def prime_suffix(path):
    suffix = path.stem.rsplit("mod", 1)[1]
    return int(suffix) if suffix.isdigit() else None


inputs = args.inputs or sorted(
    path for path in LOCAL.glob("fixed-reverse-5a1-abel-word-seeds-mod*.json")
    if prime_suffix(path) is not None and prime_suffix(path) > 1000000
)
inputs = [path if path.is_absolute() else ROOT / path for path in inputs]
output = args.output if args.output.is_absolute() else ROOT / args.output
started = time.monotonic()


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduce_qq(value, prime):
    value = QQ(value)
    return (
        int(value.numerator())
        * pow(int(value.denominator()), -1, int(prime))
    ) % int(prime)


assert inputs and WITHHELD.exists()
records = [read_json(path) for path in inputs]
withheld = read_json(WITHHELD)
primes = [ZZ(record["prime"]) for record in records]
withheld_prime = ZZ(withheld["prime"])
assert withheld_prime not in primes and len(set(primes)) == len(primes)
assert all(prime.is_prime() for prime in primes)
assert all(
    record["status"] == "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_REGULAR_SEEDS"
    for record in records + [withheld]
)
assert len({json.dumps(record["words"]) for record in records + [withheld]}) == 1
expected_poles = [15, 1, 1, 3]
assert all(
    [section["P_dot_O"] for section in record["sections"]] == expected_poles
    for record in records + [withheld]
)

modulus = math.prod(primes)
model = read_json(MODEL)
R = PolynomialRing(QQ, "t")
A = R([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])


def sizes(pole):
    return {"X": 2 * pole + 5, "Y": 3 * pole + 7, "Z": pole + 1}


def modular_values(record, section_index, coordinate, size):
    values = list(record["sections"][section_index][
        f"{coordinate}_coefficients_low_to_high"
    ])
    assert len(values) <= size
    values += [0] * (size - len(values))
    assert values[-1] % int(record["prime"])
    if coordinate == "Z":
        assert values[-1] % int(record["prime"]) == 1
    return [int(value) % int(record["prime"]) for value in values]


def projective_candidates(section_index, coordinate, size):
    residue_rows = [
        modular_values(record, section_index, coordinate, size)
        for record in records
    ]
    # Z is monic, so its last coefficient is already the affine anchor.
    unknown_size = size - 1 if coordinate == "Z" else size
    crt = [
        ZZ(CRT_list([row[index] for row in residue_rows], primes))
        for index in range(unknown_size)
    ]
    dimension = unknown_size + 1
    basis = matrix(ZZ, dimension, dimension)
    for index in range(unknown_size):
        basis[index, index] = modulus
    basis[-1] = vector(ZZ, crt + [1])
    reduced = basis.LLL(delta=0.99)
    withheld_values = modular_values(
        withheld, section_index, coordinate, size
    )
    candidates = []
    congruent_count = 0
    for row in sorted(reduced.rows(), key=lambda value: value * value):
        scale = ZZ(row[-1])
        if not scale or gcd(scale, modulus) != 1:
            continue
        if not all(
            (row[index] - scale * crt[index]) % modulus == 0
            for index in range(unknown_size)
        ):
            continue
        congruent_count += 1
        coefficients = [QQ(row[index]) / scale for index in range(unknown_size)]
        if coordinate == "Z":
            coefficients.append(QQ.one())
        if [reduce_qq(value, withheld_prime) for value in coefficients] != withheld_values:
            continue
        candidates.append({
            "polynomial": R(coefficients),
            "scale": scale,
            "primitive_max_bits": max(abs(ZZ(value)).nbits() for value in row),
        })
    return candidates, residue_rows, congruent_count


def verify_all_reductions(section_index, coordinate, polynomial, residue_rows):
    size = sizes(expected_poles[section_index])[coordinate]
    coefficients = list(polynomial) + [QQ.zero()] * (size - len(polynomial.list()))
    return all(
        [reduce_qq(value, prime) for value in coefficients] == expected
        for prime, expected in zip(primes, residue_rows)
    )


def polynomial_record(polynomial):
    return {
        "coefficients_low_to_high": list(map(str, polynomial.list())),
        "degree": int(polynomial.degree()),
    }


section_results = []
all_passed = True
for section_index, pole in enumerate(expected_poles):
    coordinate_candidates = {}
    residue_rows = {}
    congruent_counts = {}
    for coordinate, size in sizes(pole).items():
        candidates, residues, congruent_count = projective_candidates(
            section_index, coordinate, size
        )
        coordinate_candidates[coordinate] = candidates
        residue_rows[coordinate] = residues
        congruent_counts[coordinate] = congruent_count
    selected = None
    triples_tested = 0
    for Z_candidate in coordinate_candidates["Z"]:
        Z = Z_candidate["polynomial"]
        if Z.degree() != pole or not Z.is_monic():
            continue
        for X_candidate in coordinate_candidates["X"]:
            X = X_candidate["polynomial"]
            if X.degree() != 2 * pole + 4:
                continue
            for Y_candidate in coordinate_candidates["Y"]:
                Y = Y_candidate["polynomial"]
                triples_tested += 1
                if Y.degree() != 3 * pole + 6:
                    continue
                if Y ** 2 != X ** 3 + A * X * Z ** 4 + B * Z ** 6:
                    continue
                if not all(
                    verify_all_reductions(
                        section_index, coordinate, candidate["polynomial"],
                        residue_rows[coordinate],
                    )
                    for coordinate, candidate in (
                        ("X", X_candidate), ("Y", Y_candidate), ("Z", Z_candidate)
                    )
                ):
                    continue
                selected = {
                    "X": X_candidate, "Y": Y_candidate, "Z": Z_candidate
                }
                break
            if selected is not None:
                break
        if selected is not None:
            break
    passed = selected is not None
    all_passed = all_passed and passed
    result = {
        "section_index": section_index,
        "P_dot_O": pole,
        "candidate_counts_X_Y_Z": [
            len(coordinate_candidates[name]) for name in ("X", "Y", "Z")
        ],
        "congruent_LLL_row_counts_X_Y_Z": [
            congruent_counts[name] for name in ("X", "Y", "Z")
        ],
        "triples_tested": triples_tested,
        "exact_QQ_weierstrass_identity": passed,
        "reduction_to_CRT_and_withheld_primes": passed,
    }
    if passed:
        X = selected["X"]["polynomial"]
        Y = selected["Y"]["polynomial"]
        Z = selected["Z"]["polynomial"]
        result.update({
            "X": polynomial_record(X),
            "Y": polynomial_record(Y),
            "Z": polynomial_record(Z),
            "projective_scales_X_Y_Z": [
                str(selected[name]["scale"]) for name in ("X", "Y", "Z")
            ],
            "primitive_max_bits_X_Y_Z": [
                selected[name]["primitive_max_bits"] for name in ("X", "Y", "Z")
            ],
            "x": {
                "numerator_coefficients_low_to_high": list(map(str, X.list())),
                "denominator_coefficients_low_to_high": list(map(str, (Z ** 2).list())),
            },
            "y": {
                "numerator_coefficients_low_to_high": list(map(str, Y.list())),
                "denominator_coefficients_low_to_high": list(map(str, (Z ** 3).list())),
            },
        })
    section_results.append(result)
    print(
        "FIXEDREVERSE5A1HOMCRT|section={}|pole={}|candidates={}|tested={}|passed={}|seconds={:.3f}".format(
            section_index, pole, result["candidate_counts_X_Y_Z"],
            triples_tested, int(passed), time.monotonic() - started,
        ),
        flush=True,
    )

payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-homogeneous-crt-qq.v1",
    "status": (
        "PASS_EXACT_QQ_FIXED_REVERSE_5A1_SECTIONS_HOMOGENEOUS_CRT"
        if all_passed else "PARTIAL_FIXED_REVERSE_5A1_SECTIONS_HOMOGENEOUS_CRT"
    ),
    "primes": list(map(int, primes)),
    "withheld_prime": int(withheld_prime),
    "crt_modulus_bits": int(ZZ(modulus).nbits()),
    "words": records[0]["words"],
    "sections": section_results,
    "method": {
        "homogeneous_X_Y_Z_reconstruction": True,
        "largest_LLL_dimension": 3 * expected_poles[0] + 8,
        "withheld_prime_filter_before_QQ_substitution": True,
        "groebner_or_nonlinear_QQ_solve": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Exact homogeneous projective CRT reconstruction, literal QQ(t) substitution, and replay at every CRT and withheld prime. The q52 RR pencil, quartic, Jacobian and pointing remain separate gates."
        if all_passed else
        "The current CRT modulus does not yet reconstruct every homogeneous section coordinate. This is reconstruction progress only."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "withheld_path": str(WITHHELD.relative_to(ROOT)),
        "withheld_sha256": sha256(WITHHELD),
        "model_path": str(MODEL.relative_to(ROOT)),
        "model_sha256": sha256(MODEL),
    },
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1HOMCRT|primes={}|bits={}|sections={}/4|status={}|output={}".format(
        len(primes), payload["crt_modulus_bits"], sum(
            result["exact_QQ_weierstrass_identity"] for result in section_results
        ), payload["status"], output,
    ),
    flush=True,
)

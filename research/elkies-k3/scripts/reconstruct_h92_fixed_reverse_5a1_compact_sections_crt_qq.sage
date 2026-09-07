#!/usr/bin/env sage
"""CRT-reconstruct the four q52 sections on the compact 4A1 model.

Transform every modular Abel-word section through the exact cross-ratio base
change and QQ Weierstrass isomorphism before CRT.  In this normalization the
surface coefficients are at most 215 bits, so ordinary coefficientwise
rational reconstruction replaces the failed million-bit Hensel lift.  Exact
QQ(t) substitution and one withheld-prime replay are mandatory.
"""

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from sage.all import CRT_list, GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT = LOCAL / "fixed-reverse-4a1-compact-crossratio-qq.json"
WITHHELD = LOCAL / "fixed-reverse-5a1-abel-word-seeds-mod167.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "fixed-reverse-5a1-compact-sections-crt-qq.json",
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


def pad(values, size):
    assert len(values) <= size
    return list(values) + [0] * (size - len(values))


compact = read_json(COMPACT)
records = [read_json(path) for path in inputs]
withheld = read_json(WITHHELD)
primes = [ZZ(record["prime"]) for record in records]
assert compact["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_COMPACT_CROSSRATIO"
assert inputs and len(set(primes)) == len(primes)
assert all(prime.is_prime() for prime in primes)
assert ZZ(withheld["prime"]) not in primes
assert all(
    record["status"] == "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_REGULAR_SEEDS"
    for record in records + [withheld]
)
assert len({json.dumps(record["words"]) for record in records + [withheld]}) == 1
expected_poles = [15, 1, 1, 3]

modulus = math.prod(primes)
RQ = PolynomialRing(QQ, "T")
A = RQ(compact["compact_model"]["A_coefficients_low_to_high"])
B = RQ(compact["compact_model"]["B_coefficients_low_to_high"])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


def transformed_section(record, section_index):
    p = ZZ(record["prime"])
    F = GF(p)
    R = PolynomialRing(F, "T")
    K = R.fraction_field()
    base_numerator = R([
        reduce_qq(value, F) for value in
        compact["base_change"]["base_numerator_coefficients_low_to_high"]
    ])
    base_denominator = R([
        reduce_qq(value, F) for value in
        compact["base_change"]["base_denominator_coefficients_low_to_high"]
    ])
    u = reduce_qq(compact["weierstrass_isomorphism"]["u"], F)
    base = K(base_numerator) / K(base_denominator)

    def evaluate(values):
        answer = K.zero()
        for coefficient in reversed(values):
            answer = answer * base + F(coefficient)
        return answer

    section = record["sections"][section_index]
    x_old = (
        evaluate(section["x_numerator_coefficients_low_to_high"])
        / evaluate(section["x_denominator_coefficients_low_to_high"])
    )
    y_old = (
        evaluate(section["y_numerator_coefficients_low_to_high"])
        / evaluate(section["y_denominator_coefficients_low_to_high"])
    )
    x = K(base_denominator ** 4) * x_old / F(u ** 2)
    y = K(base_denominator ** 6) * y_old / F(u ** 3)
    A_mod = R([reduce_qq(value, F) for value in A])
    B_mod = R([reduce_qq(value, F) for value in B])
    assert y ** 2 == x ** 3 + K(A_mod) * x + K(B_mod)

    def normalized(value):
        numerator = R(value.numerator())
        denominator = R(value.denominator())
        scale = denominator.leading_coefficient() ** -1
        numerator *= scale
        denominator *= scale
        assert denominator.is_monic()
        return {
            "numerator": list(map(int, numerator.list())),
            "denominator": list(map(int, denominator.list())),
            "degrees": [int(numerator.degree()), int(denominator.degree())],
        }

    return {"x": normalized(x), "y": normalized(y)}


transformed = [
    [transformed_section(record, index) for index in range(4)]
    for record in records
]
withheld_transformed = [transformed_section(withheld, index) for index in range(4)]
degree_fingerprints = [
    [transformed[0][index][name]["degrees"] for name in ("x", "y")]
    for index in range(4)
]
assert all(
    [
        [row[index][name]["degrees"] for name in ("x", "y")]
        for index in range(4)
    ] == degree_fingerprints
    for row in transformed
)


def reconstruct_part(section_index, coordinate, part):
    size = max(
        len(row[section_index][coordinate][part]) for row in transformed
    )
    coefficients = []
    failures = []
    for coefficient_index in range(size):
        residues = [
            pad(row[section_index][coordinate][part], size)[coefficient_index]
            for row in transformed
        ]
        residue = ZZ(CRT_list(residues, primes))
        try:
            coefficients.append(QQ(residue.rational_reconstruction(modulus)))
        except (ArithmeticError, ValueError) as error:
            coefficients.append(None)
            failures.append({"index": coefficient_index, "reason": str(error)})
    return coefficients, failures


def reduce_polynomial(poly, prime):
    return [
        int(
            (coefficient.numerator()
             * pow(int(coefficient.denominator()), -1, int(prime)))
            % int(prime)
        )
        for coefficient in poly
    ]


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
    reconstructed = {}
    failures = []
    for coordinate in ("x", "y"):
        reconstructed[coordinate] = {}
        for part in ("numerator", "denominator"):
            values, part_failures = reconstruct_part(
                section_index, coordinate, part
            )
            reconstructed[coordinate][part] = values
            failures.extend({
                "coordinate": coordinate, "part": part, **failure
            } for failure in part_failures)
    passed = not failures
    x = y = None
    withheld_verified = False
    if passed:
        x = (
            RQ(reconstructed["x"]["numerator"])
            / RQ(reconstructed["x"]["denominator"])
        )
        y = (
            RQ(reconstructed["y"]["numerator"])
            / RQ(reconstructed["y"]["denominator"])
        )
        passed = y ** 2 == x ** 3 + A * x + B
        if passed:
            withheld_prime = ZZ(withheld["prime"])
            withheld_verified = all(
                reduce_polynomial(getattr(value, part)(), withheld_prime)
                == withheld_transformed[section_index][coordinate][part]
                for coordinate, value in (("x", x), ("y", y))
                for part in ("numerator", "denominator")
            )
            passed = withheld_verified
    all_passed = all_passed and passed
    successful = [
        value for coordinate in reconstructed.values()
        for values in coordinate.values() for value in values if value is not None
    ]
    result = {
        "section_index": section_index,
        "P_dot_O": pole,
        "compact_degree_fingerprint_x_y": degree_fingerprints[section_index],
        "rational_reconstruction_failures": failures,
        "successful_coefficient_count": len(successful),
        "exact_QQ_weierstrass_identity": bool(passed),
        "withheld_prime_replay": bool(withheld_verified),
    }
    if passed:
        result.update({
            "x": qq_record(x),
            "y": qq_record(y),
            "maximum_rational_bits": max(
                max(abs(value.numerator()).nbits(), value.denominator().nbits())
                for value in successful
            ),
        })
    section_results.append(result)
    print(
        "FIXEDREVERSE5A1COMPACTCRT|section={}|pole={}|degrees={}|failures={}|passed={}|seconds={:.3f}".format(
            section_index, pole, degree_fingerprints[section_index],
            len(failures), int(passed), time.monotonic() - started,
        ),
        flush=True,
    )

payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-compact-sections-crt-qq.v1",
    "status": (
        "PASS_EXACT_QQ_FIXED_REVERSE_5A1_COMPACT_SECTIONS_CRT"
        if all_passed else "PARTIAL_FIXED_REVERSE_5A1_COMPACT_SECTIONS_CRT"
    ),
    "primes": list(map(int, primes)),
    "withheld_prime": int(withheld["prime"]),
    "crt_modulus_bits": int(ZZ(modulus).nbits()),
    "words": records[0]["words"],
    "sections": section_results,
    "method": {
        "exact_crossratio_transform_before_CRT": True,
        "coefficientwise_rational_reconstruction": True,
        "groebner_or_nonlinear_QQ_solve": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Exact compact-model CRT reconstruction, literal QQ(t) substitution, and withheld-prime replay. The q52 RR pencil, quartic, Jacobian and pointing remain separate gates."
        if all_passed else
        "The current CRT modulus does not yet reconstruct every compact section coordinate. This is progress only."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "compact_path": str(COMPACT.relative_to(ROOT)),
        "compact_sha256": sha256(COMPACT),
        "withheld_path": str(WITHHELD.relative_to(ROOT)),
        "withheld_sha256": sha256(WITHHELD),
    },
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1COMPACTCRT|primes={}|bits={}|sections={}/4|status={}|output={}".format(
        len(primes), payload["crt_modulus_bits"], sum(
            result["exact_QQ_weierstrass_identity"] for result in section_results
        ), payload["status"], output,
    ),
    flush=True,
)

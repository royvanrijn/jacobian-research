#!/usr/bin/env sage -python
"""Interpolate the generic p=19 third-q12 Jacobian from mapped fibres."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-weierstrass-sample-batch.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"
EXPECTED = "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_MOD19_QUADRATIC"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(INPUT.read_text())
if payload.get("status") != EXPECTED:
    raise ValueError("mapped sample batch is not certified")

base_finite = GF(19)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)
r = finite.gen()
polynomial_ring = PolynomialRing(finite, "V")
V = polynomial_ring.gen()
function_field = polynomial_ring.fraction_field()


def element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


def sample_value(record):
    return element(record["new_base_coefficients_1_r"])


training = payload["training_samples"]
held_out = payload["held_out_samples"]
if len(training) < 62 or len(held_out) != 8:
    raise ArithmeticError("unexpected training/holdout split")


def interpolate(records, accessor, label, max_total_degree=60):
    samples = [(sample_value(record), element(accessor(record))) for record in records]
    for total_degree in range(max_total_degree + 1):
        for numerator_degree in range(total_degree + 1):
            denominator_degree = total_degree - numerator_degree
            column_count = numerator_degree + denominator_degree + 2
            if column_count > len(samples):
                continue
            rows = []
            for value, target in samples:
                rows.append(
                    [value**index for index in range(numerator_degree + 1)]
                    + [
                        -target * value**index
                        for index in range(denominator_degree + 1)
                    ]
                )
            kernel = Matrix(finite, rows).right_kernel()
            if kernel.dimension() != 1:
                continue
            relation = kernel.basis()[0]
            numerator = polynomial_ring(list(relation[: numerator_degree + 1]))
            denominator = polynomial_ring(list(relation[numerator_degree + 1 :]))
            if not denominator:
                continue
            common = numerator.gcd(denominator)
            numerator //= common
            denominator //= common
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            if any(
                not denominator(value) or numerator(value) != target * denominator(value)
                for value, target in samples
            ):
                continue
            candidate = function_field(numerator / denominator)
            for record in held_out:
                value = sample_value(record)
                target = element(accessor(record))
                if not denominator(value) or candidate(value) != target:
                    break
            else:
                return {
                    "label": label,
                    "numerator": numerator,
                    "denominator": denominator,
                    "numerator_degree": int(numerator.degree()),
                    "denominator_degree": int(denominator.degree()),
                    "function": candidate,
                }
    raise ArithmeticError(f"no held-out-valid interpolant for {label}")


coefficient_names = ("a1", "a2", "a3", "a4", "a6")
coefficient_records = []
for index, name in enumerate(coefficient_names):
    coefficient_records.append(
        interpolate(
            training,
            lambda record, index=index: record["a1_a2_a3_a4_a6"][index],
            name,
        )
    )

delta_record = interpolate(training, lambda record: record["discriminant"], "Delta")
j_record = interpolate(training, lambda record: record["j"], "j")
coefficients = [record["function"] for record in coefficient_records]
curve = EllipticCurve(function_field, coefficients)
if curve.discriminant() != delta_record["function"]:
    raise ArithmeticError("interpolated long model has the wrong discriminant")
if curve.j_invariant() != j_record["function"]:
    raise ArithmeticError("interpolated long model has the wrong j-invariant")

# Replay every sample, not only the held-out subset.
for record in training + held_out:
    value = sample_value(record)
    if any(
        coefficient(value) != element(record["a1_a2_a3_a4_a6"][index])
        for index, coefficient in enumerate(coefficients)
    ):
        raise ArithmeticError("generic coefficient replay failed")
    if curve.discriminant()(value) != element(record["discriminant"]):
        raise ArithmeticError("generic discriminant replay failed")
    if curve.j_invariant()(value) != element(record["j"]):
        raise ArithmeticError("generic j replay failed")


def polynomial_record(poly):
    return [
        [int((list(value.list()) + [base_finite.zero(), base_finite.zero()])[0]),
         int((list(value.list()) + [base_finite.zero(), base_finite.zero()])[1])]
        for value in poly.list()
    ]


def rational_record(record):
    return {
        "numerator_coefficients_low_to_high_1_r": polynomial_record(record["numerator"]),
        "denominator_coefficients_low_to_high_1_r": polynomial_record(record["denominator"]),
        "degrees_numerator_denominator": [
            record["numerator_degree"],
            record["denominator_degree"],
        ],
    }


delta_numerator = curve.discriminant().numerator()
delta_denominator = curve.discriminant().denominator()
output = {
    "schema": "elkies-k3.q80-third-q12-jacobian-interpolated-modp2.v1",
    "status": "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"},
    "gauge": (
        "new base V is the resolved kernel ratio; at t=1/W,xi=-6, "
        "x=t^-2+... with constant 0 and y=t^-3+0*t^-2+... with constant 0"
    ),
    "weierstrass": {
        name: rational_record(record)
        for name, record in zip(coefficient_names, coefficient_records)
    },
    "discriminant": {
        **rational_record(delta_record),
        "numerator_factorization": [
            [str(factor), int(exponent)]
            for factor, exponent in delta_numerator.factor()
        ],
        "denominator_factorization": [
            [str(factor), int(exponent)]
            for factor, exponent in delta_denominator.factor()
        ],
    },
    "j": rational_record(j_record),
    "validation": {
        "training_samples": len(training),
        "held_out_samples": len(held_out),
        "all_samples_replayed": len(training) + len(held_out),
    },
    "input": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "generic long Weierstrass Jacobian in the pinned Laurent gauge over GF(19^2)(V)",
            "exact discriminant and j identities",
            f"{len(training)} training and {len(held_out)} held-out mapped-fibre replays",
        ],
        "not_proved": [
            "generic interpolation of the forward/inverse maps",
            "global minimality and A5+A3+3A1 fibre marking",
            "Frobenius-invariant encoding or a second prime",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/interpolate_q80_third_q12_jacobian_mod19_quadratic.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12JACOBIAN|training={len(training)}|heldout={len(held_out)}|"
    f"replayed={len(training) + len(held_out)}|"
    f"degrees={tuple((record['numerator_degree'], record['denominator_degree']) for record in coefficient_records)}|"
    f"Delta_degrees={delta_record['numerator_degree']},{delta_record['denominator_degree']}|"
    f"j_degrees={j_record['numerator_degree']},{j_record['denominator_degree']}|"
    "status=PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC",
    flush=True,
)

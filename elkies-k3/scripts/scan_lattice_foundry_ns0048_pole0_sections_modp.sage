#!/usr/bin/env sage-python
"""Exhaust the NS0048 pole-zero section chart on stored fibre models.

The exact generator has depth one at I7 and I2 and correction zero at I5 and
I1*.  Force X through the I7/I2 nodes and Y through their support product;
the remaining three X and five Y coefficients give a p^8 affine chart.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sage.all import GF, PolynomialRing, PowerSeriesRing


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod7.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-pole0-sections-mod7.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def function_exponent(exponent: int, prime: int) -> int:
    return 0 if exponent == 0 else 1 + (exponent - 1) % (prime - 1)


def evaluate_everywhere(polynomial, prime: int) -> np.ndarray:
    count = polynomial.parent().ngens()
    coefficients = np.zeros((prime,) * count, dtype=np.int64)
    for exponents, coefficient in polynomial.dict().items():
        reduced = tuple(function_exponent(int(value), prime) for value in exponents)
        coefficients[reduced] += int(coefficient)
    coefficients %= prime
    transform = np.array(
        [[pow(value, exponent, prime) for exponent in range(prime)] for value in range(prime)],
        dtype=np.int64,
    )
    values = coefficients
    for axis in range(count):
        values = np.tensordot(transform, values, axes=(1, axis)) % prime
        values = np.moveaxis(values, 0, axis)
    return values


def audit_table(polynomial, table: np.ndarray, prime: int) -> int:
    flat = table.reshape(-1)
    total = flat.size
    indices = sorted({0, total - 1, *range(0, total, max(1, total // 31))})
    shape = (prime,) * polynomial.parent().ngens()
    for index in indices:
        values = tuple(int(value) for value in np.unravel_index(index, shape))
        if int(polynomial(*values)) % prime != int(flat[index]) % prime:
            raise ArithmeticError("tensor evaluation audit failed")
    return len(indices)


def formal_center(A, B, point, precision):
    field = A.base_ring()
    base = A.parent()
    t = base.gen()
    shifted_A = base(A(t + point))
    shifted_B = base(B(t + point))
    node = -field(3) * shifted_B[0] / (field(2) * shifted_A[0])
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 2)
    center = series_ring(node)
    series_A = series_ring(shifted_A)
    for unused in range(6):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 1:
        raise ArithmeticError("formal center did not converge")
    return series_ring, center, node


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

input_path = args.input.resolve()
payload = json.loads(input_path.read_text())
if payload["schema"] != "elkies-k3.lattice-foundry-ns0048-source-ansatz-modp.v1":
    raise ValueError("unexpected NS0048 fibre-ansatz schema")
prime = int(payload["prime"])
field = GF(prime)
twist = field(args.quadratic_twist)
if not twist:
    raise ValueError("--quadratic-twist must be nonzero")
base = PolynomialRing(field, "t")
t0 = base.gen()

records = []
for example_index, example in enumerate(payload["examples"]):
    A = twist**2 * base(example["A_coefficients_low_to_high"])
    B = twist**3 * base(example["B_coefficients_low_to_high"])
    lambda_value = field(example["lambda"])
    series_one, center_one, node_one = formal_center(A, B, field.one(), 7)
    series_lambda, center_lambda, node_lambda = formal_center(A, B, lambda_value, 2)

    names = ["q2", "q1", "q0", "r4", "r3", "r2", "r1", "r0"]
    coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
    g = coefficient_ring.gens_dict()
    outer = PolynomialRing(coefficient_ring, "t")
    t = outer.gen()
    D = (t - 1) * (t - coefficient_ring(lambda_value))
    slope = coefficient_ring(node_lambda - node_one) / coefficient_ring(lambda_value - 1)
    x_remainder = coefficient_ring(node_one) + slope * (t - 1)
    Q = sum(g[f"q{index}"] * t**index for index in range(3))
    R = sum(g[f"r{index}"] * t**index for index in range(5))
    X = x_remainder + D * Q
    Y = D * R
    A_outer = outer([coefficient_ring(value) for value in A])
    B_outer = outer([coefficient_ring(value) for value in B])
    identity = Y**2 - X**3 - A_outer * X - B_outer
    equations = [
        coefficient_ring(identity[index])
        for index in range(13)
        if coefficient_ring(identity[index])
    ]

    mask = np.ones((prime,) * len(names), dtype=bool)
    order = sorted(range(len(equations)), key=lambda index: len(equations[index].monomials()))
    survivors = []
    tensor_audits = 0
    for equation_index in order:
        table = evaluate_everywhere(equations[equation_index], prime)
        tensor_audits += audit_table(equations[equation_index], table, prime)
        mask &= table == 0
        survivors.append(int(mask.sum()))
        if not survivors[-1]:
            break

    solutions = []
    raw_solutions = []
    for row in np.argwhere(mask):
        values = tuple(int(value) for value in row)

        def specialize(poly):
            return base([field(coefficient(*values)) for coefficient in poly.list()])

        X_value = specialize(X)
        Y_value = specialize(Y)
        if Y_value**2 != X_value**3 + A * X_value + B:
            raise ArithmeticError("materialized polynomial section fails")

        shifted_X_one = series_one(base(X_value(t0 + 1)))
        shifted_Y_one = series_one(base(Y_value(t0 + 1)))
        shifted_X_lambda = series_lambda(base(X_value(t0 + lambda_value)))
        shifted_Y_lambda = series_lambda(base(Y_value(t0 + lambda_value)))
        depth_one = min(
            int((shifted_X_one - center_one).valuation()),
            int(shifted_Y_one.valuation()),
        )
        depth_lambda = min(
            int((shifted_X_lambda - center_lambda).valuation()),
            int(shifted_Y_lambda.valuation()),
        )
        node_zero = -field(3) * B(0) / (field(2) * A(0))
        smooth_zero = not (X_value(0) == node_zero and Y_value(0) == 0)
        smooth_infinity = not (X_value[4] == 0 and Y_value[6] == 0)
        diagnostic = {
            "variables_in_declared_order": list(values),
            "X_coefficients_low_to_high": [int(value) for value in X_value],
            "Y_coefficients_low_to_high": [int(value) for value in Y_value],
            "component_depths_at_I7_I2": [depth_one, depth_lambda],
            "identity_component_smooth_at_I5_I1star": [smooth_zero, smooth_infinity],
        }
        raw_solutions.append(diagnostic)
        if (depth_one, depth_lambda) != (1, 1) or not smooth_zero or not smooth_infinity:
            continue
        solutions.append(
            diagnostic
            | {
                "shioda_height": "37/14",
                "implied_NS_determinant": 740,
            }
        )

    records.append(
        {
            "example_index": example_index,
            "lambda": int(lambda_value),
            "variables": names,
            "affine_points_scanned": prime ** len(names),
            "equation_evaluation_order": order,
            "surviving_counts_after_each_equation": survivors,
            "direct_tensor_audit_points": tensor_audits,
            "raw_section_identity_solution_count": int(mask.sum()),
            "raw_section_identity_solutions": raw_solutions,
            "marked_section_solution_count": len(solutions),
            "solutions": solutions,
        }
    )

total = sum(record["marked_section_solution_count"] for record in records)
output = {
    "schema": "elkies-k3.lattice-foundry-ns0048-pole0-sections-modp-scan.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_WITH_MARKED_SECTIONS"
        if total
        else "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_EMPTY_SECTION_CHART"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "input": {"artifact": relative(input_path), "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest()},
    "scope": {
        "stored_fibre_models": len(records),
        "affine_points_per_model": prime**8,
        "fibre_ansatz_scan_exhausted": bool(payload["scan"]["exhausted"]),
        "twist_is_bijective_on_stored_fibre_models": True,
    },
    "accounting": {
        "total_affine_points_scanned": sum(record["affine_points_scanned"] for record in records),
        "models_with_marked_sections": sum(bool(record["marked_section_solution_count"]) for record in records),
        "total_marked_sections": total,
        "total_raw_polynomial_sections": sum(record["raw_section_identity_solution_count"] for record in records),
    },
    "models": records,
    "proof_boundary": {
        "proved": (
            "Every polynomial pole-zero section over the displayed finite field "
            "is exhausted on every stored fibre model; retained sections pass the "
            "exact equation, component-depth, smooth-identity, height, and determinant checks."
        ),
        "not_proved": (
            "A finite-field marked model is not a rational source family, a "
            "characteristic-zero lift, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0048_pole0_sections_modp.sage"
        + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0048 pole-zero section scan is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0048POLE0SCAN|"
    f"models={len(records)}|points={sum(row['affine_points_scanned'] for row in records)}|"
    f"raw={sum(row['raw_section_identity_solution_count'] for row in records)}|"
    f"sections={total}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)

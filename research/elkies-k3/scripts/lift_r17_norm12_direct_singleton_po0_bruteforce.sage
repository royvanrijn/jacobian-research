#!/usr/bin/env sage-python
"""Lift isolated finite-field sections of a direct norm-12 singleton twist.

The input is an exhaustive ``chi=3`` brute-force shell produced by
``run_twist_polynomial_sections_bruteforce.py``.  For each requested solution
the script restores the full eight-variable coefficient system, checks the
first Hensel obstruction using every residual equation, Newton-lifts an
isolated compatible point, rationally reconstructs it, transforms it back to
the original base chart, and verifies the literal characteristic-zero twist
equation.  Failed compatibility is an exact local obstruction; failed rational
reconstruction is only a bounded negative result.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import lcm
from pathlib import Path
import sys

from sage.all import GF, Integers, QQ, ZZ, PolynomialRing, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from screen_elkies_2026_quadratic_twist_ranks import (  # noqa: E402
    square_equivalent_integer_polynomial,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def convolution(left, right):
    zero = left[0].parent()(0)
    result = [zero] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return result


def coefficient(values, index):
    return values[index] if index < len(values) else values[0].parent()(0)


def residual_and_jacobian(z, coefficient_a, coefficient_b):
    """Return the complete reduced residual system, its Jacobian, and Y."""

    x_degree = (len(z) - 2)
    y_degree = 3 * x_degree // 2
    if 2 * y_degree != 3 * x_degree:
        raise ArithmeticError("incompatible section degree bounds")
    X = list(z[: x_degree + 1])
    leading_y = z[-1]
    zero = leading_y.parent()(0)
    variable_count = x_degree + 2

    X2 = convolution(X, X)
    rhs = convolution(X2, X)
    AX = convolution(coefficient_a, X)
    rhs += [zero] * (2 * y_degree + 1 - len(rhs))
    for degree, value in enumerate(AX):
        rhs[degree] += value
    for degree, value in enumerate(coefficient_b):
        rhs[degree] += value

    kernel = [
        3 * coefficient(X2, degree) + coefficient(coefficient_a, degree)
        for degree in range(max(len(X2), len(coefficient_a)))
    ]
    derivative_rhs = []
    for variable in range(x_degree + 1):
        row = [zero] * (2 * y_degree + 1)
        for degree, value in enumerate(kernel):
            if degree + variable < len(row):
                row[degree + variable] += value
        derivative_rhs.append(row)
    derivative_rhs.append([zero] * (2 * y_degree + 1))

    Y = [zero] * (y_degree + 1)
    derivative_y = [[zero] * variable_count for unused in range(y_degree + 1)]
    Y[y_degree] = leading_y
    derivative_y[y_degree][-1] = leading_y.parent()(1)
    denominator = 2 * leading_y
    if not denominator:
        raise ZeroDivisionError("the leading Y coefficient vanishes")
    for degree in range(2 * y_degree - 1, y_degree - 1, -1):
        index = degree - y_degree
        known = zero
        for left in range(index + 1, y_degree + 1):
            right = degree - left
            if index < right <= y_degree:
                known += Y[left] * Y[right]
        Y[index] = (rhs[degree] - known) / denominator
        for variable in range(variable_count):
            derivative_known = zero
            for left in range(index + 1, y_degree + 1):
                right = degree - left
                if index < right <= y_degree:
                    derivative_known += (
                        derivative_y[left][variable] * Y[right]
                        + Y[left] * derivative_y[right][variable]
                    )
            derivative_denominator = 2 if variable == variable_count - 1 else 0
            derivative_y[index][variable] = (
                derivative_rhs[variable][degree]
                - derivative_known
                - Y[index] * derivative_denominator
            ) / denominator

    Y2 = convolution(Y, Y)
    residual_degrees = [2 * y_degree, *range(y_degree)]
    residual = [Y2[degree] - rhs[degree] for degree in residual_degrees]
    jacobian = []
    for degree in residual_degrees:
        row = []
        for variable in range(variable_count):
            derivative = zero
            for left in range(y_degree + 1):
                right = degree - left
                if 0 <= right <= y_degree:
                    derivative += (
                        derivative_y[left][variable] * Y[right]
                        + Y[left] * derivative_y[right][variable]
                    )
            row.append(derivative - derivative_rhs[variable][degree])
        jacobian.append(row)
    return residual, jacobian, Y


def require_polynomial(value, ring, label):
    value = ring.fraction_field()(value)
    if value.denominator().degree() != 0:
        raise ArithmeticError(f"{label} is not polynomial in the original chart")
    return ring(value)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--bruteforce", type=Path, required=True)
parser.add_argument("--export", type=Path)
parser.add_argument("--bisections", type=Path, required=True)
parser.add_argument("--model", type=Path, required=True)
parser.add_argument("--solution-index", type=int, action="append", default=[])
parser.add_argument("--hensel-depth", type=int, default=128)
parser.add_argument("--output", type=Path, required=True)
arguments = parser.parse_args()
if arguments.hensel_depth < 2:
    parser.error("--hensel-depth must be at least two")

bruteforce_path = arguments.bruteforce.resolve()
bruteforce = json.loads(bruteforce_path.read_text())
if bruteforce.get("schema") != "elkies-k3.twist-polynomial-section-bruteforce.v1":
    raise ValueError("unexpected brute-force artifact schema")
if bruteforce["candidate"].get("kind") != "direct_singleton":
    raise ValueError("expected a direct singleton candidate")
if int(bruteforce["candidate"].get("chi", 0)) != 3:
    raise ValueError("this lift is for chi=3 singleton twists")

if arguments.export is None:
    export_inputs = [
        ROOT / path
        for path in bruteforce.get("inputs", {})
        if path.endswith("/export.json")
    ]
    if len(export_inputs) != 1:
        raise ValueError("cannot infer a unique modular export")
    export_path = export_inputs[0]
else:
    export_path = arguments.export.resolve()
export = json.loads(export_path.read_text())
if digest(export_path) != bruteforce["inputs"].get(relative(export_path)):
    raise ArithmeticError("brute-force artifact does not hash-pin this export")
if export["candidate"] != bruteforce["candidate"]:
    raise ArithmeticError("candidate mismatch between export and brute-force shell")

prime = int(export["prime"])
field = GF(prime)
label = str(export["candidate"]["key"])
solutions = bruteforce["solutions"]
selected_indices = arguments.solution_index or [
    index
    for index, solution in enumerate(solutions)
    if int(solution["full_shell_tangent_rank"]) == 8
]
if not selected_indices:
    raise ValueError("no full-rank solution was selected")
if any(index < 0 or index >= len(solutions) for index in selected_indices):
    raise ValueError("solution index is out of range")

bisections_path = arguments.bisections.resolve()
bisections = json.loads(bisections_path.read_text())
record = next(item for item in bisections["bisections"] if item["label"] == label)
q_coefficients = square_equivalent_integer_polynomial(
    record["branch"]["numerator_coefficients"]
)
model_path = arguments.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct norm-12 model")
model_coefficients = model["weierstrass_model"]

ring = PolynomialRing(QQ, "t")
t = ring.gen()
A0 = ring([QQ(value) for value in model_coefficients["A_coefficients_low_to_high"]])
B0 = ring([QQ(value) for value in model_coefficients["B_coefficients_low_to_high"]])
q = ring([QQ(value) for value in q_coefficients])
q_denominator = lcm(*(value.denominator() for value in q.list()))
q_integral = ring(q_denominator**2 * q)
A_twist = A0 * q_integral**2
B_twist = B0 * q_integral**3

x_degree = 6
y_degree = 9
chart = export["infinity_fibre"]["chart"]
chart_parameter = None
if chart == "original_infinity":
    A_chart = A_twist
    B_chart = B_twist
elif chart.startswith("original_t=") and " via " in chart:
    chart_parameter = QQ(int(chart[len("original_t=") :].split(" via ", 1)[0]))
    A_chart = ring(
        sum(
            A_twist[index]
            * (chart_parameter * t + 1) ** index
            * t ** (2 * x_degree - index)
            for index in range(A_twist.degree() + 1)
        )
    )
    B_chart = ring(
        sum(
            B_twist[index]
            * (chart_parameter * t + 1) ** index
            * t ** (2 * y_degree - index)
            for index in range(B_twist.degree() + 1)
        )
    )
else:
    raise ValueError("unrecognized exporter chart")

if [int(field(value)) for value in A_chart.list()] != export[
    "twist_A_coefficients_low_to_high"
] or [int(field(value)) for value in B_chart.list()] != export[
    "twist_B_coefficients_low_to_high"
]:
    raise ArithmeticError("exact twist chart does not reduce to the export")

lift_records = []
exact_sections = []
for solution_index in selected_indices:
    solution = solutions[solution_index]
    seed_x = [int(value) for value in solution["X_coefficients_low_to_high"]]
    seed_y = [int(value) for value in solution["Y_coefficients_low_to_high"]]
    if len(seed_x) != x_degree + 1 or len(seed_y) != y_degree + 1:
        raise ArithmeticError("brute-force solution has wrong coefficient lengths")
    seed = [*seed_x, seed_y[-1]]

    finite_residual, finite_jacobian, recovered_y = residual_and_jacobian(
        vector(field, seed),
        [field(value) for value in export["twist_A_coefficients_low_to_high"]],
        [field(value) for value in export["twist_B_coefficients_low_to_high"]],
    )
    if any(finite_residual) or [int(value) for value in recovered_y] != seed_y:
        raise ArithmeticError("brute-force seed does not replay the finite-field system")
    finite_matrix = matrix(field, finite_jacobian)
    tangent_rank = int(finite_matrix.rank())

    rational_residual, rational_jacobian, unused_y = residual_and_jacobian(
        vector(QQ, seed), list(A_chart), list(B_chart)
    )
    if any(value.valuation(prime) < 1 for value in rational_residual if value):
        raise ArithmeticError("seed residual is not divisible by the lift prime")
    first_order_target = vector(
        field, [field(-value / prime) for value in rational_residual]
    )
    augmented_rank = int(
        finite_matrix.augment(first_order_target.column()).rank()
    )
    common = {
        "solution_index": solution_index,
        "seed_X_coefficients_low_to_high": seed_x,
        "seed_Y_coefficients_low_to_high": seed_y,
        "tangent_rank": tangent_rank,
        "first_order_augmented_rank": augmented_rank,
        "lifts_mod_prime_squared": augmented_rank == tangent_rank,
    }
    if augmented_rank != tangent_rank:
        lift_records.append(
            {**common, "exact_local_obstruction": "no lift modulo p^2"}
        )
        continue
    if tangent_rank != len(seed):
        lift_records.append(
            {
                **common,
                "exact_local_obstruction": None,
                "unique_hensel_lift": False,
                "first_order_lift_dimension": len(seed) - tangent_rank,
            }
        )
        continue

    # The Jacobian has full column rank, so at each p-adic digit there is at
    # most one correction.  Test the *entire* overdetermined residual system,
    # not merely an invertible minor: failure of compatibility at level e is
    # an exact proof that the seed has no lift modulo p^(e+1).
    state = [ZZ(value) for value in seed]
    hensel_levels = []
    obstruction_exponent = None
    final_modulus = ZZ(prime)
    for exponent in range(1, arguments.hensel_depth):
        next_modulus = ZZ(prime) ** (exponent + 1)
        residue_ring = Integers(next_modulus)
        z_mod = vector(residue_ring, state)
        A_mod = [residue_ring(value.numerator()) / residue_ring(value.denominator()) for value in A_chart]
        B_mod = [residue_ring(value.numerator()) / residue_ring(value.denominator()) for value in B_chart]
        residual_mod, unused_jacobian, unused_y = residual_and_jacobian(
            z_mod, A_mod, B_mod
        )
        divisor = ZZ(prime) ** exponent
        lifted_residual = [ZZ(value.lift()) for value in residual_mod]
        if any(value % divisor for value in lifted_residual):
            raise ArithmeticError("Hensel state lost its residual invariant")
        target = vector(
            field, [field(-(value // divisor)) for value in lifted_residual]
        )
        compatible = (
            finite_matrix.augment(target.column()).rank() == finite_matrix.rank()
        )
        hensel_levels.append(
            {
                "input_exponent": exponent,
                "output_exponent": exponent + 1,
                "compatible": bool(compatible),
            }
        )
        if not compatible:
            obstruction_exponent = exponent + 1
            break
        correction = finite_matrix.solve_right(target)
        state = [
            (state[index] + divisor * ZZ(correction[index])) % next_modulus
            for index in range(len(state))
        ]
        final_modulus = next_modulus

    if obstruction_exponent is not None:
        lift_records.append(
            {
                **common,
                "exact_local_obstruction": (
                    f"no lift modulo {prime}^{obstruction_exponent}"
                ),
                "hensel_levels": hensel_levels,
            }
        )
        continue

    reconstruction = None
    try:
        reconstructed = [
            ZZ(value).rational_reconstruction(final_modulus) for value in state
        ]
        exact_residual, unused_jacobian, Y_chart_values = residual_and_jacobian(
            vector(QQ, reconstructed), list(A_chart), list(B_chart)
        )
        if not any(exact_residual):
            X_chart = ring(reconstructed[: x_degree + 1])
            Y_chart = ring(Y_chart_values)
            if chart_parameter is None:
                X_original = X_chart
                Y_original = Y_chart
            else:
                fraction_field = ring.fraction_field()
                inverse_parameter = fraction_field(1) / (t - chart_parameter)
                X_original = require_polynomial(
                    (t - chart_parameter) ** x_degree * X_chart(inverse_parameter),
                    ring,
                    "X",
                )
                Y_original = require_polynomial(
                    (t - chart_parameter) ** y_degree * Y_chart(inverse_parameter),
                    ring,
                    "Y",
                )
            if Y_original**2 != X_original**3 + A_twist * X_original + B_twist:
                raise ArithmeticError("reconstructed original-chart section fails")
            reconstruction = {
                "chart_X_coefficients_low_to_high": [str(value) for value in X_chart],
                "chart_Y_coefficients_low_to_high": [str(value) for value in Y_chart],
                "original_X_coefficients_low_to_high": [str(value) for value in X_original],
                "original_Y_coefficients_low_to_high": [str(value) for value in Y_original],
                "literal_curve_substitution": True,
            }
            exact_sections.append({"solution_index": solution_index, **reconstruction})
    except (ArithmeticError, ValueError, ZeroDivisionError):
        reconstruction = None
    lift_records.append(
        {
            **common,
            "exact_local_obstruction": None,
            "unique_hensel_lift": True,
            "hensel_levels": hensel_levels,
            "lifted_through_exponent": arguments.hensel_depth,
            "exact_rational_reconstruction": reconstruction,
        }
    )

payload = {
    "schema": "elkies-k3.r17-norm12-direct-singleton-po0-bruteforce-hensel.v1",
    "status": "PASS_EXACT_BRUTEFORCE_SEED_HENSEL_AUDIT",
    "candidate": export["candidate"],
    "prime": prime,
    "chart": chart,
    "hensel_depth": arguments.hensel_depth,
    "selected_solution_indices": selected_indices,
    "exact_local_obstruction_count": sum(
        item.get("exact_local_obstruction") is not None for item in lift_records
    ),
    "unique_hensel_lift_count": sum(
        item.get("unique_hensel_lift", False) for item in lift_records
    ),
    "exact_rational_section_count": len(exact_sections),
    "exact_sections": exact_sections,
    "lifts": lift_records,
    "proof_boundary": (
        "A first-order incompatibility is an exact local nonlifting theorem. "
        "A compatible isolated branch has a unique p-adic lift, but failure of "
        "bounded rational reconstruction is not a characteristic-zero exclusion."
    ),
    "inputs": {
        relative(path): digest(path)
        for path in (bruteforce_path, export_path, bisections_path, model_path)
    },
}
output_path = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17DIRECTPO0BRUTELIFT|label={label}|p={prime}"
    f"|selected={len(selected_indices)}"
    f"|obstructed={payload['exact_local_obstruction_count']}"
    f"|unique={payload['unique_hensel_lift_count']}"
    f"|exact={payload['exact_rational_section_count']}"
    f"|output={relative(output_path)}|status=PASS_AUDIT",
    flush=True,
)

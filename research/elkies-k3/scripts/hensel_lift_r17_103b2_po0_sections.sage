#!/usr/bin/env sage-python
"""Hensel-lift reduced p=29 P.O=0 branches for the 0x103b2 twist."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from sage.all import GF, Qp, QQ, ZZ, PolynomialRing, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BRUTEFORCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-norm12-orbit-103b2-p29-polynomial-section-bruteforce-v1.json"
)
DEFAULT_EXPORT = (
    ROOT
    / "artifacts/local/elkies-k3/twist-polynomial-sections/genus-one-norm12-orbit-103b2/p29/export.json"
)
DEFAULT_CONSTRUCTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
DEFAULT_MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-norm12-orbit-103b2-p29-hensel-lifts-v1.json"
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def convolution(left, right):
    result = [left[0].parent()(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return result


def coefficient(values, index):
    return values[index] if index < len(values) else values[0].parent()(0)


def residual_and_jacobian(z, A, B):
    X = list(z[:9])
    leading_y = z[9]
    zero = leading_y.parent()(0)
    X2 = convolution(X, X)
    rhs = convolution(X2, X)
    AX = convolution(A, X)
    rhs += [zero] * (25 - len(rhs))
    for degree, value in enumerate(AX):
        rhs[degree] += value
    for degree, value in enumerate(B):
        rhs[degree] += value

    kernel = [3 * coefficient(X2, degree) + coefficient(A, degree) for degree in range(17)]
    derivative_rhs = []
    for variable in range(9):
        derivative_rhs.append([zero] * variable + kernel + [zero] * (8 - variable))
    derivative_rhs.append([zero] * 25)

    Y = [zero] * 13
    derivative_y = [[zero] * 10 for unused in range(13)]
    Y[12] = leading_y
    derivative_y[12][9] = leading_y.parent()(1)
    denominator = 2 * leading_y
    for degree in range(23, 11, -1):
        index = degree - 12
        known = zero
        for left in range(index + 1, 13):
            right = degree - left
            if index < right <= 12:
                known += Y[left] * Y[right]
        Y[index] = (rhs[degree] - known) / denominator
        for variable in range(10):
            derivative_known = zero
            for left in range(index + 1, 13):
                right = degree - left
                if index < right <= 12:
                    derivative_known += (
                        derivative_y[left][variable] * Y[right]
                        + Y[left] * derivative_y[right][variable]
                    )
            denominator_derivative = 2 if variable == 9 else 0
            derivative_y[index][variable] = (
                derivative_rhs[variable][degree]
                - derivative_known
                - Y[index] * denominator_derivative
            ) / denominator

    Y2 = convolution(Y, Y)
    degrees = [24, *range(12)]
    residual = [Y2[degree] - rhs[degree] for degree in degrees]
    jacobian = []
    for degree in degrees:
        row = []
        for variable in range(10):
            derivative = zero
            for left in range(13):
                right = degree - left
                if 0 <= right <= 12:
                    derivative += (
                        derivative_y[left][variable] * Y[right]
                        + Y[left] * derivative_y[right][variable]
                    )
            row.append(derivative - derivative_rhs[variable][degree])
        jacobian.append(row)
    return residual, jacobian, Y


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--bruteforce", type=Path, default=DEFAULT_BRUTEFORCE)
parser.add_argument("--export", type=Path, default=DEFAULT_EXPORT)
parser.add_argument("--constructions", type=Path, default=DEFAULT_CONSTRUCTIONS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--precision", type=int, default=220)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

bruteforce = json.loads(args.bruteforce.read_text())
export = json.loads(args.export.read_text())
prime = int(bruteforce["prime"])
if prime != 29 or int(export["prime"]) != prime:
    raise ValueError("this certificate expects the p=29 shell")

model = json.loads(args.model.read_text())
constructions = json.loads(args.constructions.read_text())
record = next(
    item
    for item in constructions["construction"]["records"]
    if item["label"] == "norm12-orbit-103b2"
)
ring = PolynomialRing(QQ, "t")
t = ring.gen()
A0 = ring([QQ(value) for value in model["A_coefficients_low_to_high"]])
B0 = ring([QQ(value) for value in model["B_coefficients_low_to_high"]])
q = ring([QQ(value) for value in record["branch_polynomial_q_coefficients_low_to_high"]])
denominator = lcm(*(value.denominator() for value in q.list()))
q_integral = ring(denominator**2 * q)
A_twist = A0 * q_integral**2
B_twist = B0 * q_integral**3
chart = export["infinity_fibre"]["chart"]
prefix = "original_t="
if not chart.startswith(prefix) or " via " not in chart:
    raise ValueError("expected a finite chart in the p=29 export")
chart_parameter = QQ(int(chart[len(prefix) :].split(" via ", 1)[0]))
A_chart = ring(
    sum(
        A_twist[index] * (chart_parameter * t + 1) ** index * t ** (16 - index)
        for index in range(A_twist.degree() + 1)
    )
)
B_chart = ring(
    sum(
        B_twist[index] * (chart_parameter * t + 1) ** index * t ** (24 - index)
        for index in range(B_twist.degree() + 1)
    )
)
if [int(value) % prime for value in A_chart.list()] != export[
    "twist_A_coefficients_low_to_high"
]:
    raise ArithmeticError("exact A chart does not reduce to the export")
if [int(value) % prime for value in B_chart.list()] != export[
    "twist_B_coefficients_low_to_high"
]:
    raise ArithmeticError("exact B chart does not reduce to the export")

padics = Qp(prime, prec=args.precision, type="capped-rel")
A_padic = [padics(value) for value in A_chart.list()]
B_padic = [padics(value) for value in B_chart.list()]
field = GF(prime)
lift_records = []
exact_sections = []
first_order_records = []
for solution_index, solution in enumerate(bruteforce["solutions"]):
    seed = [*solution["X_coefficients_low_to_high"], solution["Y_coefficients_low_to_high"][12]]
    rational_residual, rational_jacobian, unused_y = residual_and_jacobian(
        vector(QQ, seed), list(A_chart), list(B_chart)
    )
    if any(value.valuation(prime) < 1 for value in rational_residual if value):
        raise ArithmeticError("mod-p shell point failed exact residual replay")
    finite_jacobian = matrix(
        field, [[field(value) for value in row] for row in rational_jacobian]
    )
    target = vector(
        field,
        [field(-value / prime) for value in rational_residual],
    )
    rank = finite_jacobian.rank()
    augmented_rank = finite_jacobian.augment(target.column()).rank()
    first_order_records.append(
        {
            "solution_index": solution_index,
            "tangent_rank": int(rank),
            "augmented_rank": int(augmented_rank),
            "lifts_to_mod_p_squared": rank == augmented_rank,
            "first_order_lift_dimension": 10 - int(rank) if rank == augmented_rank else None,
        }
    )
for solution_index, solution in enumerate(bruteforce["solutions"]):
    if int(solution["full_shell_tangent_rank"]) != 10:
        continue
    seed = [*solution["X_coefficients_low_to_high"], solution["Y_coefficients_low_to_high"][12]]
    finite_residual, finite_jacobian, unused_y = residual_and_jacobian(
        vector(field, seed),
        [field(value) for value in export["twist_A_coefficients_low_to_high"]],
        [field(value) for value in export["twist_B_coefficients_low_to_high"]],
    )
    finite_matrix = matrix(field, finite_jacobian)
    pivot_rows = list(finite_matrix.transpose().pivots())
    if len(pivot_rows) != 10:
        raise ArithmeticError("stored full-rank branch failed replay")

    z = vector(padics, seed)
    valuations = []
    for unused_iteration in range(12):
        residual, jacobian, unused_y = residual_and_jacobian(z, A_padic, B_padic)
        finite_values = [value.valuation() for value in residual if value]
        valuations.append(min(finite_values) if finite_values else args.precision)
        square_matrix = matrix(
            padics, [[jacobian[row][column] for column in range(10)] for row in pivot_rows]
        )
        correction = square_matrix.solve_right(
            -vector(padics, [residual[row] for row in pivot_rows])
        )
        z += correction
        if valuations[-1] >= args.precision // 2:
            break

    residual, unused_jacobian, Y_padic = residual_and_jacobian(z, A_padic, B_padic)
    final_valuation = min(
        [value.valuation() for value in residual if value] or [args.precision]
    )
    reconstruction = None
    try:
        modulus = ZZ(prime) ** (args.precision // 2)
        reconstructed = [ZZ(value.lift()).rational_reconstruction(modulus) for value in z]
        exact_residual, unused_jacobian, exact_Y = residual_and_jacobian(
            vector(QQ, reconstructed), list(A_chart), list(B_chart)
        )
        if not any(exact_residual):
            reconstruction = {
                "X_coefficients_low_to_high": [str(value) for value in reconstructed[:9]],
                "Y_coefficients_low_to_high": [str(value) for value in exact_Y],
                "leading_y": str(reconstructed[9]),
            }
            exact_sections.append({"solution_index": solution_index, **reconstruction})
    except (ArithmeticError, ValueError, ZeroDivisionError):
        reconstruction = None
    lift_records.append(
        {
            "solution_index": solution_index,
            "seed": seed,
            "pivot_rows": pivot_rows,
            "newton_residual_valuations": [int(value) for value in valuations],
            "final_residual_valuation_at_least": int(final_valuation),
            "exact_rational_reconstruction": reconstruction,
        }
    )

payload = {
    "schema": "elkies-k3.r17-103b2-po0-hensel-lifts.v1",
    "status": "PASS_EXACT_HENSEL_LIFT_AUDIT",
    "prime": prime,
    "precision": args.precision,
    "full_rank_branch_count": len(lift_records),
    "first_order_mod_p_squared_audit": {
        "branch_count": len(first_order_records),
        "lifting_branch_count": sum(
            record["lifts_to_mod_p_squared"] for record in first_order_records
        ),
        "records": first_order_records,
    },
    "exact_rational_section_count": len(exact_sections),
    "exact_sections": exact_sections,
    "lifts": lift_records,
    "proof_boundary": (
        "A failed rational reconstruction means only that the displayed p-adic branch did "
        "not reconstruct within this bound. This does not exclude algebraic or larger rational "
        "coordinates, and the P.O=0 shell does not bound the full Mordell-Weil rank."
    ),
    "inputs": {
        relative(path): digest(path)
        for path in (args.bruteforce, args.export, args.constructions, args.model)
    },
}
output = args.output if args.output.is_absolute() else ROOT / args.output
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17103B2HENSEL|branches={len(lift_records)}|exact={len(exact_sections)}"
    f"|output={relative(output)}|status=PASS_AUDIT"
)

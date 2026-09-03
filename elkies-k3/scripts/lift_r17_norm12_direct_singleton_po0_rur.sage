#!/usr/bin/env sage-python
"""Decode and Hensel-lift a direct norm-12 singleton-twist P.O=0 RUR.

The modular exporter eliminates the nonleading Y coefficients recursively.
This script decodes every rational point in a selected zero-dimensional msolve
RUR, replays all residual equations, and attempts a unique p-adic lift along
six independent rows.  The three unused equations are retained throughout as
compatibility checks.  Coefficientwise rational reconstruction is accepted
only after exact substitution over QQ.
"""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
from itertools import product
import json
from math import lcm
from pathlib import Path
import sys

from sage.all import GF, Integers, Qp, QQ, ZZ, PolynomialRing, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from screen_elkies_2026_quadratic_twist_ranks import (  # noqa: E402
    square_equivalent_integer_polynomial,
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
    """Return the leading and nine low residuals in the full eight variables."""

    X = list(z[:7])
    leading_y = z[7]
    zero = leading_y.parent()(0)
    x_degree = len(X) - 1
    y_degree = 3 * x_degree // 2
    X2 = convolution(X, X)
    rhs = convolution(X2, X)
    AX = convolution(A, X)
    rhs += [zero] * (2 * y_degree + 1 - len(rhs))
    for degree, value in enumerate(AX):
        rhs[degree] += value
    for degree, value in enumerate(B):
        rhs[degree] += value

    kernel = [
        3 * coefficient(X2, degree) + coefficient(A, degree)
        for degree in range(2 * x_degree + 1)
    ]
    derivative_rhs = []
    for variable in range(x_degree + 1):
        row = [zero] * (2 * y_degree + 1)
        for degree, value in enumerate(kernel):
            row[degree + variable] += value
        derivative_rhs.append(row)
    derivative_rhs.append([zero] * (2 * y_degree + 1))

    Y = [zero] * (y_degree + 1)
    variable_count = x_degree + 2
    derivative_y = [[zero] * variable_count for unused in range(y_degree + 1)]
    Y[y_degree] = leading_y
    derivative_y[y_degree][-1] = leading_y.parent()(1)
    denominator = 2 * leading_y
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
            denominator_derivative = 2 if variable == variable_count - 1 else 0
            derivative_y[index][variable] = (
                derivative_rhs[variable][degree]
                - derivative_known
                - Y[index] * denominator_derivative
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


def residual_only(z, A, B):
    X = list(z[:7])
    leading_y = z[7]
    zero = leading_y.parent()(0)
    X2 = convolution(X, X)
    rhs = convolution(X2, X)
    AX = convolution(A, X)
    rhs += [zero] * (19 - len(rhs))
    for degree, value in enumerate(AX):
        rhs[degree] += value
    for degree, value in enumerate(B):
        rhs[degree] += value
    Y = [zero] * 10
    Y[9] = leading_y
    for degree in range(17, 8, -1):
        index = degree - 9
        known = sum(
            Y[left] * Y[degree - left]
            for left in range(index + 1, 10)
            if index < degree - left <= 9
        )
        Y[index] = (rhs[degree] - known) / (2 * leading_y)
    Y2 = convolution(Y, Y)
    return [Y2[degree] - rhs[degree] for degree in [18, *range(9)]]


def reduce_rational_mod(value, residue_ring):
    value = QQ(value)
    return residue_ring(value.numerator()) / residue_ring(value.denominator())


def hensel_tree(seed, finite_matrix, A, B, prime, depth, max_states):
    """Enumerate all digit lifts through p^depth, unless the cap is reached."""

    finite_field = GF(prime)
    kernel_basis = list(finite_matrix.right_kernel().basis())
    kernel_dimension = len(kernel_basis)
    states = [tuple(int(value) % prime for value in seed)]
    levels = []
    truncated = False
    for exponent in range(1, depth):
        next_modulus = ZZ(prime) ** (exponent + 1)
        residue_ring = Integers(next_modulus)
        A_mod = [reduce_rational_mod(value, residue_ring) for value in A]
        B_mod = [reduce_rational_mod(value, residue_ring) for value in B]
        compatible_states = 0
        next_states = []
        for state in states:
            residual = residual_only(vector(residue_ring, state), A_mod, B_mod)
            divisor = ZZ(prime) ** exponent
            lifted_residuals = [ZZ(value.lift()) for value in residual]
            if any(value % divisor for value in lifted_residuals):
                raise ArithmeticError("Hensel tree state lost its residual invariant")
            target = vector(
                finite_field,
                [finite_field(-(value // divisor)) for value in lifted_residuals],
            )
            if finite_matrix.augment(target.column()).rank() != finite_matrix.rank():
                continue
            compatible_states += 1
            particular = finite_matrix.solve_right(target)
            for parameters in product(range(prime), repeat=kernel_dimension):
                delta = particular + sum(
                    (
                        finite_field(parameter) * basis_vector
                        for parameter, basis_vector in zip(parameters, kernel_basis)
                    ),
                    vector(finite_field, [0] * 8),
                )
                next_states.append(
                    tuple(
                        int((ZZ(state[index]) + divisor * ZZ(delta[index])) % next_modulus)
                        for index in range(8)
                    )
                )
                if len(next_states) >= max_states:
                    truncated = True
                    break
            if truncated:
                break
        levels.append(
            {
                "input_solution_count_mod_p_power": len(states),
                "input_exponent": exponent,
                "compatible_parent_count": compatible_states,
                "output_solution_count_mod_p_power": len(next_states),
                "output_exponent": exponent + 1,
                "truncated_at_state_cap": truncated,
            }
        )
        states = next_states
        if truncated or not states:
            break
    return {
        "target_exponent": depth,
        "kernel_dimension": kernel_dimension,
        "state_cap": max_states,
        "truncated": truncated,
        "levels": levels,
        "surviving_state_count": len(states),
        "surviving_state_sample": [list(state) for state in states[:20]],
    }


def decode_rational_rur_points(solution_path, prime, variable_names):
    text = solution_path.read_text().strip()
    solution = ast.literal_eval(text[:-1] if text.endswith(":") else text)
    if solution[0] != 0:
        raise ArithmeticError("msolve did not return a zero-dimensional RUR")
    payload = solution[1]
    if int(payload[0]) != prime:
        raise ArithmeticError("RUR characteristic disagrees with the export")
    added_separator = payload[3][-1] == "A"
    if (
        payload[3][:-1] != variable_names
        if added_separator
        else payload[3] != variable_names
    ):
        raise ArithmeticError("unexpected RUR variable order")
    parametrization = payload[5]
    if parametrization[0] != 1:
        raise ArithmeticError("unexpected number of RUR blocks")
    elimination_data, denominator_data, coordinate_data = parametrization[1]
    expected_coordinate_count = len(variable_names) if added_separator else len(variable_names) - 1
    if denominator_data != [0, [1]] or len(coordinate_data) != expected_coordinate_count:
        raise ArithmeticError("unsupported RUR denominator or coordinate count")
    if not added_separator and payload[4] != [0, 0, 0, 0, 0, 1]:
        raise ArithmeticError("unsupported native separating variable")

    field = GF(prime)
    ring = PolynomialRing(field, "R")
    elimination = ring(elimination_data[1]).squarefree_part()
    coordinate_polynomials = [ring(block[0][1]) for block in coordinate_data]
    points = []
    for factor, exponent in elimination.factor():
        if exponent != 1:
            raise ArithmeticError("square-free RUR factor has multiplicity")
        if factor.degree() != 1:
            continue
        for root in factor.roots(multiplicities=False):
            values = [
                -sum(
                    field(value) * root**index
                    for index, value in enumerate(poly.list())
                )
                for poly in coordinate_polynomials
            ]
            if not added_separator:
                values.append(root)
            points.append(
                {
                    "elimination_factor": str(factor),
                    "root": str(root),
                    "variable_values": [int(value) for value in values],
                }
            )
    return payload, elimination, points


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--export", type=Path, required=True)
parser.add_argument("--solution", type=Path, required=True)
parser.add_argument("--bisections", type=Path, required=True)
parser.add_argument("--model", type=Path, required=True)
parser.add_argument("--precision", type=int, default=240)
parser.add_argument("--tree-depth", type=int, default=4)
parser.add_argument("--max-tree-states", type=int, default=200000)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
if args.tree_depth < 2 or args.max_tree_states <= 0:
    parser.error("--tree-depth must be at least two and --max-tree-states positive")

export = json.loads(args.export.read_text())
if export.get("schema") != "elkies-k3.elkies-2026-twist-polynomial-section-msolve-export.v1":
    raise ValueError("unexpected modular export schema")
if export["candidate"].get("kind") != "direct_singleton":
    raise ValueError("this lift expects a direct singleton twist")
prime = int(export["prime"])
field = GF(prime)

system_record = next(
    item
    for item in export["systems"]
    if Path(item["path"]).stem == args.solution.stem
)
block_index = int(system_record["block_index"])
leading_x_integer, leading_y_integer = system_record["leading_x_y"]
leading_x_finite = field(leading_x_integer)
leading_y_finite = field(leading_y_integer)
variable_names = [f"x{index}" for index in range(5, -1, -1)]
payload, elimination, decoded_points = decode_rational_rur_points(
    args.solution, prime, variable_names
)

bisections = json.loads(args.bisections.read_text())
label = export["candidate"]["key"]
record = next(item for item in bisections["bisections"] if item["label"] == label)
q_coefficients = square_equivalent_integer_polynomial(
    record["branch"]["numerator_coefficients"]
)
model = json.loads(args.model.read_text())
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

chart = export["infinity_fibre"]["chart"]
if chart == "original_infinity":
    A_chart = A_twist
    B_chart = B_twist
elif chart.startswith("original_t=") and " via " in chart:
    chart_parameter = QQ(int(chart[len("original_t=") :].split(" via ", 1)[0]))
    A_chart = ring(
        sum(
            A_twist[index]
            * (chart_parameter * t + 1) ** index
            * t ** (12 - index)
            for index in range(A_twist.degree() + 1)
        )
    )
    B_chart = ring(
        sum(
            B_twist[index]
            * (chart_parameter * t + 1) ** index
            * t ** (18 - index)
            for index in range(B_twist.degree() + 1)
        )
    )
else:
    raise ValueError("unrecognized export chart")
if [int(field(QQ(value))) for value in A_chart.list()] != export[
    "twist_A_coefficients_low_to_high"
]:
    raise ArithmeticError("exact A chart does not reduce to the export")
if [int(field(QQ(value))) for value in B_chart.list()] != export[
    "twist_B_coefficients_low_to_high"
]:
    raise ArithmeticError("exact B chart does not reduce to the export")

padics = Qp(prime, prec=args.precision, type="capped-rel")
A_padic = [padics(value) for value in A_chart.list()]
B_padic = [padics(value) for value in B_chart.list()]

lift_records = []
exact_sections = []
for point_index, point in enumerate(decoded_points):
    # RUR variables are x5,...,x0; the recurrence uses low-to-high order.
    seed_x = [*reversed(point["variable_values"]), leading_x_integer]
    seed = [*seed_x, leading_y_integer]
    residual_finite, jacobian_finite, Y_finite = residual_and_jacobian(
        vector(field, seed),
        [field(value) for value in export["twist_A_coefficients_low_to_high"]],
        [field(value) for value in export["twist_B_coefficients_low_to_high"]],
    )
    if any(residual_finite):
        raise ArithmeticError("decoded RUR point fails the exported equations")
    finite_matrix = matrix(field, jacobian_finite)
    tangent_rank = int(finite_matrix.rank())
    residual_rational, jacobian_rational, unused_Y_rational = residual_and_jacobian(
        vector(QQ, seed),
        list(A_chart),
        list(B_chart),
    )
    if any(value and value.valuation(prime) < 1 for value in residual_rational):
        raise ArithmeticError("decoded point has a nonintegral first-order residual")
    first_order_target = vector(
        field, [field(-value / prime) for value in residual_rational]
    )
    augmented_rank = int(
        finite_matrix.augment(first_order_target.column()).rank()
    )
    lifts_mod_prime_squared = augmented_rank == tangent_rank
    tree = hensel_tree(
        seed,
        finite_matrix,
        list(A_chart),
        list(B_chart),
        prime,
        args.tree_depth,
        args.max_tree_states,
    )
    pivot_rows = list(finite_matrix.transpose().pivots())
    common_record = {
        **point,
        "point_index": point_index,
        "seed_X_coefficients_low_to_high": seed_x,
        "seed_leading_Y": leading_y_integer,
        "tangent_rank": tangent_rank,
        "first_order_augmented_rank": augmented_rank,
        "lifts_mod_prime_squared": lifts_mod_prime_squared,
        "hensel_tree": tree,
    }
    if not lifts_mod_prime_squared:
        lift_records.append(
            {
                **common_record,
                "unique_hensel_lift": False,
                "obstruction": "augmented_rank_exceeds_tangent_rank",
            }
        )
        continue
    if tangent_rank != 8 or len(pivot_rows) != 8:
        lift_records.append(
            {
                **common_record,
                "unique_hensel_lift": False,
                "first_order_lift_dimension": 8 - tangent_rank,
                "obstruction": None,
            }
        )
        continue

    z = vector(padics, seed)
    valuations = []
    compatibility_valuations = []
    for unused_iteration in range(12):
        residual, jacobian, unused_Y = residual_and_jacobian(
            z, A_padic, B_padic
        )
        all_valuations = [value.valuation() for value in residual if value]
        valuations.append(min(all_valuations) if all_valuations else args.precision)
        unused_rows = [row for row in range(10) if row not in pivot_rows]
        unused_valuations = [
            residual[row].valuation() for row in unused_rows if residual[row]
        ]
        compatibility_valuations.append(
            min(unused_valuations) if unused_valuations else args.precision
        )
        square_matrix = matrix(
            padics,
            [[jacobian[row][column] for column in range(8)] for row in pivot_rows],
        )
        correction = square_matrix.solve_right(
            -vector(padics, [residual[row] for row in pivot_rows])
        )
        z += correction
        if valuations[-1] >= args.precision // 2:
            break

    residual, unused_jacobian, Y_padic = residual_and_jacobian(
        z, A_padic, B_padic
    )
    final_valuation = min(
        [value.valuation() for value in residual if value] or [args.precision]
    )
    reconstruction = None
    try:
        modulus = ZZ(prime) ** (args.precision // 2)
        reconstructed = [
            ZZ(value.lift()).rational_reconstruction(modulus) for value in z
        ]
        exact_residual, unused_jacobian, exact_Y = residual_and_jacobian(
            vector(QQ, reconstructed),
            list(A_chart),
            list(B_chart),
        )
        if not any(exact_residual):
            reconstruction = {
                "X_coefficients_low_to_high": [
                    str(value) for value in reconstructed[:7]
                ],
                "Y_coefficients_low_to_high": [str(value) for value in exact_Y],
                "literal_curve_substitution": True,
            }
            exact_sections.append({"point_index": point_index, **reconstruction})
    except (ArithmeticError, ValueError, ZeroDivisionError):
        reconstruction = None
    lift_records.append(
        {
            **common_record,
            "unique_hensel_lift": True,
            "pivot_rows": pivot_rows,
            "newton_residual_valuations": [int(value) for value in valuations],
            "unused_equation_valuations": [
                int(value) for value in compatibility_valuations
            ],
            "final_all_equation_valuation_at_least": int(final_valuation),
            "exact_rational_reconstruction": reconstruction,
        }
    )

payload_out = {
    "schema": "elkies-k3.r17-norm12-direct-singleton-po0-rur-hensel.v1",
    "status": "PASS_EXACT_HENSEL_LIFT_AUDIT",
    "candidate": export["candidate"],
    "prime": prime,
    "block_index": block_index,
    "leading_x_y": [leading_x_integer, leading_y_integer],
    "precision": args.precision,
    "rur": {
        "quotient_dimension_with_multiplicity": int(payload[2]),
        "squarefree_elimination_polynomial": str(elimination),
        "rational_point_count": len(decoded_points),
    },
    "lift_count": sum(record.get("unique_hensel_lift", False) for record in lift_records),
    "exact_rational_section_count": len(exact_sections),
    "exact_sections": exact_sections,
    "lifts": lift_records,
    "proof_boundary": (
        "This exactly audits the rational support of the selected modular RUR. "
        "A failed compatibility lift rules out that nonsingular mod-p point as a "
        "characteristic-zero section in this integral model. Failed rational "
        "reconstruction alone is only a bounded negative result."
    ),
    "inputs": {
        relative(path): digest(path)
        for path in (args.export, args.solution, args.bisections, args.model)
    },
}
output = args.output if args.output.is_absolute() else ROOT / args.output
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload_out, indent=2, sort_keys=True) + "\n")
print(
    f"R17DIRECTPO0HENSEL|label={label}|p={prime}|block={block_index}"
    f"|points={len(decoded_points)}|unique={payload_out['lift_count']}"
    f"|exact={len(exact_sections)}|output={relative(output)}|status=PASS_AUDIT",
    flush=True,
)

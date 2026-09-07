#!/usr/bin/env sage-python
"""Hensel-lift one rational P.O=0 RUR block for a 074d9 twist."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, Qp, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
COVERS = ROOT / "artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"


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

    kernel = [
        3 * coefficient(X2, degree) + coefficient(A, degree)
        for degree in range(13)
    ]
    derivative_rhs = []
    for variable in range(7):
        row = [zero] * 19
        for degree, value in enumerate(kernel):
            row[degree + variable] += value
        derivative_rhs.append(row)
    derivative_rhs.append([zero] * 19)

    Y = [zero] * 10
    derivative_y = [[zero] * 8 for unused in range(10)]
    Y[9] = leading_y
    derivative_y[9][7] = leading_y.parent()(1)
    denominator = 2 * leading_y
    for degree in range(17, 8, -1):
        index = degree - 9
        known = zero
        for left in range(index + 1, 10):
            right = degree - left
            if index < right <= 9:
                known += Y[left] * Y[right]
        Y[index] = (rhs[degree] - known) / denominator
        for variable in range(8):
            derivative_known = zero
            for left in range(index + 1, 10):
                right = degree - left
                if index < right <= 9:
                    derivative_known += (
                        derivative_y[left][variable] * Y[right]
                        + Y[left] * derivative_y[right][variable]
                    )
            denominator_derivative = 2 if variable == 7 else 0
            derivative_y[index][variable] = (
                derivative_rhs[variable][degree]
                - derivative_known
                - Y[index] * denominator_derivative
            ) / denominator

    Y2 = convolution(Y, Y)
    residual_degrees = [18, *range(9)]
    residual = [Y2[degree] - rhs[degree] for degree in residual_degrees]
    jacobian = []
    for degree in residual_degrees:
        row = []
        for variable in range(8):
            derivative = zero
            for left in range(10):
                right = degree - left
                if 0 <= right <= 9:
                    derivative += (
                        derivative_y[left][variable] * Y[right]
                        + Y[left] * derivative_y[right][variable]
                    )
            row.append(derivative - derivative_rhs[variable][degree])
        jacobian.append(row)
    return residual, jacobian, Y


def decode_rational_points(solution_path, prime):
    text = solution_path.read_text().strip()
    solution = ast.literal_eval(text[:-1] if text.endswith(":") else text)
    if solution[0] != 0:
        raise ArithmeticError("msolve output is not a zero-dimensional RUR")
    payload = solution[1]
    variable_names = [f"x{index}" for index in range(5, -1, -1)]
    added_separator = payload[3][-1] == "A"
    reported_names = payload[3][:-1] if added_separator else payload[3]
    if int(payload[0]) != prime or reported_names != variable_names:
        raise ArithmeticError("RUR field or variable order mismatch")
    if not added_separator and payload[4] != [0, 0, 0, 0, 0, 1]:
        raise ArithmeticError("unsupported separating variable")
    parametrization = payload[5]
    elimination_data, denominator_data, coordinate_data = parametrization[1]
    expected_coordinates = 6 if added_separator else 5
    if (
        parametrization[0] != 1
        or denominator_data != [0, [1]]
        or len(coordinate_data) != expected_coordinates
    ):
        raise ArithmeticError("unsupported RUR parametrization")
    field = GF(prime)
    ring = PolynomialRing(field, "R")
    elimination = ring(elimination_data[1]).squarefree_part()
    coordinates = [ring(block[0][1]) for block in coordinate_data]
    points = []
    for factor, multiplicity in elimination.factor():
        if multiplicity != 1:
            raise ArithmeticError("RUR eliminant is not squarefree")
        if factor.degree() != 1:
            continue
        root = factor.roots(multiplicities=False)[0]
        values = [
            -sum(value * root**index for index, value in enumerate(poly.list()))
            for poly in coordinates
        ]
        if not added_separator:
            values.append(root)
        points.append([int(value) for value in values])
    return int(payload[2]), str(elimination), points


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--export", type=Path, required=True)
parser.add_argument("--solution", type=Path, required=True)
parser.add_argument("--precision", type=int, default=1600)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

export = json.loads(args.export.read_text())
if (
    export.get("schema") != "elkies-k3.r17-074d9-twist-section-ladder-msolve-export.v1"
    or export["intersection_P_dot_O"] != 0
):
    raise ValueError("this lifter requires a 074d9 P.O=0 export")
prime = int(export["prime"])
field = GF(prime)
system = next(
    row for row in export["systems"]
    if Path(row["path"]).stem == args.solution.stem
)
if digest(ROOT / system["path"]) != system["sha256"]:
    raise ArithmeticError("stale msolve input")
multiplicity, elimination, decoded = decode_rational_points(args.solution, prime)

covers = json.loads(COVERS.read_text())
source = next(
    row for fibre in covers["fibres"] for row in fibre["records"]
    if row["label"] == export["label"]
)
model = json.loads(MODEL.read_text())["representative"]
ring = PolynomialRing(QQ, "u")
u = ring.gen()
q = ring([QQ(value) for value in source["branch_quadratic_coefficients_low_to_high"]])
A0 = ring([QQ(value) for value in model["A_coefficients_low_to_high"]])
B0 = ring([QQ(value) for value in model["B_coefficients_low_to_high"]])
A = q**2 * A0
B = q**3 * B0
chart = export["charts"][system["chart_index"]]["original_fibre"]
if chart == "original_infinity":
    chart_parameter = None
    A_chart, B_chart = A, B
elif chart.startswith("t="):
    chart_parameter = QQ(int(chart[2:]))
    A_chart = ring(sum(A[index] * (chart_parameter * u + 1) ** index * u ** (12-index) for index in range(13)))
    B_chart = ring(sum(B[index] * (chart_parameter * u + 1) ** index * u ** (18-index) for index in range(19)))
else:
    raise ValueError("unrecognized chart")

padics = Qp(prime, prec=args.precision, type="capped-rel")
A_padic = [padics(value) for value in A_chart.list()]
B_padic = [padics(value) for value in B_chart.list()]
leading_x, leading_y = system["leading_x_y"]
lift_rows = []
exact_sections = []
for point_index, variable_values in enumerate(decoded):
    seed = [*reversed(variable_values), leading_x, leading_y]
    residual_mod_p, jacobian_mod_p, unused_Y = residual_and_jacobian(
        vector(field, seed), [field(value) for value in A_chart.list()],
        [field(value) for value in B_chart.list()]
    )
    if any(residual_mod_p):
        raise ArithmeticError("decoded point fails the section equations")
    finite_matrix = matrix(field, jacobian_mod_p)
    tangent_rank = int(finite_matrix.rank())
    row = {
        "point_index": point_index,
        "seed_X_coefficients_low_to_high": seed[:7],
        "seed_leading_Y": leading_y,
        "tangent_rank": tangent_rank,
    }
    if tangent_rank != 8:
        row["outcome"] = "NONUNIQUE_OR_SINGULAR_HENSEL_BRANCH"
        lift_rows.append(row)
        continue
    pivot_rows = list(finite_matrix.transpose().pivots())
    z = vector(padics, seed)
    valuations = []
    for unused_iteration in range(16):
        residual, jacobian, unused_Y = residual_and_jacobian(z, A_padic, B_padic)
        valuations.append(min([value.valuation() for value in residual if value] or [args.precision]))
        square = matrix(padics, [[jacobian[r][c] for c in range(8)] for r in pivot_rows])
        z += square.solve_right(-vector(padics, [residual[r] for r in pivot_rows]))
        if valuations[-1] >= args.precision // 2:
            break
    residual, unused_jacobian, Y_padic = residual_and_jacobian(z, A_padic, B_padic)
    final_valuation = min([value.valuation() for value in residual if value] or [args.precision])
    if final_valuation < args.precision // 2:
        row.update(
            {
                "outcome": "OBSTRUCTED_MOD_PRIME_POWER",
                "pivot_rows": pivot_rows,
                "newton_residual_valuations": [int(value) for value in valuations],
                "first_incompatible_modulus_exponent": int(final_valuation) + 1,
                "exact_rational_reconstruction": None,
            }
        )
        lift_rows.append(row)
        continue
    reconstruction = None
    try:
        modulus = ZZ(prime) ** (args.precision // 2)
        reconstructed = [ZZ(value.lift()).rational_reconstruction(modulus) for value in z]
        exact_residual, unused_jacobian, exact_Y = residual_and_jacobian(
            vector(QQ, reconstructed), list(A_chart), list(B_chart)
        )
        if not any(exact_residual):
            X_chart = ring(reconstructed[:7])
            Y_chart = ring(exact_Y)
            if chart_parameter is None:
                X_original, Y_original = X_chart, Y_chart
            else:
                fraction_field = ring.fraction_field()
                v = fraction_field.gen()
                X_original = ring((v-chart_parameter)**6 * X_chart(1/(v-chart_parameter)))
                Y_original = ring((v-chart_parameter)**9 * Y_chart(1/(v-chart_parameter)))
            if Y_original**2 != X_original**3 + A*X_original + B:
                raise ArithmeticError("inverse chart section fails exact substitution")
            reconstruction = {
                "X_coefficients_low_to_high": [str(value) for value in X_original.list()],
                "Y_coefficients_low_to_high": [str(value) for value in Y_original.list()],
                "literal_curve_substitution": True,
            }
            exact_sections.append(reconstruction)
    except (ArithmeticError, ValueError, ZeroDivisionError):
        reconstruction = None
    row.update(
        {
            "outcome": "UNIQUE_HENSEL_LIFT",
            "pivot_rows": pivot_rows,
            "newton_residual_valuations": [int(value) for value in valuations],
            "final_all_equation_valuation_at_least": int(final_valuation),
            "exact_rational_reconstruction": reconstruction,
        }
    )
    lift_rows.append(row)

payload = {
    "schema": "elkies-k3.r17-074d9-twist-po0-rur-hensel.v1",
    "status": "PASS_EXACT_HENSEL_LIFT_AUDIT",
    "label": export["label"],
    "prime": prime,
    "chart": chart,
    "block_index": system["block_index"],
    "leading_x_y": system["leading_x_y"],
    "precision": args.precision,
    "rur_multiplicity": multiplicity,
    "squarefree_eliminant": elimination,
    "rational_point_count": len(decoded),
    "lifts": lift_rows,
    "exact_rational_sections": exact_sections,
    "proof_boundary": (
        "Every rational point of the selected zero-dimensional RUR is decoded. "
        "Nonsingular branches are Hensel-lifted; only literal QQ substitution is "
        "accepted as an exact section. Failed rational reconstruction is bounded "
        "negative evidence only."
    ),
    "inputs": {
        relative(path): digest(path)
        for path in (args.export, args.solution, COVERS, MODEL)
    },
}
output = args.output if args.output.is_absolute() else ROOT / args.output
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17074D9PO0HENSEL|label={export['label']}|block={system['block_index']}"
    f"|points={len(decoded)}|exact={len(exact_sections)}|output={relative(output)}",
    flush=True,
)

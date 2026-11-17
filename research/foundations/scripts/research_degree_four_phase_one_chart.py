#!/usr/bin/env python3
"""Explore the eight-dimensional two-direction phase-one moment chart.

The phase +1 and phase -1 tau-even direction spaces both have dimension
two.  On the open set where the first positive coefficient is nonzero,
the residual diagonal torus gives the birational section

    diag(a0,...,a4) + B0 + u*B1 + v*C0 + w*C1.

Thus the quotient chart has coordinates a0,...,a4,u,v,w.  The apolar
involution reverses the a_i and fixes u,v,w.

This script constructs the restricted moments exactly and asks msolve for
the fiber of a selected moment prefix over a prime field.  It is an
exploratory modular screen, not a characteristic-zero degree certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import tempfile
from fractions import Fraction
from functools import reduce
from math import factorial, gcd
from operator import mul
from pathlib import Path

import sympy as sp

import verify_two_pair_counterexample_missing_invariant as invariant_base


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_phase_one_chart_modular.json"
)
VARIABLES = ("a0", "a1", "a2", "a3", "a4", "u", "v", "w")
DEFAULT_POINT = (2, 3, 5, 7, 11, 2, 3, 5)
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, int]
Direction = dict[tuple[int, int], int]
EVEN_CUBIC_TRIPLES = tuple(
    triple
    for triple in invariant_base.CUBIC_TRIPLES
    if triple != (2, 3, 4)
)
EVEN_QUARTIC_WORDS = (
    (1, 1, 1, 3),
    (1, 1, 2, 2),
    (1, 1, 2, 4),
    (1, 1, 3, 3),
    (1, 1, 4, 4),
    (1, 2, 2, 3),
    (1, 2, 3, 2),
    (1, 2, 3, 4),
)


def tau_position(position: tuple[int, int]) -> tuple[int, int]:
    row, column = position
    return 4 - column, 4 - row


def tau_even_directions(phase: int, positive: bool) -> list[Direction]:
    positions = [
        (
            (index + phase, index)
            if positive
            else (index, index + phase)
        )
        for index in range(5 - phase)
    ]
    seen: set[int] = set()
    result = []
    for index, position in enumerate(positions):
        if index in seen:
            continue
        partner = 4 - phase - index
        seen.update((index, partner))
        direction = {position: 1}
        if partner != index:
            direction[positions[partner]] = (-1) ** phase
        for target, coefficient in direction.items():
            source = tau_position(target)
            expected = (-1) ** sum(target) * direction.get(source, 0)
            assert coefficient == expected
        result.append(direction)
    return result


def chart_moments(cutoff: int) -> list[Polynomial]:
    positive = tau_even_directions(1, True)
    negative = tau_even_directions(1, False)
    assert len(positive) == len(negative) == 2

    # A term is (row contribution, column contribution, chart variable,
    # coefficient).  Variable None denotes the normalized coefficient of B0.
    terms: list[tuple[int, int, int | None, int]] = [
        (index, index, index, 1) for index in range(5)
    ]
    terms.extend(
        (row, column, None, coefficient)
        for (row, column), coefficient in positive[0].items()
    )
    terms.extend(
        (row, column, 5, coefficient)
        for (row, column), coefficient in positive[1].items()
    )
    terms.extend(
        (row, column, 6, coefficient)
        for (row, column), coefficient in negative[0].items()
    )
    terms.extend(
        (row, column, 7, coefficient)
        for (row, column), coefficient in negative[1].items()
    )

    state: dict[tuple[int, int, Exponent], int] = {
        (0, 0, (0,) * len(VARIABLES)): 1
    }
    moments = []
    factorials = [1]
    for value in range(1, 4 * cutoff + 1):
        factorials.append(factorials[-1] * value)

    for order in range(1, cutoff + 1):
        updated: dict[tuple[int, int, Exponent], int] = {}
        for (left, right, exponents), value in state.items():
            for delta_left, delta_right, variable, coefficient in terms:
                new_exponents = list(exponents)
                if variable is not None:
                    new_exponents[variable] += 1
                key = (
                    left + delta_left,
                    right + delta_right,
                    tuple(new_exponents),
                )
                updated[key] = updated.get(key, 0) + value * coefficient
        state = {key: value for key, value in updated.items() if value}

        moment: Polynomial = {}
        for (left, right, exponents), value in state.items():
            if left != right:
                continue
            moment[exponents] = (
                moment.get(exponents, 0)
                + factorials[left] * factorials[4 * order - left] * value
            )
        moments.append(
            {
                exponents: value
                for exponents, value in moment.items()
                if value
            }
        )
    return moments


def evaluate(polynomial: Polynomial, point: tuple[int, ...]) -> int:
    return sum(
        coefficient
        * reduce(
            mul,
            (
                value**exponent
                for value, exponent in zip(point, exponents, strict=True)
            ),
            1,
        )
        for exponents, coefficient in polynomial.items()
    )


def jacobian_determinant(
    polynomials: list[Polynomial],
    point: tuple[object, ...],
) -> sp.Expr:
    rows = []
    for polynomial in polynomials:
        row = []
        for differentiated in range(len(VARIABLES)):
            value = 0
            for exponents, coefficient in polynomial.items():
                exponent = exponents[differentiated]
                if not exponent:
                    continue
                term = coefficient * exponent
                for index, power in enumerate(exponents):
                    term *= point[index] ** (
                        power - (1 if index == differentiated else 0)
                    )
                value += term
            row.append(value)
        rows.append(row)
    return sp.det(sp.Matrix(rows))


def polynomial_string(
    polynomial: Polynomial,
    prime: int,
    variables: tuple[str, ...] = VARIABLES,
) -> str:
    terms = []
    for exponents, coefficient in polynomial.items():
        reduced = coefficient if prime == 0 else coefficient % prime
        if reduced == 0:
            continue
        negative = reduced < 0
        absolute = abs(reduced)
        factors = []
        if absolute != 1 or not any(exponents):
            factors.append(str(absolute))
        for variable, exponent in zip(variables, exponents, strict=True):
            if exponent == 1:
                factors.append(variable)
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        term = "*".join(factors) if factors else "0"
        terms.append(("-" if negative else "+") + term)
    if not terms:
        return "0"
    return "".join(terms).lstrip("+")


def primitive_polynomial(polynomial: Polynomial) -> Polynomial:
    content = reduce(
        gcd,
        (abs(coefficient) for coefficient in polynomial.values()),
        0,
    )
    if content <= 1:
        return polynomial
    return {
        exponents: coefficient // content
        for exponents, coefficient in polynomial.items()
    }


def matrix_multiply(
    left: list[list[int]],
    right: list[list[int]],
    prime: int,
) -> list[list[int]]:
    return [
        [
            sum(
                left[row][middle] * right[middle][column]
                for middle in range(5)
            )
            % prime
            for column in range(5)
        ]
        for row in range(5)
    ]


def trace_word(
    word: tuple[int, ...],
    components: dict[int, list[list[int]]],
    prime: int,
) -> int:
    product = components[word[0]]
    for component in word[1:]:
        product = matrix_multiply(product, components[component], prime)
    return sum(product[index][index] for index in range(5)) % prime


def chart_matrix(point: tuple[int, ...], prime: int) -> list[list[int]]:
    positive = tau_even_directions(1, True)
    negative = tau_even_directions(1, False)
    matrix = [[0] * 5 for _ in range(5)]
    for index, value in enumerate(point[:5]):
        matrix[index][index] = value if prime == 0 else value % prime
    for direction, coefficient in (
        (positive[0], 1),
        (positive[1], point[5]),
        (negative[0], point[6]),
        (negative[1], point[7]),
    ):
        for (row, column), value in direction.items():
            updated = matrix[row][column] + coefficient * value
            matrix[row][column] = (
                updated if prime == 0 else updated % prime
            )
    return matrix


def operator_matrix(point: tuple[int, ...], prime: int) -> sp.Matrix:
    coefficients = sp.Matrix(chart_matrix(point, prime))
    diagonal = sp.diag(
        *[
            factorial(index) * factorial(4 - index)
            for index in range(5)
        ]
    )
    return coefficients.T * diagonal


def symmetric_fourth_matrix(
    group_variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    a, b, c, d = group_variables
    x, y = sp.symbols("x y")
    result = sp.zeros(5)
    for column in range(5):
        image = sp.Poly(
            (a * x + c * y) ** (4 - column)
            * (b * x + d * y) ** column,
            x,
            y,
        )
        for row in range(5):
            result[row, column] = image.coeff_monomial(
                x ** (4 - row) * y**row
            )
    return result


def sympy_polynomial_string(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
) -> str:
    polynomial = sp.Poly(expression, *variables)
    terms = []
    for exponents, coefficient in polynomial.terms():
        residue = int(coefficient) % prime
        if not residue:
            continue
        factors = []
        if residue != 1 or not any(exponents):
            factors.append(str(residue))
        for variable, exponent in zip(variables, exponents, strict=True):
            if exponent == 1:
                factors.append(str(variable))
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        terms.append("*".join(factors) if factors else "0")
    return "+".join(terms) if terms else "0"


def orbit_witness_equations(
    source: tuple[int, ...],
    target: tuple[int, ...],
    prime: int,
) -> tuple[tuple[sp.Symbol, ...], list[str]]:
    group_variables = sp.symbols("ga gb gc gd")
    ga, gb, gc, gd = group_variables
    representation = symmetric_fourth_matrix(group_variables)
    source_operator = operator_matrix(source, prime)
    target_operator = operator_matrix(target, prime)
    intertwining = sp.expand(
        representation * source_operator
        - target_operator * representation
    )
    expressions = [ga * gd - gb * gc - 1]
    expressions.extend(
        intertwining[row, column]
        for row in range(5)
        for column in range(5)
        if intertwining[row, column] != 0
    )
    return group_variables, [
        sympy_polynomial_string(expression, group_variables, prime)
        for expression in expressions
    ]


def component_matrices(
    point: tuple[int, ...],
    prime: int,
) -> dict[int, list[list[int]]]:
    projectors = invariant_base.casimir_projectors()
    matrix = chart_matrix(point, prime)
    diagonal = [
        factorial(index) * factorial(4 - index)
        for index in range(5)
    ]
    vector = [
        matrix[column][row] * diagonal[column] % prime
        for column in range(5)
        for row in range(5)
    ]
    result = {}
    for component, projector in projectors.items():
        projected = []
        for row in projector.tolist():
            value = 0
            for entry, coordinate in zip(row, vector, strict=True):
                numerator, denominator = sp.fraction(entry)
                residue = (
                    int(numerator)
                    * pow(int(denominator) % prime, -1, prime)
                    % prime
                )
                value += residue * coordinate
            projected.append(value % prime)
        result[component] = [
            [projected[row + 5 * column] for column in range(5)]
            for row in range(5)
        ]
    return result


def tau_even_parameter_values(
    point: tuple[int, ...],
    prime: int,
) -> list[tuple[str, int]]:
    components = component_matrices(point, prime)
    values = [("tr(A_0)", trace_word((0,), components, prime))]
    values.extend(
        (
            f"tr(A_{2 * component}^2)",
            trace_word((component, component), components, prime),
        )
        for component in range(1, 5)
    )
    for triple in EVEN_CUBIC_TRIPLES:
        permutations = sorted(set(itertools.permutations(triple)))
        values.append(
            (
                "symtr(" + ",".join(map(str, triple)) + ")",
                sum(
                    trace_word(word, components, prime)
                    for word in permutations
                )
                % prime,
            )
        )
    values.extend(
        (
            "tr(" + "".join(map(str, word)) + ")",
            trace_word(word, components, prime),
        )
        for word in EVEN_QUARTIC_WORDS
    )
    assert len(values) == 22
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument(
        "--orders",
        type=int,
        nargs="+",
        default=list(range(1, 10)),
        help="one-based moment orders used in the fiber ideal",
    )
    parser.add_argument(
        "--point",
        type=int,
        nargs=8,
        default=DEFAULT_POINT,
        metavar=VARIABLES,
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument(
        "--linear-algebra",
        type=int,
        choices=(1, 2, 42, 44),
        default=2,
    )
    parser.add_argument(
        "--lifting-mulmat",
        type=int,
        choices=(0, 1),
        default=0,
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="truncate large msolve basis output to its header and tail",
    )
    parser.add_argument(
        "--leading-homogeneous",
        action="store_true",
        help=(
            "replace each fiber equation by its top standard-degree "
            "homogeneous part to test the boundary at infinity"
        ),
    )
    parser.add_argument(
        "--projective-boundary",
        action="store_true",
        help=(
            "homogenize the fiber equations, saturate by the "
            "homogenizing variable in Singular, and inspect the boundary"
        ),
    )
    parser.add_argument(
        "--groebner-basis",
        type=int,
        choices=(0, 1, 2),
        default=1,
    )
    parser.add_argument(
        "--parametrization",
        type=int,
        choices=(0, 1),
        default=0,
    )
    parser.add_argument(
        "--compare-even-parameters",
        type=int,
        nargs=8,
        metavar=VARIABLES,
        help=(
            "compare the 22 known tau-even trace parameters at this "
            "second chart point"
        ),
    )
    parser.add_argument(
        "--test-apolar-orbit",
        type=int,
        nargs=8,
        metavar=VARIABLES,
        help=(
            "test whether this chart point lies in the SL2 orbit of the "
            "apolar reversal of --point"
        ),
    )
    parser.add_argument(
        "--certify-example",
        action="store_true",
        help=(
            "verify and record the fixed F_101 first-ten-moment example"
        ),
    )
    parser.add_argument(
        "--centered-singular",
        action="store_true",
        help=(
            "run a characteristic-zero Singular basis after the affine "
            "change centered on the multimodularly reconstructed fiber"
        ),
    )
    parser.add_argument(
        "--singular-algorithm",
        choices=("liftstd", "modstd", "slimgb", "std"),
        default="modstd",
        help="standard-basis algorithm used with --centered-singular",
    )
    parser.add_argument(
        "--singular-task",
        choices=("basis", "residual-saturation"),
        default="basis",
    )
    parser.add_argument(
        "--modstd-exactness",
        type=int,
        choices=(0, 1),
        default=1,
    )
    parser.add_argument(
        "--verify-lift",
        action="store_true",
        help=(
            "after reconstructing a centered basis, verify forward "
            "reduction and exact reverse containment with lift(I,G)"
        ),
    )
    parser.add_argument(
        "--radical-generator",
        choices=(
            "x0",
            "x1",
            "x2",
            "x3",
            "t",
            "y",
            "r_quadratic",
            "s_quadratic",
        ),
        help=(
            "in centered Singular mode, adjoin inverse*g-1 for the "
            "selected predicted-fiber generator"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    executable = shutil.which("msolve")
    if executable is None:
        raise SystemExit("msolve is required on PATH")
    if min(args.orders) < 1:
        raise SystemExit("moment orders are one-based")

    point = tuple(args.point)
    comparison_payload = None
    if args.compare_even_parameters is not None:
        if args.prime == 0:
            raise SystemExit(
                "--compare-even-parameters currently requires a prime field"
            )
        comparison = tuple(args.compare_even_parameters)
        left = tau_even_parameter_values(point, args.prime)
        right = tau_even_parameter_values(comparison, args.prime)
        differences = [
            {
                "label": left_value[0],
                "left": left_value[1],
                "right": right_value[1],
            }
            for left_value, right_value in zip(left, right, strict=True)
            if left_value[1] != right_value[1]
        ]
        left_components = component_matrices(point, args.prime)
        right_components = component_matrices(comparison, args.prime)
        left_odd = trace_word((2, 3, 4), left_components, args.prime)
        right_odd = trace_word((2, 3, 4), right_components, args.prime)
        print("tau-even parameter differences:", differences)
        print("odd cubic values:", left_odd, right_odd)
        comparison_payload = {
            "point": list(comparison),
            "tau_even_parameter_differences": differences,
            "odd_cubic_values": [left_odd, right_odd],
        }

    cutoff = max(args.orders)
    moments = chart_moments(cutoff)
    targets = [evaluate(moment, point) for moment in moments]
    equations = []
    equation_polynomials = []
    equation_contents = []
    for order in args.orders:
        if args.leading_homogeneous:
            polynomial = {
                exponents: coefficient
                for exponents, coefficient in moments[order - 1].items()
                if sum(exponents) == order
            }
        else:
            polynomial = dict(moments[order - 1])
            zero = (0,) * len(VARIABLES)
            polynomial[zero] = (
                polynomial.get(zero, 0) - targets[order - 1]
            )
        content = reduce(
            gcd,
            (abs(coefficient) for coefficient in polynomial.values()),
            0,
        )
        equation_contents.append(content)
        if args.prime == 0:
            polynomial = primitive_polynomial(polynomial)
        equation_polynomials.append(polynomial)
        equations.append(polynomial_string(polynomial, args.prime))
    print("equation contents:", equation_contents)

    if args.projective_boundary:
        singular = shutil.which("Singular")
        if singular is None:
            raise SystemExit("Singular is required on PATH")
        homogenized_variables = VARIABLES + ("H",)
        homogenized_equations = []
        for order, polynomial in zip(
            args.orders,
            equation_polynomials,
            strict=True,
        ):
            homogenized: Polynomial = {}
            for exponents, coefficient in polynomial.items():
                homogenizing_exponent = order - sum(exponents)
                assert homogenizing_exponent >= 0
                new_exponents = exponents + (homogenizing_exponent,)
                homogenized[new_exponents] = coefficient
            homogenized_equations.append(
                polynomial_string(
                    homogenized,
                    args.prime,
                    homogenized_variables,
                )
            )
        code = (
            'LIB "elim.lib";\n'
            + f"ring rb={args.prime},("
            + ",".join(homogenized_variables)
            + "),dp;\n"
            + "ideal I="
            + ",\n".join(homogenized_equations)
            + ";\n"
            + "list saturation_data=sat(I,ideal(H));\n"
            + "ideal S=saturation_data[1];\n"
            + "ideal B=S+ideal(H);\n"
            + "ideal G=std(B);\n"
            + 'print("BOUNDARY_DIM");\n'
            + "dim(G);\n"
            + 'print("BOUNDARY_SIZE");\n'
            + "size(G);\n"
            + 'print("BOUNDARY_BASIS");\n'
            + "G;\n"
        )
        completed = subprocess.run(
            [singular, "-q"],
            input=code,
            text=True,
            capture_output=True,
            timeout=args.timeout,
            check=False,
        )
        print(
            "saturated projective boundary:",
            f"prime={args.prime}",
            f"orders={args.orders}",
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        print("Singular exit code:", completed.returncode)
        return

    if args.centered_singular:
        singular = shutil.which("Singular")
        if singular is None:
            raise SystemExit("Singular is required on PATH")
        replacements = {
            "a0": "(x0-s+13)",
            "a1": "((x1-4*s+71)/9)",
            "a2": "(x2+5)",
            "a3": "((x3+4*s+19)/9)",
            "a4": "(s)",
            "u": "((t-r+11)/3)",
            "v": "(y+3)",
            "w": "(r)",
        }
        centered_equations = []
        for equation in equations:
            centered = equation
            for variable in VARIABLES:
                centered = centered.replace(variable, replacements[variable])
            centered_equations.append(centered)
        transformation_report = ""
        if args.singular_task == "residual-saturation":
            library = 'LIB "elim.lib";\n'
            basis_setup = (
                "ideal J=x0,x1,x2,x3,t,y,"
                "(r-5)*(r-6),(s-2)*(s-11);\n"
                "list saturation_data=sat(I,J);\n"
                "ideal G=saturation_data[1];\n"
            )
        elif args.singular_algorithm == "liftstd":
            library = ""
            basis_setup = (
                "matrix transformation;\n"
                'ideal G=liftstd(I,transformation,"slimgb");\n'
            )
            transformation_report = (
                'print("TRANSFORMATION_ROWS");\n'
                "nrows(transformation);\n"
                'print("TRANSFORMATION_COLUMNS");\n'
                "ncols(transformation);\n"
                'print("TRANSFORMATION_TERM_COUNTS");\n'
                "int tr;\n"
                "int tc;\n"
                "for (tr=1;tr<=nrows(transformation);tr++)\n"
                "{\n"
                "  for (tc=1;tc<=ncols(transformation);tc++)\n"
                "  {\n"
                "    size(transformation[tr,tc]);\n"
                "  }\n"
                "}\n"
            )
        elif args.singular_algorithm == "modstd":
            library = 'LIB "modstd.lib";\n'
            basis_setup = (
                f"ideal G=modStd(I,{args.modstd_exactness});\n"
            )
        else:
            library = ""
            basis_setup = (
                f"ideal G={args.singular_algorithm}(I);\n"
            )
        lift_verification = ""
        if args.verify_lift:
            lift_verification = (
                'print("FORWARD_REDUCTION_SIZE");\n'
                "ideal verified_basis=std(G);\n"
                "size(reduce(I,verified_basis));\n"
                "matrix reverse_lift=lift(I,G);\n"
                "int lift_row;\n"
                "int lift_column;\n"
                "poly reconstructed_generator;\n"
                "for (lift_column=1;lift_column<=size(G);lift_column++)\n"
                "{\n"
                "  reconstructed_generator=0;\n"
                "  for (lift_row=1;lift_row<=size(I);lift_row++)\n"
                "  {\n"
                "    reconstructed_generator="
                "reconstructed_generator"
                "+I[lift_row]*reverse_lift[lift_row,lift_column];\n"
                "  }\n"
                "  if (reconstructed_generator-G[lift_column]!=0)\n"
                "  {\n"
                '    print("BAD_REVERSE_LIFT");\n'
                "    exit(1);\n"
                "  }\n"
                "}\n"
                'print("PASS_EXACT_REVERSE_LIFT");\n'
            )
        centered_ring_variables = (
            "x0,x1,x2,x3,t,y,r,s"
            if args.radical_generator is None
            else "inverse,x0,x1,x2,x3,t,y,r,s"
        )
        centered_ideal_entries = list(centered_equations)
        if args.radical_generator is not None:
            radical_generators = {
                "x0": "x0",
                "x1": "x1",
                "x2": "x2",
                "x3": "x3",
                "t": "t",
                "y": "y",
                "r_quadratic": "(r-5)*(r-6)",
                "s_quadratic": "(s-2)*(s-11)",
            }
            selected_generator = radical_generators[
                args.radical_generator
            ]
            centered_ideal_entries.append(
                f"inverse*({selected_generator})-1"
            )
        code = (
            library
            + f"ring rr={args.prime},("
            + centered_ring_variables
            + "),dp;\n"
            + "ideal I="
            + ",\n".join(centered_ideal_entries)
            + ";\n"
            + basis_setup
            + 'print("CENTERED_DIM");\n'
            + "dim(G);\n"
            + 'print("CENTERED_VDIM");\n'
            + "vdim(G);\n"
            + 'print("CENTERED_SIZE");\n'
            + "size(G);\n"
            + transformation_report
            + lift_verification
            + 'print("CENTERED_BASIS");\n'
            + "G;\n"
        )
        completed = subprocess.run(
            [singular, "-q"],
            input=code,
            text=True,
            capture_output=True,
            timeout=args.timeout,
            check=False,
        )
        print(
            "centered characteristic-zero chart:",
            f"orders={args.orders}",
            f"point={point}",
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        print("Singular exit code:", completed.returncode)
        return

    with tempfile.TemporaryDirectory(
        prefix="degree-four-phase-one-",
        dir=ROOT / "artifacts" / "generated-results",
    ) as temporary:
        input_path = Path(temporary) / "fiber.ms"
        output_path = Path(temporary) / "fiber.out"
        msolve_variables = VARIABLES
        msolve_equations = list(equations)
        if args.radical_generator is not None:
            original_radical_generators = {
                "x0": "a0+a4-13",
                "x1": "9*a1+4*a4-71",
                "x2": "a2-5",
                "x3": "9*a3-4*a4-19",
                "t": "3*u+w-11",
                "y": "v-3",
                "r_quadratic": "(w-5)*(w-6)",
                "s_quadratic": "(a4-2)*(a4-11)",
            }
            msolve_variables = VARIABLES + ("inverse",)
            msolve_equations.append(
                "inverse*("
                + original_radical_generators[args.radical_generator]
                + ")-1"
            )
        input_path.write_text(
            ",".join(msolve_variables)
            + "\n"
            + str(args.prime)
            + "\n"
            + ",\n".join(msolve_equations)
            + "\n"
        )
        command = [
            executable,
            "-f",
            str(input_path),
            "-o",
            str(output_path),
            "-t",
            str(args.threads),
            "-v",
            "1",
            "-l",
            str(args.linear_algebra),
            "-L",
            str(args.lifting_mulmat),
            "-g",
            str(args.groebner_basis),
            "-P",
            str(args.parametrization),
            "--random-seed",
            "1",
        ]
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=args.timeout,
            check=False,
        )
        print(
            "phase-one chart:",
            f"prime={args.prime}",
            f"orders={args.orders}",
            f"point={point}",
        )
        print(
            "moment term counts:",
            [len(moments[order - 1]) for order in args.orders],
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        print("msolve exit code:", completed.returncode)
        if output_path.exists():
            print("msolve output:")
            fiber_text = output_path.read_text().strip()
            if args.summary_only and len(fiber_text.splitlines()) > 60:
                fiber_lines = fiber_text.splitlines()
                print("\n".join(fiber_lines[:35]))
                print(
                    f"... {len(fiber_lines) - 55} output lines omitted ..."
                )
                print("\n".join(fiber_lines[-20:]))
            else:
                print(fiber_text)
        else:
            fiber_text = ""

        orbit_text = ""
        if args.test_apolar_orbit is not None:
            if args.prime == 0:
                raise SystemExit(
                    "--test-apolar-orbit currently requires a prime field"
                )
            reversed_point = point[4::-1] + point[5:]
            orbit_target = tuple(args.test_apolar_orbit)
            group_variables, orbit_equations = orbit_witness_equations(
                reversed_point,
                orbit_target,
                args.prime,
            )
            orbit_input = Path(temporary) / "orbit.ms"
            orbit_output = Path(temporary) / "orbit.out"
            orbit_input.write_text(
                ",".join(map(str, group_variables))
                + "\n"
                + str(args.prime)
                + "\n"
                + ",\n".join(orbit_equations)
                + "\n"
            )
            orbit_command = [
                executable,
                "-f",
                str(orbit_input),
                "-o",
                str(orbit_output),
                "-t",
                str(args.threads),
                "-v",
                "1",
                "-g",
                "2",
                "--random-seed",
                "1",
            ]
            orbit_completed = subprocess.run(
                orbit_command,
                text=True,
                capture_output=True,
                timeout=args.timeout,
                check=False,
            )
            print(
                "apolar-orbit test:",
                f"source={reversed_point}",
                f"target={orbit_target}",
            )
            print(orbit_completed.stdout.strip())
            if orbit_completed.stderr.strip():
                print(orbit_completed.stderr.strip())
            print("orbit msolve exit code:", orbit_completed.returncode)
            if orbit_output.exists():
                print("orbit msolve output:")
                orbit_text = orbit_output.read_text().strip()
                print(orbit_text)

        if args.certify_example:
            expected_comparison = (2, 3, 5, 7, 11, 69, 3, 6)
            assert args.prime == 101
            assert args.orders == list(range(1, 11))
            assert point == DEFAULT_POINT
            assert tuple(args.compare_even_parameters or ()) == (
                expected_comparison
            )
            assert tuple(args.test_apolar_orbit or ()) == expected_comparison
            assert completed.returncode == 0
            expected_fiber_basis = [
                "1*v^1+98",
                "1*u^1+34*w^1+30",
                "1*a3^1+22*a4^1+54",
                "1*a2^1+96",
                "1*a1^1+79*a4^1+37",
                "1*a0^1+1*a4^1+88",
                "1*w^2+90*w^1+30",
                "1*a4^2+88*a4^1+22",
            ]
            assert all(entry in fiber_text for entry in expected_fiber_basis)
            assert comparison_payload is not None
            assert not comparison_payload["tau_even_parameter_differences"]
            assert comparison_payload["odd_cubic_values"] == [11, 90]
            expected_orbit_basis = [
                "1*gd^1",
                "1*gb^1+34*gc^1",
                "1*ga^1",
                "1*gc^2+98",
            ]
            assert all(entry in orbit_text for entry in expected_orbit_basis)

            moments_through_eleven = chart_moments(11)
            moment_differences = [
                (
                    evaluate(moment, point)
                    - evaluate(moment, expected_comparison)
                )
                % 101
                for moment in moments_through_eleven
            ]
            assert moment_differences == [0] * 11

            rational_branch = (
                2,
                3,
                5,
                7,
                11,
                Fraction(5, 3),
                3,
                6,
            )
            assert all(
                evaluate(moment, point)
                == evaluate(moment, rational_branch)
                for moment in moments_through_eleven
            )
            projectors = invariant_base.casimir_projectors()

            def exact_odd_cubic(branch: tuple[object, ...]) -> sp.Expr:
                coefficients = sp.Matrix(chart_matrix(branch, 0))
                _, vector = invariant_base.vectorized_operator(coefficients)
                components = invariant_base.component_matrices(
                    vector,
                    projectors,
                )
                return invariant_base.cubic_invariant(
                    components,
                    (2, 3, 4),
                )

            assert exact_odd_cubic(point) == 1728
            assert exact_odd_cubic(rational_branch) == -1728

            rational_points = [
                tuple(point),
                rational_branch,
                point[4::-1] + point[5:],
                rational_branch[4::-1] + rational_branch[5:],
            ]
            jacobian_determinants = [
                jacobian_determinant(
                    moments_through_eleven[:8],
                    branch,
                )
                for branch in rational_points
            ]
            jacobian_absolute_value = int(
                abs(jacobian_determinants[0])
            )
            assert jacobian_absolute_value == int(
                "664886859697950537008868405994595186881355199398116360"
                "363695342053926137270301114156567309246654270561514082"
                "878820219617280000000000000000000000"
            )
            assert jacobian_determinants == [
                -jacobian_absolute_value,
                jacobian_absolute_value,
                -jacobian_absolute_value,
                jacobian_absolute_value,
            ]

            square_root_three = sp.sqrt(3)
            exact_group = (
                sp.Integer(0),
                -square_root_three / 3,
                square_root_three,
                sp.Integer(0),
            )
            representation = symmetric_fourth_matrix(exact_group)
            exact_reversal = point[4::-1] + point[5:]
            intertwining = sp.simplify(
                representation * operator_matrix(exact_reversal, 0)
                - operator_matrix(rational_branch, 0) * representation
            )
            assert intertwining == sp.zeros(5)
            assert sp.simplify(
                exact_group[0] * exact_group[3]
                - exact_group[1] * exact_group[2]
            ) == 1

            payload = {
                "format": "degree-four-phase-one-chart-modular-v3",
                "status": (
                    "exact modular fiber completeness plus exact "
                    "characteristic-zero branch and orbit identities"
                ),
                "prime": 101,
                "chart_coordinates": list(VARIABLES),
                "point": list(point),
                "moment_orders_in_fiber_ideal": list(range(1, 11)),
                "reduced_groebner_basis": expected_fiber_basis,
                "fiber_length": 4,
                "fiber_reduced": True,
                "second_branch_same_diagonal": list(expected_comparison),
                "moments_one_through_eleven_agree_on_branches": True,
                "known_tau_even_parameters_agree_on_branches": True,
                "odd_cubic_values_on_branches": [11, 90],
                "characteristic_zero_second_branch": [
                    "2",
                    "3",
                    "5",
                    "7",
                    "11",
                    "5/3",
                    "3",
                    "6",
                ],
                "characteristic_zero_moments_one_through_eleven_agree": (
                    True
                ),
                "characteristic_zero_odd_cubic_values": [1728, -1728],
                "characteristic_zero_sl2_orbit_matrix": [
                    ["0", "-sqrt(3)/3"],
                    ["sqrt(3)", "0"],
                ],
                "characteristic_zero_orbit_identity_verified": True,
                "characteristic_zero_reconstructed_points": [
                    [str(coordinate) for coordinate in branch]
                    for branch in rational_points
                ],
                "first_eight_moment_jacobian_determinants": [
                    str(value) for value in jacobian_determinants
                ],
                "characteristic_zero_reconstructed_points_reduced_and_isolated": (
                    True
                ),
                "apolar_reversal_of_point": list(
                    point[4::-1] + point[5:]
                ),
                "sl2_orbit_witness_basis": expected_orbit_basis,
                "orbit_witness_interpretation": (
                    "ga=gd=0, gb=-34*gc, gc^2=3; the second branch "
                    "is SL2-conjugate to the apolar reversal"
                ),
                "quotient_interpretation": (
                    "the four raw chart points pair into two candidate "
                    "SL2 quotient points exchanged by the apolar involution"
                ),
                "scope_warning": (
                    "fiber completeness is proved only over F_101; the "
                    "displayed characteristic-zero branches and their "
                    "orbit identity are exact, but exclusion of further "
                    "characteristic-zero points and the generic degree "
                    "remain open"
                ),
            }
            OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
            print("PASS certified modular phase-one chart example")


if __name__ == "__main__":
    main()

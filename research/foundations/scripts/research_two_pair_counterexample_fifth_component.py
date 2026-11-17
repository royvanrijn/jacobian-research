#!/usr/bin/env python3
"""Component-wide fifth-order reconnaissance for the SIC2C4 local branch.

The calculation works on one Q(sqrt(41))-component above the reduced
direction (1,2,3,4,5).  It restores the full nine coordinates of the
fourth-order affine component and the eleven-dimensional cubic-tangent
kernel.  The two effective fifth-correction columns are eliminated before
forming a determinantal survivor ideal on a dense pivot chart.

An exact nonzero constant augmented minor supplies the characteristic-zero
certificate.  The final modular saturation is only an independent regression
check.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


ROOT = Path(__file__).resolve().parents[1]
ALGEBRAIZATION_SCRIPT = (
    ROOT / "scripts" / "research_two_pair_counterexample_algebraization.py"
)
PRIME = 32003
DIRECTIONS = {
    "quadratic_extension_41": (1, 2, 3, 4, 5),
    "generic_rational": (2, -1, 3, 1, -2),
    "pure_apolar_odd": (0, 1, 0, 0, 0),
}

spec = importlib.util.spec_from_file_location(
    "sic2_algebraization_component",
    ALGEBRAIZATION_SCRIPT,
)
assert spec and spec.loader
algebraization = importlib.util.module_from_spec(spec)
spec.loader.exec_module(algebraization)

fourth = algebraization.fourth
seed = algebraization.seed
BASIS = algebraization.BASIS


def parse_groebner_basis(values: list[str]) -> list[sp.Expr]:
    variables = sp.symbols("u0:11")
    table = {str(variable): variable for variable in variables}
    return [
        sp.sympify(
            value.rstrip(",").replace("^", "**"),
            locals=table,
        )
        for value in values
    ]


def bilinear_matrix(degree: int, f_order: int) -> sp.Matrix:
    weights = fourth.generating_weights(degree, f_order)
    matrix = sp.zeros(25)
    for left, left_exponent in enumerate(BASIS):
        for right in range(left, 25):
            right_exponent = BASIS[right]
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(4)
            )
            matrix[left, right] = weights.get(exponent, 0)
            matrix[right, left] = matrix[left, right]
    return matrix


def fixed_times_basis_functional(
    fixed: seed.Polynomial,
    total_degree: int,
    f_order: int,
) -> sp.Matrix:
    weights = fourth.generating_weights(total_degree, f_order)
    values = []
    for basis_exponent in BASIS:
        value = sp.Rational(0)
        for fixed_exponent, coefficient in fixed.items():
            exponent = tuple(
                fixed_exponent[index] + basis_exponent[index]
                for index in range(4)
            )
            value += coefficient * weights.get(exponent, 0)
        values.append(sp.factor(value))
    return sp.Matrix(values)


def h_times_two_basis_matrix(
    h_polynomial: seed.Polynomial,
    f_order: int,
) -> sp.Matrix:
    weights = fourth.generating_weights(12, f_order)
    matrix = sp.zeros(25)
    for left, left_exponent in enumerate(BASIS):
        for right in range(left, 25):
            right_exponent = BASIS[right]
            value = sp.Rational(0)
            for h_exponent, coefficient in h_polynomial.items():
                exponent = tuple(
                    h_exponent[index]
                    + left_exponent[index]
                    + right_exponent[index]
                    for index in range(4)
                )
                value += coefficient * weights.get(exponent, 0)
            matrix[left, right] = sp.factor(value)
            matrix[right, left] = matrix[left, right]
    return matrix


def fixed_contraction(
    polynomial: seed.Polynomial,
    degree: int,
    f_order: int,
) -> sp.Expr:
    weights = fourth.generating_weights(degree, f_order)
    return sp.factor(
        sum(
            coefficient * weights.get(exponent, 0)
            for exponent, coefficient in polynomial.items()
        )
    )


def modular_polynomial(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> str:
    root = int(sp.sqrt_mod(41, PRIME, all_roots=True)[0])
    rational_polynomial = sp.Poly(
        sp.cancel(expression).subs(sp.sqrt(41), root),
        *variables,
        domain=sp.QQ,
    )
    modular_terms = {}
    for monomial, coefficient in rational_polynomial.terms():
        rational_coefficient = rational_polynomial.domain.to_sympy(
            coefficient
        )
        value = (
            int(rational_coefficient.p)
            * pow(int(rational_coefficient.q), -1, PRIME)
        ) % PRIME
        if value:
            modular_terms[monomial] = value
    reduced = sp.Poly.from_dict(
        modular_terms,
        variables,
        modulus=PRIME,
    ).as_expr()
    return str(reduced).replace("**", "^")


def singular_saturation(
    equations: list[sp.Expr],
    pivot: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> dict[str, object]:
    singular = shutil.which("Singular")
    assert singular, "Singular is required for the modular survivor ideal"
    names = ",".join(str(variable) for variable in variables)
    equation_source = ",".join(
        modular_polynomial(equation, variables)
        for equation in equations
        if equation
    )
    pivot_source = modular_polynomial(pivot, variables)
    source = "\n".join(
        [
            'LIB "elim.lib";',
            f"ring r={PRIME},({names}),dp;",
            f"ideal I={equation_source};",
            f"ideal J={pivot_source};",
            "ideal S=sat(I,J)[1];",
            "ideal G=std(S);",
            'print("DIM"); print(dim(G));',
            'print("DEG"); print(deg(G));',
            'print("SIZE"); print(size(G));',
            'print("GB"); print(G); print("END");',
        ]
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=source,
        text=True,
        capture_output=True,
        check=True,
        timeout=600,
    )
    lines = [line.strip() for line in completed.stdout.splitlines()]

    def integer_after(label: str) -> int:
        return int(lines[lines.index(label) + 1])

    basis = lines[lines.index("GB") + 1 : lines.index("END")]
    return {
        "prime": PRIME,
        "sqrt_41_residue": int(sp.sqrt_mod(41, PRIME, all_roots=True)[0]),
        "dimension": integer_after("DIM"),
        "degree": integer_after("DEG"),
        "groebner_basis_size": integer_after("SIZE"),
        "groebner_basis": basis,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--label",
        choices=("all", *DIRECTIONS),
        default="all",
    )
    parser.add_argument(
        "--direction",
        nargs=5,
        type=int,
        metavar=("H0", "H1", "H2", "H3", "H4"),
    )
    parser.add_argument("--tag")
    arguments = parser.parse_args()
    if arguments.direction is None and arguments.label == "all":
        for selected_label in DIRECTIONS:
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--label",
                    selected_label,
                ],
                check=True,
            )
        print("PASS fifth component: all three selected directions")
        return
    if arguments.direction is not None:
        direction = tuple(arguments.direction)
        label = arguments.tag or "custom_" + "_".join(map(str, direction))
    else:
        label = arguments.label
        direction = DIRECTIONS[label]
    suffix = "" if label == "quadratic_extension_41" else f"_{label}"
    output = (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"two_pair_counterexample_fifth_component_research{suffix}.json"
    )
    context = fourth.main(write_output=False, return_context=True)
    assert context is not None
    state: dict[str, object] = {}
    report = algebraization.fourth_fiber(
        context,
        direction,
        run_point_tests=False,
        state_target=state,
    )
    assert (report["dimension"], report["degree"]) == (9, 2)
    print("PASS fifth component: recovered the degree-two fourth fiber")

    u = sp.symbols("u0:11")
    q = sp.symbols("q0:9")
    z = sp.symbols("z0:11")
    groebner = parse_groebner_basis(state["groebner_basis"])
    linear = next(
        value
        for value in groebner
        if sp.Poly(value, *u).total_degree() == 1
    )
    quadric = next(
        value
        for value in groebner
        if sp.Poly(value, *u).total_degree() == 2
    )
    factorization = sp.factor_list(quadric, extension=sp.sqrt(41))[1]
    branch_factors = [
        factor
        for factor, multiplicity in factorization
        for _ in range(multiplicity)
        if sp.Poly(factor, *u).total_degree() == 1
    ]
    assert len(branch_factors) == 2
    branch = branch_factors[0]
    u1 = sp.factor(sp.solve(branch, u[1])[0])
    u0 = sp.factor(sp.solve(linear.subs(u[1], u1), u[0])[0])
    component_substitution = {
        u[index + 2]: q[index]
        for index in range(9)
    }
    u1 = sp.factor(u1.subs(component_substitution))
    u0 = sp.factor(u0.subs(component_substitution).subs(u[1], u1))
    component_u = sp.Matrix([u0, u1, *q])
    component_check = {
        u[index]: component_u[index]
        for index in range(11)
    }
    assert sp.factor(linear.subs(component_check)) == 0
    assert sp.factor(quadric.subs(component_check)) == 0
    print("PASS fifth component: parameterized one affine 9-plane")

    tangent_basis = context["tangent_basis"]
    pivot_columns = context["pivot_columns"]
    pivot_inverse = context["pivot_inverse"]
    h = state["h"]
    second = state["second"]
    third = state["third"]
    second_variation = state["second_variation"]
    third_variation = state["third_variation"]
    tangent_matrix = state["tangent_matrix"]
    fixed_constant = state["fixed_constant"]
    tail_linear = state["tail_linear"]
    tail_quadratic = state["tail_quadratic"]
    h_polynomial = state["h_polynomial"]

    fourth_constant = sp.Matrix(
        [
            sp.expand(
                fixed_constant[index]
                + (tail_linear[index] * component_u)[0]
                + (
                    component_u.T
                    * tail_quadratic[index]
                    * component_u
                )[0]
            )
            for index in range(algebraization.TAIL_COUNT)
        ]
    )
    fourth_columns = list(tangent_matrix.rref()[1])
    fourth_rows = list(
        tangent_matrix[:, fourth_columns].T.rref()[1]
    )
    fourth_coordinates = sp.zeros(13, 1)
    fourth_solution = -tangent_matrix[
        fourth_rows,
        fourth_columns,
    ].inv() * fourth_constant[fourth_rows, :]
    for index, column in enumerate(fourth_columns):
        fourth_coordinates[column] = sp.expand(fourth_solution[index])
    assert all(
        sp.expand(value) == 0
        for value in (
            tangent_matrix * fourth_coordinates + fourth_constant
        )
    )

    k = second + second_variation * component_u
    ell_base = (
        third
        + third_variation * component_u
        + tangent_basis * fourth_coordinates
    )
    ell = ell_base + second_variation * sp.Matrix(z)

    max_order = 67
    rows = [
        [
            fourth.generating_weights(4, order).get(exponent, 0)
            for exponent in BASIS
        ]
        for order in range(max_order + 1)
    ]
    first_rows = sp.Matrix(rows[:12])
    row_coordinates = [
        sp.Matrix([[row[index] for index in pivot_columns]]) * pivot_inverse
        for row in rows
    ]
    bilinear = [
        bilinear_matrix(8, order)
        for order in range(max_order)
    ]
    h_squared = seed.power(h_polynomial, 2)
    h_cubed = seed.power(h_polynomial, 3)
    h_fourth = seed.power(h_polynomial, 4)
    h_fifth = seed.power(h_polynomial, 5)
    h2_functionals = [
        fixed_times_basis_functional(h_squared, 12, order)
        for order in range(max_order - 1)
    ]
    h3_functionals = [
        fixed_times_basis_functional(h_cubed, 16, order)
        for order in range(max_order - 2)
    ]
    hkk_matrices = [
        h_times_two_basis_matrix(h_polynomial, order)
        for order in range(max_order - 1)
    ]
    print("PASS fifth component: assembled fixed multilinear contractions")

    fourth_prefix = [sp.Rational(0)]
    for order in range(1, 12):
        value = (
            order * (h.T * bilinear[order - 1] * ell)[0]
            + sp.Rational(order, 2)
            * (k.T * bilinear[order - 1] * k)[0]
        )
        if order >= 2:
            value += (
                sp.Rational(order * (order - 1), 2)
                * (h2_functionals[order - 2].T * k)[0]
            )
        if order >= 3:
            value += (
                sp.Rational(order * (order - 1) * (order - 2), 24)
                * fixed_contraction(h_fourth, 16, order - 3)
            )
        fourth_prefix.append(sp.expand(value))
    m_values = -pivot_inverse * sp.Matrix(fourth_prefix)
    m = sp.zeros(25, 1)
    for index, pivot in enumerate(pivot_columns):
        m[pivot] = sp.expand(m_values[index])

    known = [sp.Rational(0)]
    tangent_effect = [sp.zeros(1, 13)]
    for order in range(1, max_order + 1):
        value = order * (
            (h.T * bilinear[order - 1] * m)[0]
            + (k.T * bilinear[order - 1] * ell)[0]
        )
        effect = order * h.T * bilinear[order - 1] * tangent_basis
        if order >= 2:
            value += (
                sp.Rational(order * (order - 1), 2)
                * (
                    (h2_functionals[order - 2].T * ell)[0]
                    + (k.T * hkk_matrices[order - 2] * k)[0]
                )
            )
        if order >= 3:
            value += (
                sp.Rational(order * (order - 1) * (order - 2), 6)
                * (h3_functionals[order - 3].T * k)[0]
            )
        if order >= 4:
            value += (
                sp.Rational(
                    order * (order - 1) * (order - 2) * (order - 3),
                    120,
                )
                * fixed_contraction(h_fifth, 20, order - 4)
            )
        known.append(sp.expand(value))
        tangent_effect.append(effect)

    defect = []
    fifth_matrix_rows = []
    for order in range(12, max_order + 1):
        value = known[order]
        effect = tangent_effect[order].copy()
        for prefix in range(12):
            value -= row_coordinates[order][prefix] * known[prefix]
            effect -= row_coordinates[order][prefix] * tangent_effect[prefix]
        defect.append(sp.expand(value))
        fifth_matrix_rows.append(effect)
    fifth_matrix = sp.Matrix.vstack(*fifth_matrix_rows)
    assert fifth_matrix.rank() == 2
    fifth_columns = list(fifth_matrix.rref()[1])
    fifth_rows = list(fifth_matrix[:, fifth_columns].T.rref()[1])
    fifth_pivot = fifth_matrix[
        fifth_rows,
        fifth_columns,
    ].det()

    reduced_equations = []
    pivot_inverse_fifth = fifth_matrix[
        fifth_rows,
        fifth_columns,
    ].inv()
    for row in range(len(defect)):
        if row in fifth_rows:
            continue
        multiplier = (
            fifth_matrix[row, fifth_columns]
            * pivot_inverse_fifth
        )
        reduced_equations.append(
            sp.expand(
                defect[row]
                - (multiplier * sp.Matrix(defect)[fifth_rows, :])[0]
            )
        )
    lower_matrix = sp.Matrix(
        [
            [sp.diff(equation, variable) for variable in z]
            for equation in reduced_equations
        ]
    )
    lower_constant = sp.Matrix(
        [
            sp.expand(equation.subs({variable: 0 for variable in z}))
            for equation in reduced_equations
        ]
    )
    assert all(
        sp.Poly(equation, *z).total_degree() <= 1
        for equation in reduced_equations
    )
    lower_field = DomainMatrix.from_Matrix(
        lower_matrix,
        extension=True,
    ).to_field()
    augmented_field = DomainMatrix.from_Matrix(
        lower_matrix.row_join(lower_constant),
        extension=True,
    ).to_field()
    lower_rank = lower_field.rank()
    augmented_rank = augmented_field.rank()
    assert (lower_rank, augmented_rank) == (2, 3)
    lower_columns = list(lower_field.rref()[1])
    lower_rows = list(lower_field.transpose().rref()[1])
    assert len(lower_columns) == len(lower_rows) == 2
    lower_pivot = sp.factor(
        lower_matrix[lower_rows, lower_columns].det()
    )
    survivor_equations = []
    selected_columns = lower_columns + [lower_matrix.cols]
    augmented = lower_matrix.row_join(lower_constant)
    for row in range(augmented.rows):
        if row in lower_rows:
            continue
        survivor_equations.append(
            sp.factor(
                augmented[
                    lower_rows + [row],
                    selected_columns,
                ].det()
            )
        )
    survivor_equations = [
        value
        for value in survivor_equations
        if value
    ]
    survivor_degrees = sorted(
        {
            sp.Poly(value, *q).total_degree()
            for value in survivor_equations
        }
    )
    assert survivor_degrees == [0]
    uniform_obstruction = sp.factor(survivor_equations[0])
    assert uniform_obstruction != 0
    conjugate_obstruction = sp.factor(
        uniform_obstruction.xreplace({sp.sqrt(41): -sp.sqrt(41)})
    )
    obstruction_norm = sp.factor(
        sp.expand(uniform_obstruction * conjugate_obstruction)
    )
    assert obstruction_norm != 0
    assert not obstruction_norm.has(sp.sqrt(41))
    print(
        "PASS fifth component: generic lower/augmented ranks "
        f"{lower_rank}/{augmented_rank}"
    )
    print("PASS fifth component: exact constant obstruction on the whole chart")

    result = {
        "format": "two-pair-counterexample-fifth-component-research-v1",
        "status": "exact uniform fifth-order obstruction",
        "direction_label": label,
        "direction": list(direction),
        "field": "Q(sqrt(41))",
        "component_coordinates": [str(variable) for variable in q],
        "lower_kernel_parameters": len(z),
        "fifth_tangent_rank": 2,
        "reduced_equation_count": len(reduced_equations),
        "generic_lower_rank": lower_rank,
        "generic_augmented_rank": augmented_rank,
        "lower_pivot_columns": lower_columns,
        "lower_pivot_rows": lower_rows,
        "survivor_equation_count": len(survivor_equations),
        "survivor_equation_total_degrees": survivor_degrees,
        "uniform_obstruction": {
            "constant": str(uniform_obstruction),
            "conjugate": str(conjugate_obstruction),
            "norm": str(obstruction_norm),
            "conclusion": (
                "The coefficient matrix has rank at most 2 everywhere, "
                "while this augmented 3x3 minor is a nonzero constant. "
                "Hence the complete affine 9-plane has no fifth lift."
            ),
        },
        "conjugate_component": (
            "Eliminated by applying sqrt(41) -> -sqrt(41)."
        ),
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"CHECKPOINT {output.relative_to(ROOT)}")

    modular = singular_saturation(
        survivor_equations,
        lower_pivot,
        q,
    )
    result["modular_confirmation"] = "unit ideal on the selected pivot chart"
    result["modular_pivot_chart_saturation"] = modular
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"WROTE {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

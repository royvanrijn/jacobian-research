#!/usr/bin/env python3
"""Exploratory exact algebraization tests for the local SIC2C4 five-plane.

This is deliberately a research script rather than a status checker.  It
restores the complete eleven-parameter second-correction freedom at selected
directions, computes the fourth-order compatibility fiber over Q, and tests
one point on each conjugate component after restoring the eleven omitted
cubic-tangent parameters at fifth order.  Directions obstructed at fifth
order cannot be continued to orders 6--12.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


ROOT = Path(__file__).resolve().parents[1]
FOURTH_SCRIPT = ROOT / "scripts" / "verify_two_pair_counterexample_fourth_order.py"
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_algebraization_research.json"
)

spec = importlib.util.spec_from_file_location("sic2_fourth_research", FOURTH_SCRIPT)
assert spec and spec.loader
fourth = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fourth)

seed = fourth.seed
BASIS = fourth.BASIS
MAX_ORDER = fourth.MAX_ORDER
TAIL_COUNT = fourth.TAIL_ORDER_COUNT


def evaluate_vector_map(
    polynomial: fourth.VectorMap,
    values: tuple[int, int, int, int, int],
) -> sp.Matrix:
    result = sp.zeros(25, 1)
    for exponent, vector in polynomial.items():
        result += vector * sp.prod(
            values[index] ** exponent[index] for index in range(5)
        )
    return result


def vector_polynomial(vector: sp.Matrix) -> seed.Polynomial:
    return {
        BASIS[index]: sp.factor(value)
        for index, value in enumerate(vector)
        if value
    }


def contraction(
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


def singular_groebner(equations: list[sp.Expr]) -> dict[str, object]:
    singular = shutil.which("Singular")
    assert singular, "Singular is required for exact fourth-order fibers"
    nonzero = [sp.factor(equation) for equation in equations if equation]
    source = "\n".join(
        [
            "ring r=0,(u0,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10),dp;",
            f"ideal I={','.join(str(eq).replace('**', '^') for eq in nonzero)};",
            "ideal G=std(I);",
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
        timeout=300,
    )
    lines = [line.strip() for line in completed.stdout.splitlines()]

    def integer_after(label: str) -> int:
        return int(lines[lines.index(label) + 1])

    basis = lines[lines.index("GB") + 1 : lines.index("END")]
    return {
        "dimension": integer_after("DIM"),
        "degree": integer_after("DEG"),
        "size": integer_after("SIZE"),
        "groebner_basis": basis,
    }


def selected_component_points(
    groebner_basis: list[str],
) -> tuple[list[sp.Matrix], int]:
    variables = sp.symbols("u0:11")
    table = {str(variable): variable for variable in variables}
    polynomials = [
        sp.sympify(
            value.rstrip(",").replace("^", "**"),
            locals=table,
        )
        for value in groebner_basis
    ]
    linear = next(
        value
        for value in polynomials
        if sp.Poly(value, *variables).total_degree() == 1
    )
    quadric = next(
        value
        for value in polynomials
        if sp.Poly(value, *variables).total_degree() == 2
    )
    specialization = {variable: 0 for variable in variables[2:]}
    univariate = sp.Poly(quadric.subs(specialization), variables[1])
    discriminant = int(sp.discriminant(univariate.as_expr()))
    square_class = 1
    for prime, exponent in sp.factorint(discriminant).items():
        if exponent % 2:
            square_class *= prime
    points = []
    for root in sp.solve(univariate.as_expr(), variables[1]):
        point = sp.zeros(11, 1)
        point[1] = sp.factor(root)
        point[0] = sp.factor(
            sp.solve(
                linear.subs(
                    {
                        variables[index]: point[index]
                        for index in range(1, 11)
                    }
                ),
                variables[0],
            )[0]
        )
        substitution = {
            variables[index]: point[index]
            for index in range(11)
        }
        assert all(sp.factor(value.subs(substitution)) == 0 for value in polynomials)
        points.append(point)
    assert len(points) == 2
    return points, square_class


def fifth_test(
    context: dict[str, object],
    h: sp.Matrix,
    k: sp.Matrix,
    ell: sp.Matrix,
    m: sp.Matrix,
    lower_variations: list[tuple[sp.Matrix, sp.Matrix]] | None = None,
) -> dict[str, object]:
    tangent_basis = context["tangent_basis"]
    max_order = 67
    rows = [
        [
            fourth.generating_weights(4, order).get(exponent, 0)
            for exponent in BASIS
        ]
        for order in range(max_order + 1)
    ]
    first_rows = sp.Matrix(rows[:12])
    pivot_columns = list(first_rows.rref()[1])
    pivot_inverse = first_rows[:, pivot_columns].inv()
    row_coordinates = [
        sp.Matrix([[row[index] for index in pivot_columns]]) * pivot_inverse
        for row in rows
    ]
    h_poly, k_poly, ell_poly, m_poly = (
        vector_polynomial(vector) for vector in (h, k, ell, m)
    )
    hm = seed.multiply(h_poly, m_poly)
    kell = seed.multiply(k_poly, ell_poly)
    h2ell = seed.multiply(seed.power(h_poly, 2), ell_poly)
    hk2 = seed.multiply(h_poly, seed.power(k_poly, 2))
    h3k = seed.multiply(seed.power(h_poly, 3), k_poly)
    h5 = seed.power(h_poly, 5)
    known = []
    tangent_effect = []
    for order in range(max_order + 1):
        value = sp.Rational(0)
        effect = sp.zeros(1, 13)
        if order >= 1:
            value += order * (
                contraction(hm, 8, order - 1)
                + contraction(kell, 8, order - 1)
            )
            weights = fourth.generating_weights(8, order - 1)
            bilinear = sp.zeros(25)
            for left, left_exponent in enumerate(BASIS):
                for right, right_exponent in enumerate(BASIS):
                    exponent = tuple(
                        left_exponent[index] + right_exponent[index]
                        for index in range(4)
                    )
                    bilinear[left, right] = weights.get(exponent, 0)
            effect = order * h.T * bilinear * tangent_basis
        if order >= 2:
            value += sp.Rational(order * (order - 1), 2) * (
                contraction(h2ell, 12, order - 2)
                + contraction(hk2, 12, order - 2)
            )
        if order >= 3:
            value += (
                sp.Rational(order * (order - 1) * (order - 2), 6)
                * contraction(h3k, 16, order - 3)
            )
        if order >= 4:
            value += (
                sp.Rational(
                    order * (order - 1) * (order - 2) * (order - 3),
                    120,
                )
                * contraction(h5, 20, order - 4)
            )
        known.append(sp.factor(value))
        tangent_effect.append(effect)
    tail_constant = []
    tail_rows = []
    for order in range(12, max_order + 1):
        constant = known[order]
        effect = tangent_effect[order].copy()
        for prefix in range(12):
            constant -= row_coordinates[order][prefix] * known[prefix]
            effect -= row_coordinates[order][prefix] * tangent_effect[prefix]
        tail_constant.append(sp.factor(constant))
        tail_rows.append([sp.factor(value) for value in effect])
    matrix = sp.Matrix(tail_rows)
    constant = sp.Matrix(tail_constant)
    frozen_domain = DomainMatrix.from_Matrix(matrix, extension=True)
    frozen_augmented_domain = DomainMatrix.from_Matrix(
        matrix.row_join(-constant),
        extension=True,
    )
    frozen_rank = frozen_domain.rank()
    frozen_augmented_rank = frozen_augmented_domain.rank()
    variation_columns = []
    for ell_delta, m_delta in lower_variations or []:
        ell_delta_poly = vector_polynomial(ell_delta)
        m_delta_poly = vector_polynomial(m_delta)
        hm_delta = seed.multiply(h_poly, m_delta_poly)
        kell_delta = seed.multiply(k_poly, ell_delta_poly)
        h2ell_delta = seed.multiply(
            seed.power(h_poly, 2),
            ell_delta_poly,
        )
        raw_delta = []
        for order in range(max_order + 1):
            value = sp.Rational(0)
            if order >= 1:
                value += order * (
                    contraction(hm_delta, 8, order - 1)
                    + contraction(kell_delta, 8, order - 1)
                )
            if order >= 2:
                value += (
                    sp.Rational(order * (order - 1), 2)
                    * contraction(h2ell_delta, 12, order - 2)
                )
            raw_delta.append(sp.factor(value))
        tail_delta = []
        for order in range(12, max_order + 1):
            value = raw_delta[order]
            for prefix in range(12):
                value -= (
                    row_coordinates[order][prefix]
                    * raw_delta[prefix]
                )
            tail_delta.append(sp.factor(value))
        variation_columns.append(sp.Matrix(tail_delta))
    if variation_columns:
        matrix = sp.Matrix.hstack(*variation_columns, matrix)
    domain_matrix = DomainMatrix.from_Matrix(matrix, extension=True)
    augmented_domain = DomainMatrix.from_Matrix(
        matrix.row_join(-constant),
        extension=True,
    )
    rank = domain_matrix.rank()
    augmented_rank = augmented_domain.rank()
    report: dict[str, object] = {
        "frozen_coefficient_rank": frozen_rank,
        "frozen_augmented_rank": frozen_augmented_rank,
        "coefficient_rank": rank,
        "augmented_rank": augmented_rank,
        "consistent": rank == augmented_rank,
        "restored_lower_kernel_parameters": len(variation_columns),
    }
    if rank == augmented_rank:
        reduced, pivots = augmented_domain.rref()
        reduced_matrix = reduced.to_Matrix()
        point = [sp.Rational(0) for _ in range(matrix.cols)]
        for row, pivot in enumerate(pivots):
            if pivot < matrix.cols:
                point[pivot] = sp.factor(reduced_matrix[row, matrix.cols])
        report["one_solution"] = [str(value) for value in point]
        report["solution_space_dimension"] = matrix.cols - rank
        return report
    return report


def fourth_fiber(
    context: dict[str, object],
    values: tuple[int, int, int, int, int],
    *,
    run_point_tests: bool = True,
    state_target: dict[str, object] | None = None,
) -> dict[str, object]:
    tangent_basis = context["tangent_basis"]
    residual_vectors = context["residual_vectors"]
    pivot_columns = context["pivot_columns"]
    pivot_inverse = context["pivot_inverse"]
    row_coordinates = context["row_coordinates"]
    bilinear_matrices = context["bilinear_matrices"]
    second_map = context["second_map"]
    third_map = context["third_map"]
    symbolic_tangent_matrix = context["symbolic_tangent_matrix"]
    symbolic_fourth_constant = context["symbolic_fourth_constant"]

    h_symbols = sp.symbols("h0:5")
    substitution = dict(zip(h_symbols, values))
    h = residual_vectors * sp.Matrix(values)
    second = evaluate_vector_map(second_map, values)
    tangent_matrix = symbolic_tangent_matrix.subs(substitution)
    tangent_rank = tangent_matrix.rank()
    tangent_kernel = sp.Matrix.hstack(*tangent_matrix.nullspace())
    second_variation = tangent_basis * tangent_kernel

    third_variation_base = sp.zeros(12, tangent_kernel.cols)
    for order in range(1, 12):
        third_variation_base[order, :] = (
            order
            * h.T
            * bilinear_matrices[order - 1]
            * second_variation
        )
    third_variation = sp.zeros(25, tangent_kernel.cols)
    third_pivot_variation = -pivot_inverse * third_variation_base
    for index, pivot in enumerate(pivot_columns):
        third_variation[pivot, :] = third_pivot_variation[index, :]

    h_polynomial = vector_polynomial(h)
    h_squared = seed.multiply(h_polynomial, h_polynomial)
    h2_functionals = []
    for order in range(MAX_ORDER - 1):
        weights = fourth.generating_weights(12, order)
        functional = []
        for basis_exponent in BASIS:
            numerator = seed.multiply(
                h_squared,
                seed.monomial(basis_exponent),
            )
            functional.append(
                sum(
                    coefficient * weights.get(exponent, 0)
                    for exponent, coefficient in numerator.items()
                )
            )
        h2_functionals.append(sp.Matrix(functional))

    raw_linear = [sp.zeros(1, tangent_kernel.cols)]
    raw_quadratic = [sp.zeros(tangent_kernel.cols)]
    for order in range(1, MAX_ORDER + 1):
        bilinear = bilinear_matrices[order - 1]
        linear = (
            order * h.T * bilinear * third_variation
            + order * second.T * bilinear * second_variation
        )
        if order >= 2:
            linear += (
                sp.Rational(order * (order - 1), 2)
                * h2_functionals[order - 2].T
                * second_variation
            )
        raw_linear.append(linear)
        raw_quadratic.append(
            sp.Rational(order, 2)
            * second_variation.T
            * bilinear
            * second_variation
        )

    tail_linear = []
    tail_quadratic = []
    for order in range(12, MAX_ORDER + 1):
        linear = raw_linear[order].copy()
        quadratic = raw_quadratic[order].copy()
        for prefix in range(12):
            linear -= row_coordinates[order][prefix] * raw_linear[prefix]
            quadratic -= (
                row_coordinates[order][prefix] * raw_quadratic[prefix]
            )
        tail_linear.append(linear)
        tail_quadratic.append(quadratic)

    left_kernel = sp.Matrix.hstack(*tangent_matrix.T.nullspace()).T
    fixed_constant = symbolic_fourth_constant.subs(substitution)
    free_variables = sp.symbols(f"u0:{tangent_kernel.cols}")
    free_vector = sp.Matrix(free_variables)
    equations = []
    for row in range(left_kernel.rows):
        equation = sp.Integer(0)
        for tail_index in range(TAIL_COUNT):
            equation += left_kernel[row, tail_index] * (
                fixed_constant[tail_index]
                + (tail_linear[tail_index] * free_vector)[0]
                + (
                    free_vector.T
                    * tail_quadratic[tail_index]
                    * free_vector
                )[0]
            )
        equations.append(sp.factor(equation))

    report = {
        "direction": list(values),
        "fourth_tail_rank": tangent_rank,
        "free_second_correction_parameters": tangent_kernel.cols,
        "compatibility_equation_count": sum(bool(eq) for eq in equations),
    }
    if tangent_kernel.cols != 11:
        report["status"] = "exceptional tangent-rank chart; not decomposed"
        return report
    groebner = singular_groebner(equations)
    report.update(groebner)
    report["status"] = (
        "fourth lifts exist"
        if groebner["dimension"] >= 0
        else "obstructed at fourth order"
    )
    if state_target is not None:
        state_target.update(
            {
                "h": h,
                "second": second,
                "third": evaluate_vector_map(third_map, values),
                "tangent_matrix": tangent_matrix,
                "tangent_kernel": tangent_kernel,
                "second_variation": second_variation,
                "third_variation": third_variation,
                "h_polynomial": h_polynomial,
                "h2_functionals": h2_functionals,
                "fixed_constant": fixed_constant,
                "tail_linear": tail_linear,
                "tail_quadratic": tail_quadratic,
                "groebner_basis": groebner["groebner_basis"],
            }
        )
    if not run_point_tests:
        return report
    if groebner["dimension"] == 9 and groebner["degree"] == 2:
        points, square_class = selected_component_points(
            groebner["groebner_basis"]
        )
        report["discriminant_square_class"] = square_class
        selected_reports = []
        third = evaluate_vector_map(third_map, values)
        h_cubed = seed.power(h_polynomial, 3)
        h_fourth = seed.power(h_polynomial, 4)
        for point in points[:1]:
            explicit_constant = sp.Matrix(
                [
                    fixed_constant[index]
                    + (tail_linear[index] * point)[0]
                    + (
                        point.T
                        * tail_quadratic[index]
                        * point
                    )[0]
                    for index in range(TAIL_COUNT)
                ]
            )
            tangent_columns = list(tangent_matrix.rref()[1])
            tangent_rows = list(
                tangent_matrix[:, tangent_columns].T.rref()[1]
            )
            tangent_coordinates = sp.zeros(13, 1)
            solution = -tangent_matrix[
                tangent_rows,
                tangent_columns,
            ].inv() * explicit_constant[tangent_rows, :]
            for index, column in enumerate(tangent_columns):
                tangent_coordinates[column] = sp.factor(solution[index])
            explicit_second = sp.simplify(second + second_variation * point)
            explicit_third = sp.simplify(
                third
                + third_variation * point
                + tangent_basis * tangent_coordinates
            )
            fourth_prefix = []
            for order in range(12):
                if order == 0:
                    fourth_prefix.append(sp.Rational(0))
                    continue
                value = (
                    order
                    * (
                        h.T
                        * bilinear_matrices[order - 1]
                        * explicit_third
                    )[0]
                    + sp.Rational(order, 2)
                    * (
                        explicit_second.T
                        * bilinear_matrices[order - 1]
                        * explicit_second
                    )[0]
                )
                if order >= 2:
                    value += (
                        sp.Rational(order * (order - 1), 2)
                        * (
                            h2_functionals[order - 2].T
                            * explicit_second
                        )[0]
                    )
                if order >= 3:
                    value += (
                        sp.Rational(
                            order * (order - 1) * (order - 2),
                            24,
                        )
                        * contraction(h_fourth, 16, order - 3)
                    )
                fourth_prefix.append(sp.factor(value))
            pivot_values = -pivot_inverse * sp.Matrix(fourth_prefix)
            explicit_fourth = sp.zeros(25, 1)
            for index, pivot in enumerate(pivot_columns):
                explicit_fourth[pivot] = sp.factor(pivot_values[index])
            selected_reports.append(
                {
                    "free_second_coordinates": [
                        str(sp.factor(value)) for value in point
                    ],
                    "conjugate_component_checked_by_field_automorphism": True,
                }
            )
            lower_variations = []
            for kernel_index in range(tangent_kernel.cols):
                ell_delta = second_variation[:, kernel_index]
                fourth_delta_prefix = [sp.Rational(0)]
                for order in range(1, 12):
                    fourth_delta_prefix.append(
                        sp.factor(
                            order
                            * (
                                h.T
                                * bilinear_matrices[order - 1]
                                * ell_delta
                            )[0]
                        )
                    )
                m_delta_values = -pivot_inverse * sp.Matrix(
                    fourth_delta_prefix
                )
                m_delta = sp.zeros(25, 1)
                for index, pivot in enumerate(pivot_columns):
                    m_delta[pivot] = sp.factor(m_delta_values[index])
                lower_variations.append((ell_delta, m_delta))
            selected_reports[-1]["fifth_order_complete_at_point"] = fifth_test(
                context,
                h,
                explicit_second,
                explicit_third,
                explicit_fourth,
                lower_variations=lower_variations,
            )
        report["selected_component_points"] = selected_reports
    return report


def main() -> None:
    context = fourth.main(write_output=False, return_context=True)
    assert context is not None
    directions = {
        "generic_rational": (2, -1, 3, 1, -2),
        "quadratic_extension_41": (1, 2, 3, 4, 5),
        "pure_apolar_odd": (0, 1, 0, 0, 0),
        "apolar_even_control": (-210, 0, 18, 69, 64),
    }
    fibers = {
        label: fourth_fiber(context, values)
        for label, values in directions.items()
    }
    result = {
        "format": "two-pair-counterexample-algebraization-research-v2",
        "scope": (
            "Exact fourth-order fibers and complete pointwise fifth tests "
            "for selected reduced five-plane directions; exploratory, not "
            "a theorem-status artifact."
        ),
        "known_exact_family_control": {
            "parameterization": "F_(1+s,1)",
            "maximum_parameter_degree": 3,
            "coefficient_formulas": {
                "A0": "R+Z",
                "B0": "2*W*(R+Z)^2-2*R^3-R^2*Z",
                "B1": "4*W*(R+Z)*R-2*R^3",
                "B2": "2*W*R^2",
                "order_0": "A0*B0/2",
                "order_1": "(R*B0+A0*B1)/2",
                "order_2": "(R*B1+A0*B2)/2",
                "order_3": "R*B2/2=W*R^3",
                "orders_4_through_12": "0",
            },
            "all_higher_coefficients": "zero from order 4 onward",
            "all_order_identity": "E_2(F_(1+s,1)^m)=0 for every m>=1",
            "reconstruction": (
                "exact polynomial (hence rational and algebraic) "
                "reconstruction of degree 3"
            ),
            "coefficient_recurrence": "C_n=0 for every n>=4",
        },
        "directions": fibers,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    for label, report in fibers.items():
        print(
            f"{label}: fourth fiber dim/deg "
            f"{report.get('dimension')}/{report.get('degree')}"
        )
    print(f"WROTE {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

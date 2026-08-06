#!/usr/bin/env python3
"""Exact nonlinear lift and order-seven test for the degree-eight survivor.

The characteristic-zero component is read from the independently reconstructed
pivot-chart artifact.  This verifier recomputes the relative ranks over its
residue field, recognizes the reduced affine-five-space inside the nonlinear
order-five lift, retains all four free order-five correction directions, and
tests the complete restricted order-seven obstruction ideal.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import json

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_quantum_residue import (
    add,
    column_rank,
    pi_power,
    poisson,
    scale,
    split_pair,
)
from interpolate_degree_eight_order_five_chart import PIVOT_COLUMNS
from screen_degree_seven_order_five_survivors import (
    projected_order_five_equations,
)
from verify_degree_eight_relative_quantization_obstruction import (
    S4_SUPPORT,
    T4_SUPPORT,
    family_presentation,
)
from verify_degree_seven_order_five_survivor import (
    binomial_term,
    coefficient_polynomial_to_field,
    coefficient_text,
    fifth_defect,
    graded_indices,
    newton_interpolation,
    newton_value,
    run_singular,
    singular_polynomial,
)
from verify_degree_seven_relative_quantization_obstruction import weight_monomials


def deserialize_lex(data):
    tau, sigma = sp.symbols("tau sigma")
    polynomials = []
    for terms in data["exact_zero_scheme"]["lexicographic_basis"]:
        expression = sum(
            sp.Rational(term["coefficient"])
            * tau ** term["tau_degree"]
            * sigma ** term["sigma_degree"]
            for term in terms
        )
        polynomials.append(sp.Poly(expression, tau, sigma, domain=QQ))
    return tau, sigma, polynomials


def component_context(path: Path):
    data = json.loads(path.read_text())
    tau_symbol, sigma_symbol, lex = deserialize_lex(data)
    minimum = sp.Poly(
        sp.sympify(
            data["exact_zero_scheme"]["primitive_sigma_polynomial"],
            locals={"sigma": sigma_symbol},
        ),
        sigma_symbol,
        domain=QQ,
    )
    field = QQ.alg_field_from_poly(minimum, alias="alpha")
    sigma = field.convert(field.ext)
    tau_relation = next(polynomial for polynomial in lex if polynomial.degree(tau_symbol) == 1)
    tau_coefficient = sp.Poly(
        tau_relation.as_expr().coeff(tau_symbol, 1), sigma_symbol, domain=QQ
    )
    constant = sp.Poly(
        tau_relation.as_expr().coeff(tau_symbol, 0), sigma_symbol, domain=QQ
    )
    tau = -coefficient_polynomial_to_field(
        constant.as_expr(), sigma_symbol, field, sigma
    ) / coefficient_polynomial_to_field(
        tau_coefficient.as_expr(), sigma_symbol, field, sigma
    )
    return data, minimum, field, sigma, tau


def minimum_text(minimum, alias: str):
    return str(minimum.as_expr()).replace(
        str(minimum.gens[0]), alias
    ).replace("**", "^")


def verify_reconstruction_chart(presentation, field):
    strong_columns = presentation["strong_columns"]
    monomials = sorted(
        set(presentation["constant"]).union(
            *(set(column) for column in strong_columns)
        )
    )
    pivot_monomials = monomials[:104]
    rows = {
        row_index: {
            column_index: coefficient
            for column_index, strong_index in enumerate(PIVOT_COLUMNS)
            if (
                coefficient := strong_columns[strong_index].get(
                    monomial, field.zero
                )
            )
        }
        for row_index, monomial in enumerate(pivot_monomials)
    }
    _, pivots, _ = sdm_irref(rows)
    assert pivots == list(range(104))
    return {
        "pivot_columns": len(PIVOT_COLUMNS),
        "pivot_rows": len(pivot_monomials),
        "rank_over_residue_field": len(pivots),
        "component_inside_chart": True,
    }


def parse_basis_polynomial(text, variable_count, field, sigma):
    a_symbol = sp.symbols("a")
    variables = sp.symbols(f"z0:{variable_count}")
    expression = sp.sympify(
        text.replace("^", "**"),
        locals={"a": a_symbol, **{str(variable): variable for variable in variables}},
    )
    terms = {}
    for exponent, coefficient in sp.Poly(expression, *variables).terms():
        terms[exponent] = coefficient_polynomial_to_field(
            coefficient, a_symbol, field, sigma
        )
    return terms


def polynomial_total_degree(terms):
    return max((sum(exponent) for exponent in terms), default=-1)


def active_indices(terms):
    return {
        index
        for exponent in terms
        for index, degree in enumerate(exponent)
        if degree
    }


def univariate_add(left, right, field):
    result = dict(left)
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, field.zero) + coefficient
        if not result[degree]:
            del result[degree]
    return result


def univariate_scale(polynomial, coefficient):
    return {
        degree: coefficient * value
        for degree, value in polynomial.items()
        if coefficient * value
    }


def univariate_multiply(left, right, field):
    result = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            result[degree] = result.get(degree, field.zero) + left_coefficient * right_coefficient
            if not result[degree]:
                del result[degree]
    return result


def univariate_power(polynomial, exponent, field):
    result = {0: field.one}
    for _ in range(exponent):
        result = univariate_multiply(result, polynomial, field)
    return result


def substitute_reduced_roots(terms, root7, root9, field):
    result = {}
    for exponent, coefficient in terms.items():
        assert all(not exponent[index] for index in range(7))
        term = univariate_power(root7, exponent[7], field)
        term = univariate_multiply(term, {exponent[8]: field.one}, field)
        term = univariate_scale(term, coefficient * root9 ** exponent[9])
        result = univariate_add(result, term, field)
    return result


def reduced_lift_structure(basis_text, field, sigma):
    polynomials = [
        parse_basis_polynomial(text, 10, field, sigma) for text in basis_text
    ]
    linear = [terms for terms in polynomials if polynomial_total_degree(terms) == 1]
    quadratic = [terms for terms in polynomials if polynomial_total_degree(terms) == 2]
    assert len(linear) == len(quadratic) == 3
    q9 = next(terms for terms in quadratic if active_indices(terms) == {9})
    q78 = next(
        terms
        for terms in quadratic
        if active_indices(terms).issubset({7, 8})
        and any(exponent[7] == 2 for exponent in terms)
    )
    bridge = next(terms for terms in quadratic if terms is not q9 and terms is not q78)

    zero = (0,) * 10
    e9 = tuple(1 if index == 9 else 0 for index in range(10))
    e99 = tuple(2 if index == 9 else 0 for index in range(10))
    leading9 = q9[e99]
    linear9 = q9.get(e9, field.zero)
    constant9 = q9.get(zero, field.zero)
    assert linear9**2 - field(4) * leading9 * constant9 == field.zero
    root9 = -linear9 / (field(2) * leading9)

    leading7 = next(
        coefficient
        for exponent, coefficient in q78.items()
        if exponent[7] == 2
    )
    linear7 = {
        exponent[8]: coefficient
        for exponent, coefficient in q78.items()
        if exponent[7] == 1
    }
    root7 = univariate_scale(linear7, -field.one / (field(2) * leading7))
    assert substitute_reduced_roots(q9, root7, root9, field) == {}
    assert substitute_reduced_roots(q78, root7, root9, field) == {}
    assert substitute_reduced_roots(bridge, root7, root9, field) == {}

    relations = {}
    for dependent in (5, 3, 1):
        unit = tuple(1 if index == dependent else 0 for index in range(10))
        relation = next((terms for terms in linear if unit in terms), None)
        assert relation is not None
        relations[dependent] = relation
    return {
        "relations": relations,
        "root7": root7,
        "root9": root9,
        "free_coordinates": (0, 2, 4, 6, 8),
    }


def exact_order_five(presentation, field, minimum, sigma):
    equations, correction_rank, dual_rank = projected_order_five_equations(
        presentation, field
    )
    equation_monomials = sorted(set().union(*(set(equation) for equation in equations)))
    equation_index = {
        monomial: index for index, monomial in enumerate(equation_monomials)
    }
    equation_rows = {
        equation_index[monomial]: {
            column_index: coefficient
            for column_index, equation in enumerate(equations)
            if (coefficient := equation.get(monomial, field.zero))
        }
        for monomial in equation_monomials
    }
    _, equation_pivots, _ = sdm_irref(equation_rows)
    independent_equations = [equations[index] for index in equation_pivots]
    assert len(independent_equations) == 7
    variables = tuple(f"z{index}" for index in range(10))
    ideal = ",".join(
        singular_polynomial(field, equation, variables, "a")
        for equation in independent_equations
    )
    program = f"""ring r=(0,a),({','.join(variables)}),dp;
minpoly={minimum_text(minimum, 'a')};
option(redSB);
ideal I={ideal};
ideal G=std(I);
print("SIZE="+string(size(G)));
print("DIMENSION="+string(dim(G)));
print("VDIM="+string(vdim(G)));
print("REDUCED_ONE="+string(reduce(1,G)));
for (int i=1; i<=size(G); i++) {{ print("G_"+string(i)+"="+string(G[i])); }}
"""
    output = run_singular(program)
    values = {}
    basis = []
    for line in output.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.startswith("G_"):
            basis.append(value)
        elif key in {"SIZE", "DIMENSION", "VDIM", "REDUCED_ONE"}:
            values[key] = value
    assert values == {
        "SIZE": "6",
        "DIMENSION": "5",
        "VDIM": "-1",
        "REDUCED_ONE": "1",
    }
    assert correction_rank == 97 and dual_rank == 172
    return (
        equations,
        basis,
        reduced_lift_structure(basis, field, sigma),
        correction_rank,
        dual_rank,
        len(independent_equations),
    )


def lower_lift_at(point, presentation, field, structure):
    parameters = [field.zero] * 10
    for index, value in zip(structure["free_coordinates"], point, strict=True):
        parameters[index] = field(value)
    parameters[9] = structure["root9"]
    parameters[7] = sum(
        (
            coefficient * parameters[8] ** degree
            for degree, coefficient in structure["root7"].items()
        ),
        field.zero,
    )
    zero = (0,) * 10
    for dependent in (5, 3, 1):
        terms = structure["relations"][dependent]
        unit = tuple(1 if index == dependent else 0 for index in range(10))
        value = terms.get(zero, field.zero)
        for exponent, coefficient in terms.items():
            if exponent == zero or exponent == unit:
                continue
            assert sum(exponent) == 1
            value += coefficient * parameters[exponent.index(1)]
        parameters[dependent] = -value / terms[unit]

    s2, t2 = presentation["base_pair"]
    for coefficient, (basis_s, basis_t) in zip(
        parameters, presentation["kernel_pairs"], strict=True
    ):
        s2 = add(s2, basis_s, coefficient)
        t2 = add(t2, basis_t, coefficient)
    return s2, t2


def solve_affine_many(columns, right_hand_sides, field):
    """Solve one fixed correction map against many exact right-hand sides."""

    monomials = sorted(
        set().union(
            *(set(column) for column in columns),
            *(set(rhs) for rhs in right_hand_sides),
        )
    )
    monomial_index = {monomial: index for index, monomial in enumerate(monomials)}
    homogeneous_rows = {}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            homogeneous_rows.setdefault(monomial_index[monomial], {})[
                column_index
            ] = coefficient
    homogeneous_reduced, homogeneous_pivots, homogeneous_nonzero = sdm_irref(
        homogeneous_rows
    )
    kernel, _ = sdm_nullspace_from_rref(
        homogeneous_reduced,
        field.one,
        len(columns),
        homogeneous_pivots,
        homogeneous_nonzero,
    )

    rows = {row: dict(entries) for row, entries in homogeneous_rows.items()}
    for rhs_index, rhs in enumerate(right_hand_sides):
        augmented_column = len(columns) + rhs_index
        for monomial, coefficient in rhs.items():
            rows.setdefault(monomial_index[monomial], {})[
                augmented_column
            ] = -coefficient
    reduced, pivots, _ = sdm_irref(rows)
    assert pivots == homogeneous_pivots
    solutions = []
    for rhs_index in range(len(right_hand_sides)):
        augmented_column = len(columns) + rhs_index
        solution = {}
        for reduced_row, pivot in enumerate(pivots):
            value = reduced.get(reduced_row, {}).get(
                augmented_column, field.zero
            )
            if value:
                solution[pivot] = -value
        solutions.append(solution)
    return solutions, kernel, len(pivots)


def fifth_correction_interpolation(presentation, field, structure):
    indices = graded_indices(5, 2)
    width = len(presentation["correction_five"])
    holdout = (3, 1, 2, 4, 2)
    points = [*indices, holdout]
    right_hand_sides = []
    for point in points:
        s2, t2 = lower_lift_at(point, presentation, field, structure)
        right_hand_sides.append(
            scale(fifth_defect(presentation, s2, t2, field), -field.one)
        )
    solutions, kernel, rank = solve_affine_many(
        presentation["correction_five"], right_hand_sides, field
    )
    assert rank == 97 and len(kernel) == 4
    values = {}
    for point, vector in zip(indices, solutions[:-1], strict=True):
        values[point] = [vector.get(index, field.zero) for index in range(width)]
    coefficients = newton_interpolation(values, indices, field)
    actual = solutions[-1]
    assert newton_value(coefficients, holdout, field) == [
        actual.get(index, field.zero) for index in range(width)
    ]
    return coefficients, kernel, holdout


def order_seven_defect_at(
    point, presentation, field, structure, fifth_coefficients, fifth_kernel
):
    lower_point = point[:5]
    s2, t2 = lower_lift_at(lower_point, presentation, field, structure)
    correction = newton_value(fifth_coefficients, lower_point, field)
    for parameter, kernel_vector in zip(point[5:], fifth_kernel, strict=True):
        correction = [
            value + field(parameter) * kernel_vector.get(index, field.zero)
            for index, value in enumerate(correction)
        ]
    correction_vector = {
        index: value for index, value in enumerate(correction) if value
    }
    s4, t4 = split_pair(correction_vector, S4_SUPPORT, T4_SUPPORT)
    S, T = presentation["S"], presentation["T"]
    result = add(poisson(s2, t4), poisson(s4, t2))
    for left, right in ((s4, T), (S, t4), (s2, t2)):
        result = add(result, pi_power(left, right, 3), field.one / field(24))
    for left, right in ((s2, T), (S, t2)):
        result = add(result, pi_power(left, right, 5), field.one / field(1920))
    return add(result, pi_power(S, T, 7), field.one / field(322560))


def exact_order_seven(presentation, field, structure, minimum):
    fifth_coefficients, fifth_kernel, fifth_holdout = fifth_correction_interpolation(
        presentation, field, structure
    )
    S, T = presentation["S"], presentation["T"]
    s6_support = weight_monomials(38, 2, 16)
    t6_support = weight_monomials(34, 1, 17)
    assert len(s6_support) == 28 and len(t6_support) == 16
    correction_seven = [
        poisson({monomial: field.one}, T) for monomial in s6_support
    ] + [poisson(S, {monomial: field.one}) for monomial in t6_support]

    indices = graded_indices(9, 3)
    defects = {
        point: order_seven_defect_at(
            point,
            presentation,
            field,
            structure,
            fifth_coefficients,
            fifth_kernel,
        )
        for point in indices
    }
    monomials = sorted(
        set().union(
            *(set(defect) for defect in defects.values()),
            *(set(column) for column in correction_seven),
        )
    )
    output_index = {monomial: index for index, monomial in enumerate(monomials)}
    rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(correction_seven)
    }
    reduced, pivots, nonzero = sdm_irref(rows)
    dual, _ = sdm_nullspace_from_rref(
        reduced, field.one, len(monomials), pivots, nonzero
    )
    functionals = [
        {
            monomials[index]: coefficient
            for index, coefficient in vector.items()
            if coefficient
        }
        for vector in dual
    ]
    values = {
        point: [
            sum(
                (
                    coefficient * defect.get(monomial, field.zero)
                    for monomial, coefficient in functional.items()
                ),
                field.zero,
            )
            for functional in functionals
        ]
        for point, defect in defects.items()
    }
    coefficients = newton_interpolation(values, indices, field)
    holdout = (4, 1, 2, 3, 2, 1, 3, 2, 1)
    holdout_defect = order_seven_defect_at(
        holdout,
        presentation,
        field,
        structure,
        fifth_coefficients,
        fifth_kernel,
    )
    assert set(holdout_defect).issubset(set(monomials))
    actual = [
        sum(
            (
                coefficient * holdout_defect.get(monomial, field.zero)
                for monomial, coefficient in functional.items()
            ),
            field.zero,
        )
        for functional in functionals
    ]
    assert newton_value(coefficients, holdout, field) == actual

    coefficient_rows = {}
    for row_index, exponent in enumerate(indices):
        row = {
            equation_index: coefficient
            for equation_index, coefficient in enumerate(coefficients[exponent])
            if coefficient
        }
        if row:
            coefficient_rows[row_index] = row
    _, equation_pivots, _ = sdm_irref(coefficient_rows)
    variable_names = tuple(f"u{index}" for index in (0, 2, 4, 6, 8)) + tuple(
        f"y{index}" for index in range(4)
    )
    polynomials = []
    for equation_index in equation_pivots:
        terms = []
        for exponent in indices:
            coefficient = coefficients[exponent][equation_index]
            if not coefficient:
                continue
            factors = [
                binomial_term(variable, degree)
                for variable, degree in zip(variable_names, exponent, strict=True)
                if degree
            ]
            term = coefficient_text(field, coefficient, "b")
            if factors:
                term += "*" + "*".join(factors)
            terms.append(term)
        polynomials.append("+".join(terms) if terms else "0")
    program = f"""ring r=(0,b),({','.join(variable_names)}),dp;
minpoly={minimum_text(minimum, 'b')};
option(redSB);
ideal I={','.join(polynomials)};
ideal G=std(I);
print("SIZE="+string(size(G)));
print("DIMENSION="+string(dim(G)));
print("REDUCED_ONE="+string(reduce(1,G)));
"""
    output = run_singular(program)
    singular_result = {}
    for line in output.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in {"SIZE", "DIMENSION", "REDUCED_ONE"}:
                singular_result[key] = value
    return {
        "fifth_correction_rank": 97,
        "fifth_correction_kernel_dimension": len(fifth_kernel),
        "fifth_correction_interpolation_points": len(graded_indices(5, 2)),
        "fifth_correction_holdout": list(fifth_holdout),
        "current_correction_support_sizes": {
            "S6": len(s6_support),
            "T6": len(t6_support),
        },
        "current_correction_rank": len(pivots),
        "output_dimension": len(monomials),
        "projected_equations": len(functionals),
        "independent_polynomial_generators": len(equation_pivots),
        "equation_degree": 3,
        "interpolation_points": len(indices),
        "independent_interpolation_check": list(holdout),
        "singular_result": singular_result,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--component",
        type=Path,
        default=Path(
            "artifacts/generated-results/degree_eight_order_five_rational_chart.json"
        ),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--order-five-only", action="store_true")
    args = parser.parse_args()

    component, minimum, field, sigma, tau = component_context(args.component)
    factors = sp.factor_list(minimum.as_expr())[1]
    assert factors == [(minimum.as_expr(), 1)]
    presentation = family_presentation(field, sigma, tau)
    pivot_chart = verify_reconstruction_chart(presentation, field)
    augmented_rank = column_rank(
        presentation["strong_columns"] + [presentation["constant"]]
    )
    (
        _,
        order_five_basis,
        structure,
        correction_rank,
        dual_rank,
        independent_order_five_generators,
    ) = exact_order_five(presentation, field, minimum, sigma)
    ranks = {
        "h3_rank": presentation["rank_three"],
        "h3_kernel_dimension": len(presentation["kernel_pairs"]),
        "h5_correction_rank": correction_rank,
        "h5_correction_kernel_dimension": len(
            presentation["correction_five"]
        )
        - correction_rank,
        "h5_strong_span_rank": pivot_chart["rank_over_residue_field"],
        "h5_augmented_rank": augmented_rank,
        "h5_output_dimension": len(presentation["output_support"]),
    }
    assert dual_rank == 172
    assert ranks == {
        "h3_rank": 168,
        "h3_kernel_dimension": 10,
        "h5_correction_rank": 97,
        "h5_correction_kernel_dimension": 4,
        "h5_strong_span_rank": 104,
        "h5_augmented_rank": 104,
        "h5_output_dimension": 269,
    }
    certificate = {
        "scope": (
            "exact nonlinear lifting in the parity-preserving inherited "
            "root-weight filtration; no unrestricted or DC_2 claim"
        ),
        "order_five_component": {
            "residue_degree": minimum.degree(),
            "minimal_polynomial": str(minimum.as_expr()),
            "factor_degrees_over_Q": component["exact_zero_scheme"][
                "factor_degrees_over_Q"
            ],
            "relative_ranks": ranks,
            "reconstruction_chart": pivot_chart,
        },
        "nonlinear_order_five_lift": {
            "projected_equations": 172,
            "independent_polynomial_generators": independent_order_five_generators,
            "groebner_basis_size": len(order_five_basis),
            "dimension_over_residue_field": 5,
            "shape": (
                "three affine-linear relations, two repeated-root quadratics, "
                "and one bridge relation"
            ),
            "reduced_lift_residue_degree": minimum.degree(),
            "reduced_scheme": "affine five-space over the residue field",
            "free_lower_lift_coordinates": ["z0", "z2", "z4", "z6", "z8"],
            "nilpotent_root_coordinates": ["z7", "z9"],
        },
    }
    if not args.order_five_only:
        order_seven = exact_order_seven(presentation, field, structure, minimum)
        certificate["order_seven"] = order_seven
        unit = order_seven["singular_result"] == {
            "SIZE": "1",
            "DIMENSION": "-1",
            "REDUCED_ONE": "0",
        }
        certificate["terminal_status"] = {
            "order_seven_unit_ideal": unit,
            "nonreduced_order_five_thickening": (
                "also excluded: a unit ideal modulo a nilpotent ideal lifts to a unit ideal"
                if unit
                else "requires a separate nilpotent-direction test"
            ),
            "ore_localized_quantization": (
                "inconsistent at hbar^7" if unit else "survives hbar^7"
            ),
            "conductor_gluing": (
                "not applicable: no hbar^7 local quantization"
                if unit
                else "required next"
            ),
            "weyl_relations_and_nonsurjectivity": (
                "not applicable: no quantized homomorphism survives"
                if unit
                else "required next"
            ),
        }
        certificate["root_at_infinity_weights"] = {
            "variables": {"X": 1, "Q": -1, "Z": -2},
            "classical": {"S": -2, "T": -1},
            "order_two": {"S2": 4, "T2": 5},
            "order_four": {"S4": 10, "T4": 11},
            "order_six": {"S6": 16, "T6": 17},
        }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")

    print(
        f"PASS: reconstructed the irreducible degree-{minimum.degree()} order-five component"
    )
    print("PASS: exact ranks over its residue field are 168, 97, 104, 104")
    print("PASS: the reduced nonlinear lift is affine five-space")
    if "order_seven" in certificate:
        print("ORDER7:", certificate["order_seven"]["singular_result"])


if __name__ == "__main__":
    main()

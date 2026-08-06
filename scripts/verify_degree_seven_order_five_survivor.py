#!/usr/bin/env python3
"""Exact nonlinear lift and order-seven test for the degree-seven survivor.

The characteristic-zero strong order-five locus was reconstructed separately
as one irreducible closed point of degree eight.  This verifier embeds that
point in its residue field, recomputes all relative ranks, and constructs the
genuine nonlinear order-five Kuranishi scheme.  Its quadratic equation is a
perfect square, so the reduced lift stays over the octic field.  The complete
order-seven calculation retains both free order-five correction directions.

All assertions concern only the parity-preserving, root-weight-homogeneous
PBW filtration used by the parent search.  They are not a DC_2 theorem.
"""

from __future__ import annotations

import argparse
from itertools import product
import json
from math import comb
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_quantum_residue import (
    add,
    column_rank,
    pi_power,
    poisson,
    scale,
    solve_affine,
    split_pair,
)
from screen_degree_seven_order_five_survivors import (
    projected_order_five_equations,
)
from reconstruct_degree_seven_order_five_zero_scheme import PIVOT_COLUMNS
from verify_degree_seven_relative_quantization_obstruction import (
    S4_SUPPORT,
    T4_SUPPORT,
    family_presentation,
    rank_record,
    weight_monomials,
)


X = sp.symbols("x")
SIGMA_POLYNOMIAL = sp.Poly(
    1687500 * X**8
    + 41047500 * X**7
    + 462666750 * X**6
    + 3259667250 * X**5
    + 15781954748 * X**4
    + 53969799492 * X**3
    + 126253770468 * X**2
    + 183369004011 * X
    + 119142437697,
    X,
    domain=QQ,
)
TAU_DENOMINATOR = 383180852409815403888
TAU_COEFFICIENTS_ASCENDING = (
    269809615785764223981,
    840178458581066844732,
    599309543339487802104,
    337020658576478999692,
    100324562410939815000,
    19181170262239043250,
    2282786985300315000,
    123670522877062500,
)


def evaluate_coefficients(field, element, coefficients):
    value = field.zero
    for coefficient in reversed(coefficients):
        value = value * element + field(coefficient)
    return value


def octic_context():
    field = QQ.alg_field_from_poly(SIGMA_POLYNOMIAL, alias="alpha")
    sigma = field.convert(field.ext)
    tau = -evaluate_coefficients(
        field, sigma, TAU_COEFFICIENTS_ASCENDING
    ) / field(TAU_DENOMINATOR)
    return field, sigma, tau


def coefficient_text(field, coefficient, alias: str) -> str:
    text = str(field.to_sympy(coefficient)).replace(str(field.ext), alias)
    text = re.sub(
        rf"{re.escape(alias)}\*\*(\d+)",
        lambda match: "*".join([alias] * int(match.group(1))),
        text,
    )
    return f"({text})"


def singular_polynomial(field, equation, variables, alias: str) -> str:
    rendered = []
    for exponent, coefficient in sorted(equation.items()):
        if not coefficient:
            continue
        term = coefficient_text(field, coefficient, alias)
        for variable, degree in zip(variables, exponent):
            if degree == 1:
                term += f"*{variable}"
            elif degree > 1:
                term += f"*{variable}^{degree}"
        rendered.append(term)
    return "+".join(rendered) if rendered else "0"


def run_singular(program: str):
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the exact survivor test")
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if completed.stderr.strip() or "   ?" in completed.stdout:
        raise AssertionError(completed.stdout + completed.stderr)
    return completed.stdout


def exact_order_five(presentation, field):
    equations, correction_rank, dual_rank = projected_order_five_equations(
        presentation, field
    )
    variables = tuple(f"z{index}" for index in range(8))
    ideal = ",".join(
        singular_polynomial(field, equation, variables, "a")
        for equation in equations
    )
    minimum = str(SIGMA_POLYNOMIAL.as_expr()).replace("x", "a").replace(
        "**", "^"
    )
    program = f"""ring r=(0,a),({','.join(variables)}),dp;
minpoly={minimum};
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
        "SIZE": "4",
        "DIMENSION": "4",
        "VDIM": "-1",
        "REDUCED_ONE": "1",
    }
    assert correction_rank == 61 and dual_rank == 117
    return equations, basis


def verify_reconstruction_chart(presentation, field):
    """Check that the reconstructed octic point is inside the Q pivot chart."""

    strong_columns = presentation["strong_columns"]
    monomials = sorted(
        set(presentation["constant"]).union(
            *(set(column) for column in strong_columns)
        )
    )
    pivot_monomials = monomials[:68]
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
    assert pivots == list(range(68))
    return {
        "pivot_columns": len(PIVOT_COLUMNS),
        "pivot_rows": len(pivot_monomials),
        "rank_over_octic_residue_field": len(pivots),
        "component_inside_chart": True,
    }


def coefficient_polynomial_to_field(expression, symbol, field, element):
    polynomial = sp.Poly(expression, symbol, domain=QQ)
    return sum(
        (
            field(coefficient.p) / field(coefficient.q) * element ** monomial[0]
            for monomial, coefficient in polynomial.terms()
        ),
        field.zero,
    )


def quadratic_root(quadratic_relation: str, field, sigma):
    """Return the unique reduced root of the nonreduced quadratic equation."""

    a_symbol, z_symbol = sp.symbols("a z7")
    expression = sp.sympify(
        quadratic_relation.replace("^", "**"),
        locals={"a": a_symbol, "z7": z_symbol},
    )
    polynomial = sp.Poly(expression, z_symbol)
    coefficients = [
        coefficient_polynomial_to_field(coefficient, a_symbol, field, sigma)
        for coefficient in polynomial.all_coeffs()
    ]
    assert len(coefficients) == 3
    leading, linear, constant = coefficients
    assert leading and linear**2 - field(4) * leading * constant == field.zero
    root = -linear / (field(2) * leading)
    assert leading * root**2 + linear * root + constant == field.zero
    return root


def triangular_relations(basis_text, field, sigma):
    a_symbol = sp.symbols("a")
    z_symbols = sp.symbols("z0:8")
    relations = {}
    for text in basis_text:
        expression = sp.sympify(
            text.replace("^", "**"),
            locals={"a": a_symbol, **{str(z): z for z in z_symbols}},
        )
        polynomial = sp.Poly(expression, *z_symbols)
        terms = {}
        for exponent, coefficient in polynomial.terms():
            terms[exponent] = coefficient_polynomial_to_field(
                coefficient, a_symbol, field, sigma
            )
        active = {
            index
            for exponent in terms
            for index, degree in enumerate(exponent)
            if degree
        }
        for dependent in (5, 3, 1):
            unit = tuple(1 if index == dependent else 0 for index in range(8))
            if dependent in active and unit in terms:
                assert all(sum(exponent) <= 1 for exponent in terms)
                relations[dependent] = terms
                break
    assert set(relations) == {1, 3, 5}
    return relations


def lower_lift_at(point, presentation, field, beta, relations):
    parameters = [field.zero] * 8
    parameters[0], parameters[2], parameters[4], parameters[6] = map(
        field, point
    )
    parameters[7] = beta
    zero = (0,) * 8
    for dependent in (5, 3, 1):
        terms = relations[dependent]
        unit = tuple(1 if index == dependent else 0 for index in range(8))
        value = terms.get(zero, field.zero)
        for exponent, coefficient in terms.items():
            if exponent == zero or exponent == unit:
                continue
            index = exponent.index(1)
            value += coefficient * parameters[index]
        parameters[dependent] = -value / terms[unit]

    s2, t2 = presentation["base_pair"]
    for coefficient, (basis_s, basis_t) in zip(
        parameters, presentation["kernel_pairs"], strict=True
    ):
        s2 = add(s2, basis_s, coefficient)
        t2 = add(t2, basis_t, coefficient)
    return s2, t2


def fifth_defect(presentation, s2, t2, field):
    S, T = presentation["S"], presentation["T"]
    result = poisson(s2, t2)
    result = add(result, pi_power(s2, T, 3), field.one / field(24))
    result = add(result, pi_power(S, t2, 3), field.one / field(24))
    return add(result, pi_power(S, T, 5), field.one / field(1920))


def graded_indices(variable_count: int, degree: int):
    return [
        exponent
        for total in range(degree + 1)
        for exponent in product(range(total + 1), repeat=variable_count)
        if sum(exponent) == total
    ]


def newton_interpolation(values, indices, field):
    coefficients = {}
    for alpha in indices:
        row = list(values[alpha])
        for gamma, earlier in coefficients.items():
            factor = 1
            for alpha_value, gamma_value in zip(alpha, gamma, strict=True):
                factor *= comb(alpha_value, gamma_value)
            if factor:
                row = [
                    value - field(factor) * coefficient
                    for value, coefficient in zip(row, earlier, strict=True)
                ]
        coefficients[alpha] = row
    return coefficients


def newton_value(coefficients, point, field):
    width = len(next(iter(coefficients.values())))
    result = [field.zero] * width
    for exponent, row in coefficients.items():
        factor = 1
        for value, degree in zip(point, exponent, strict=True):
            factor *= comb(value, degree)
        if factor:
            result = [
                value + field(factor) * coefficient
                for value, coefficient in zip(result, row, strict=True)
            ]
    return result


def fifth_correction_interpolation(presentation, field, beta, relations):
    indices = graded_indices(4, 2)
    values = {}
    expected_kernel = None
    width = len(presentation["correction_five"])
    for point in indices:
        s2, t2 = lower_lift_at(point, presentation, field, beta, relations)
        vector, kernel, rank = solve_affine(
            presentation["correction_five"],
            scale(fifth_defect(presentation, s2, t2, field), -field.one),
            field,
        )
        assert rank == 61 and len(kernel) == 2
        if expected_kernel is None:
            expected_kernel = kernel
        else:
            assert kernel == expected_kernel
        values[point] = [
            vector.get(index, field.zero) for index in range(width)
        ]
    coefficients = newton_interpolation(values, indices, field)
    holdout = (3, 1, 2, 4)
    s2, t2 = lower_lift_at(
        holdout, presentation, field, beta, relations
    )
    actual, kernel, rank = solve_affine(
        presentation["correction_five"],
        scale(fifth_defect(presentation, s2, t2, field), -field.one),
        field,
    )
    assert rank == 61 and kernel == expected_kernel
    actual_values = [actual.get(index, field.zero) for index in range(width)]
    assert newton_value(coefficients, holdout, field) == actual_values
    return coefficients, expected_kernel, holdout


def order_seven_defect_at(
    point,
    presentation,
    field,
    beta,
    relations,
    fifth_coefficients,
    fifth_kernel,
):
    lower_point = point[:4]
    s2, t2 = lower_lift_at(
        lower_point, presentation, field, beta, relations
    )
    correction = newton_value(fifth_coefficients, lower_point, field)
    for parameter, kernel_vector in zip(
        point[4:], fifth_kernel, strict=True
    ):
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
        result = add(
            result, pi_power(left, right, 3), field.one / field(24)
        )
    for left, right in ((s2, T), (S, t2)):
        result = add(
            result, pi_power(left, right, 5), field.one / field(1920)
        )
    return add(result, pi_power(S, T, 7), field.one / field(322560))


def binomial_term(variable: str, degree: int) -> str:
    if degree == 0:
        return ""
    factors = [variable] + [
        f"({variable}-{offset})" for offset in range(1, degree)
    ]
    denominator = 1
    for integer in range(2, degree + 1):
        denominator *= integer
    value = "*".join(factors)
    return f"({value}/{denominator})" if denominator > 1 else value


def exact_order_seven(
    presentation,
    field,
    beta,
    relations,
    minimum,
):
    fifth_coefficients, fifth_kernel, fifth_holdout = (
        fifth_correction_interpolation(
            presentation, field, beta, relations
        )
    )
    S, T = presentation["S"], presentation["T"]
    s6_support = weight_monomials(31, 1, 16)
    t6_support = weight_monomials(27, 0, 17)
    assert len(s6_support) == 14 and len(t6_support) == 6
    correction_seven = [
        poisson({monomial: field.one}, T) for monomial in s6_support
    ] + [poisson(S, {monomial: field.one}) for monomial in t6_support]
    assert column_rank(correction_seven) == 20

    indices = graded_indices(6, 3)
    defects = {
        point: order_seven_defect_at(
            point,
            presentation,
            field,
            beta,
            relations,
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
    output_index = {
        monomial: index for index, monomial in enumerate(monomials)
    }
    rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(correction_seven)
    }
    reduced, pivots, nonzero = sdm_irref(rows)
    dual, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(monomials),
        pivots,
        nonzero,
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
                    coefficient
                    * defect.get(monomial, field.zero)
                    for monomial, coefficient in functional.items()
                ),
                field.zero,
            )
            for functional in functionals
        ]
        for point, defect in defects.items()
    }
    coefficients = newton_interpolation(values, indices, field)
    holdout = (4, 1, 2, 3, 2, 1)
    holdout_defect = order_seven_defect_at(
        holdout,
        presentation,
        field,
        beta,
        relations,
        fifth_coefficients,
        fifth_kernel,
    )
    assert set(holdout_defect).issubset(set(monomials))
    actual = [
        sum(
            (
                coefficient
                * holdout_defect.get(monomial, field.zero)
                for monomial, coefficient in functional.items()
            ),
            field.zero,
        )
        for functional in functionals
    ]
    assert newton_value(coefficients, holdout, field) == actual

    variable_names = ("u0", "u2", "u4", "u6", "y0", "y1")
    polynomials = []
    for equation_index in range(len(functionals)):
        terms = []
        for exponent in indices:
            coefficient = coefficients[exponent][equation_index]
            if not coefficient:
                continue
            factors = [
                binomial_term(variable, degree)
                for variable, degree in zip(
                    variable_names, exponent, strict=True
                )
                if degree
            ]
            term = coefficient_text(field, coefficient, "b")
            if factors:
                term += "*" + "*".join(factors)
            terms.append(term)
        polynomials.append("+".join(terms) if terms else "0")
    minimum_text = str(minimum.as_expr()).replace(
        str(minimum.gens[0]), "b"
    ).replace("**", "^")
    program = f"""ring r=(0,b),({','.join(variable_names)}),dp;
minpoly={minimum_text};
option(redSB);
ideal I={','.join(polynomials)};
ideal G=std(I);
print("SIZE="+string(size(G)));
print("DIMENSION="+string(dim(G)));
print("REDUCED_ONE="+string(reduce(1,G)));
"""
    output = run_singular(program)
    values = {}
    for line in output.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in {"SIZE", "DIMENSION", "REDUCED_ONE"}:
                values[key] = value
    return {
        "fifth_correction_rank": 61,
        "fifth_correction_kernel_dimension": len(fifth_kernel),
        "fifth_correction_interpolation_points": len(
            graded_indices(4, 2)
        ),
        "fifth_correction_holdout": list(fifth_holdout),
        "current_correction_support_sizes": {
            "S6": len(s6_support),
            "T6": len(t6_support),
        },
        "current_correction_rank": len(pivots),
        "output_dimension": len(monomials),
        "projected_equations": len(functionals),
        "equation_degree": 3,
        "interpolation_points": len(indices),
        "independent_interpolation_check": list(holdout),
        "singular_result": values,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--order-five-only",
        action="store_true",
        help="stop after reconstructing the genuine nonlinear lift field",
    )
    args = parser.parse_args()

    assert sp.factor_list(SIGMA_POLYNOMIAL.as_expr())[1] == [
        (SIGMA_POLYNOMIAL.as_expr(), 1)
    ]
    octic, sigma, tau = octic_context()
    ranks = rank_record(octic, sigma, tau)
    assert ranks == {
        "h3_rank": 118,
        "h3_kernel_dimension": 8,
        "h5_correction_rank": 61,
        "h5_strong_span_rank": 68,
        "h5_augmented_rank": 68,
        "h5_output_dimension": 178,
    }
    octic_presentation = family_presentation(octic, sigma, tau)
    pivot_chart = verify_reconstruction_chart(octic_presentation, octic)
    _, order_five_basis = exact_order_five(octic_presentation, octic)
    quadratic = next(text for text in order_five_basis if "z7^2" in text)
    beta = quadratic_root(quadratic, octic, sigma)
    relations = triangular_relations(order_five_basis, octic, sigma)

    certificate = {
        "scope": (
            "exact nonlinear lifting in the parity-preserving inherited "
            "root-weight filtration; no unrestricted or DC_2 claim"
        ),
        "order_five_component": {
            "residue_degree": 8,
            "minimal_polynomial": str(SIGMA_POLYNOMIAL.as_expr()),
            "relative_ranks": ranks,
            "reconstruction_chart": pivot_chart,
        },
        "nonlinear_order_five_lift": {
            "projected_equations": 117,
            "groebner_basis_size": 4,
            "dimension_over_residue_field": 4,
            "shape": (
                "three affine-linear relations and one perfect-square "
                "quadratic relation"
            ),
            "reduced_lift_residue_degree": 8,
            "scheme_structure": (
                "a doubled affine four-space over the octic field"
            ),
            "free_lower_lift_coordinates": ["z0", "z2", "z4", "z6"],
        },
    }
    if not args.order_five_only:
        presentation = family_presentation(octic, sigma, tau)
        order_seven = exact_order_seven(
            presentation,
            octic,
            beta,
            relations,
            SIGMA_POLYNOMIAL,
        )
        certificate["order_seven"] = order_seven
        unit = order_seven["singular_result"] == {
            "SIZE": "1",
            "DIMENSION": "-1",
            "REDUCED_ONE": "0",
        }
        certificate["terminal_status"] = {
            "order_seven_unit_ideal": unit,
            "nonreduced_order_five_thickening": (
                "also excluded: a unit ideal modulo a nilpotent ideal "
                "lifts to a unit ideal"
                if unit
                else "requires a separate nilpotent-direction test"
            ),
            "ore_localized_quantization": (
                "inconsistent at hbar^7" if unit else "survives hbar^7"
            ),
            "root_at_infinity_valuations": (
                "all correction summands through hbar^6 recorded"
                if unit
                else "required beyond hbar^6"
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

    print("PASS: reconstructed the irreducible degree-eight order-five component")
    print("PASS: exact ranks over its residue field are 118, 61, 68, 68")
    print("PASS: genuine order-five lifts form a doubled affine four-space")
    print("PASS: the reduced nonlinear lift stays over the octic residue field")
    if "order_seven" in certificate:
        print("ORDER7:", certificate["order_seven"]["singular_result"])


if __name__ == "__main__":
    main()

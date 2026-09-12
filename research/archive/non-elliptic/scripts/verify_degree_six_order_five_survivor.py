#!/usr/bin/env python3
"""Exact reconstruction and order-seven closure of the degree-six survivor.

Modular scans of the two-parameter degree-six marked-root family find a
residue-degree-four component on which the known strong order-five cocycle
disappears.  This verifier reconstructs that component over Q, checks the
relative ranks over its quartic residue field, solves the genuine nonlinear
order-five Kuranishi equations, and proves that the complete inherited
order-seven equations generate the unit ideal.

The order-seven calculation retains the full affine four-space of
order-five lifts.  Consequently this component does not reach boundary
valuation gluing or a Weyl nonsurjectivity test.  The result is restricted to
the parity-preserving root-weight-homogeneous filtration used by the parent
degree-six calculation; it is not a result about DC_2.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_quantum_residue import (
    add,
    pi_power,
    poisson,
    scale,
    solve_affine,
    split_pair,
)
from verify_degree_six_relative_quantization_obstruction import (  # noqa: E402
    S2_SUPPORT,
    S4_SUPPORT,
    T2_SUPPORT,
    T4_SUPPORT,
    family_presentation,
    pairing,
    rank_record,
    weight_monomials,
)


MINIMAL_VARIABLE = sp.symbols("minimal_variable")
MINIMAL_POLYNOMIAL = sp.Poly(
    4 * MINIMAL_VARIABLE**4
    + 66 * MINIMAL_VARIABLE**3
    + 561 * MINIMAL_VARIABLE**2
    + 1260 * MINIMAL_VARIABLE
    + 900,
    MINIMAL_VARIABLE,
    domain=QQ,
)

EXPECTED_PRIME_POINTS = {
    17: [(2, 16), (9, 15)],
    19: [],
    23: [(2, 10), (8, 20), (13, 9), (18, 15)],
    29: [],
    31: [(14, 22), (25, 28)],
    37: [],
    41: [],
    139: [(56, 94), (62, 46), (75, 117), (138, 110)],
    167: [(6, 79), (25, 74), (96, 20), (107, 97)],
    197: [(94, 101), (95, 13), (102, 60), (185, 141)],
}


def tau_numerator(sigma):
    return 52 * sigma**3 + 498 * sigma**2 + 12693 * sigma + 13050


def rank_drop_quartic(sigma, tau):
    """Seven-prime reconstructed candidate for the strong rank-drop divisor."""

    return (
        2563 * sigma**4
        + 3954 * sigma**3 * tau
        + 2319 * sigma**2 * tau**2
        + 608 * sigma * tau**3
        + 60 * tau**4
        + 7240 * sigma**3
        + 7200 * sigma**2 * tau
        + 2280 * sigma * tau**2
        + 200 * tau**3
        + 6970 * sigma**2
        + 3400 * sigma * tau
        + 250 * tau**2
        + 2250 * sigma
    )


def component_prime_profile() -> dict[str, list[list[int]]]:
    profile = {}
    for prime, expected in EXPECTED_PRIME_POINTS.items():
        denominator = 4260 % prime
        points = []
        for sigma in range(prime):
            if int(MINIMAL_POLYNOMIAL.eval(sigma)) % prime:
                continue
            tau = (
                -tau_numerator(sigma) * pow(denominator, -1, prime)
            ) % prime
            points.append((sigma, tau))
        assert points == expected
        profile[str(prime)] = [list(point) for point in points]
    return profile


def algebraic_context():
    field = QQ.alg_field_from_poly(MINIMAL_POLYNOMIAL, alias="alpha")
    alpha = field.convert(field.ext)
    tau = -(
        field(52) * alpha**3
        + field(498) * alpha**2
        + field(12693) * alpha
        + field(13050)
    ) / field(4260)
    return field, alpha, tau


def field_polynomial(field, alpha, coefficients: tuple[int, int, int, int]):
    c3, c2, c1, c0 = coefficients
    return (
        field(c3) * alpha**3
        + field(c2) * alpha**2
        + field(c1) * alpha
        + field(c0)
    )


def order_five_relations(field, alpha) -> dict[str, object]:
    """Return the two exact affine-linear Kuranishi relations."""

    c3 = field(8512859013120)
    c4 = field(34051436052480)
    c5 = field_polynomial(
        field,
        alpha,
        (
            -7367026540544,
            -56394136944640,
            -75351648632832,
            3094945908916224,
        ),
    )
    c0 = field_polynomial(
        field,
        alpha,
        (
            -17965535935901400,
            -241898405077104300,
            -1786202422216954875,
            -363073339332742500,
        ),
    )
    d1 = field(7708564093560422400)
    d2 = field(25695213645201408000)
    d4 = field_polynomial(
        field,
        alpha,
        (
            2678585076744192000,
            20332895809830912000,
            24776911959883776000,
            -1050177685520646144000,
        ),
    )
    d5 = field_polynomial(
        field,
        alpha,
        (
            77255392911211626496,
            26038156777754722304,
            -3777441846376651554816,
            -124635212514728293171200,
        ),
    )
    d0 = field_polynomial(
        field,
        alpha,
        (
            728117581511193796616400,
            10300294838928867005910600,
            75334467310505203322549475,
            14139581352857038758393000,
        ),
    )
    return {
        "first": (c3, c4, c5, c0),
        "second": (d1, d2, d4, d5, d0),
    }


def singular_coefficient(field, coefficient) -> str:
    expression = str(field.to_sympy(coefficient)).replace("alpha", "a")
    # Singular parses powers of coefficient-field parameters as repeated
    # products, while powers of polynomial variables use the usual syntax.
    expression = expression.replace("a**3", "a*a*a")
    expression = expression.replace("a**2", "a*a")
    return f"({expression})"


def singular_polynomial(
    field,
    terms: dict[tuple[int, ...], object],
    variables: tuple[str, ...],
) -> str:
    rendered = []
    for exponent, coefficient in sorted(terms.items()):
        if not coefficient:
            continue
        monomial_factors = []
        for variable, degree in zip(variables, exponent):
            if degree == 1:
                monomial_factors.append(variable)
            elif degree > 1:
                monomial_factors.append(f"{variable}^{degree}")
        monomial = "*".join(monomial_factors)
        term = singular_coefficient(field, coefficient)
        if monomial:
            term += f"*{monomial}"
        rendered.append(term)
    return "+".join(rendered) if rendered else "0"


def run_singular(program: str) -> tuple[dict[str, str], str]:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the exact survivor certificate")
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if result.stderr.strip() or "   ?" in result.stdout:
        raise AssertionError(result.stdout + result.stderr)
    values = {}
    for line in result.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key == key.upper() and key.replace("_", "").isalnum():
            values[key] = value
    return values, hashlib.sha256(program.encode()).hexdigest()


def projected_order_five_equations(presentation, field):
    columns = presentation["correction_five"]
    constant = presentation["constant"]
    variations = presentation["lower_variations"]
    monomials = sorted(
        set(constant).union(
            *(set(column) for column in columns),
            *(set(column) for column in variations),
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
        for column_index, column in enumerate(columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(rows)
    dual, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(monomials),
        pivots,
        nonzero,
    )
    linear = variations[0:12:2]
    diagonal = variations[1:12:2]
    cross = variations[12:]
    assert len(linear) == len(diagonal) == 6
    assert len(cross) == 15

    equations = []
    for vector in dual:
        functional = {
            monomials[index]: coefficient
            for index, coefficient in vector.items()
            if coefficient
        }
        terms = {(0, 0, 0, 0, 0, 0): pairing(functional, constant, field)}
        for index in range(6):
            linear_exponent = tuple(
                1 if position == index else 0 for position in range(6)
            )
            square_exponent = tuple(
                2 if position == index else 0 for position in range(6)
            )
            terms[linear_exponent] = pairing(functional, linear[index], field)
            terms[square_exponent] = pairing(
                functional, diagonal[index], field
            )
        cross_index = 0
        for left in range(6):
            for right in range(left + 1, 6):
                exponent = tuple(
                    1 if position in (left, right) else 0
                    for position in range(6)
                )
                terms[exponent] = pairing(
                    functional,
                    cross[cross_index],
                    field,
                )
                cross_index += 1
        equations.append(
            {exponent: coefficient for exponent, coefficient in terms.items() if coefficient}
        )
    return equations, len(reduced), len(dual)


def verify_order_five_kuranishi(presentation, field, relations):
    equations, correction_rank, dual_rank = projected_order_five_equations(
        presentation, field
    )
    variables = ("z0", "z1", "z2", "z3", "z4", "z5")
    ideal = ",".join(
        singular_polynomial(field, equation, variables)
        for equation in equations
    )
    c3, c4, c5, c0 = relations["first"]
    d1, d2, d4, d5, d0 = relations["second"]
    first = (
        f"{singular_coefficient(field, c3)}*z3"
        f"+{singular_coefficient(field, c4)}*z4"
        f"+{singular_coefficient(field, c5)}*z5"
        f"+{singular_coefficient(field, c0)}"
    )
    second = (
        f"{singular_coefficient(field, d1)}*z1"
        f"+{singular_coefficient(field, d2)}*z2"
        f"+{singular_coefficient(field, d4)}*z4"
        f"+{singular_coefficient(field, d5)}*z5"
        f"+{singular_coefficient(field, d0)}"
    )
    program = f"""ring r=(0,a),({','.join(variables)}),dp;
minpoly=4*a*a*a*a+66*a*a*a+561*a*a+1260*a+900;
option(redSB);
ideal I={ideal};
ideal C={first},{second};
ideal GI=std(I);
ideal GC=std(C);
int input_in_candidate=1;
for (int k=1; k<=size(I); k++) {{ if (reduce(I[k],GC)!=0) {{ input_in_candidate=0; }} }}
int candidate_in_input=1;
for (int j=1; j<=size(C); j++) {{ if (reduce(C[j],GI)!=0) {{ candidate_in_input=0; }} }}
print("ORDER5_INPUT_IN_CANDIDATE="+string(input_in_candidate));
print("ORDER5_CANDIDATE_IN_INPUT="+string(candidate_in_input));
print("ORDER5_GB_SIZE="+string(size(GI)));
print("ORDER5_DIMENSION="+string(dim(GI)));
print("ORDER5_REDUCED_ONE="+string(reduce(1,GI)));
"""
    values, digest = run_singular(program)
    expected_values = {
        "ORDER5_INPUT_IN_CANDIDATE": "1",
        "ORDER5_CANDIDATE_IN_INPUT": "1",
        "ORDER5_GB_SIZE": "2",
        "ORDER5_DIMENSION": "4",
        "ORDER5_REDUCED_ONE": "1",
    }
    assert values == expected_values, (values, expected_values)
    return {
        "correction_rank": correction_rank,
        "dual_rank": dual_rank,
        "projected_equations": len(equations),
        "groebner_basis_size": 2,
        "dimension": 4,
        "singular_program_sha256": digest,
    }


def lower_lift_at(
    point: tuple[int, int, int, int],
    presentation,
    field,
    relations,
):
    u0, u2, u4, u5 = map(field, point)
    c3, _, c5, c0 = relations["first"]
    d1, d2, d4, d5, d0 = relations["second"]
    z3 = -field(4) * u4 - (c5 / c3) * u5 - c0 / c3
    z1 = -(
        (d2 / d1) * u2
        + (d4 / d1) * u4
        + (d5 / d1) * u5
        + d0 / d1
    )
    parameters = (u0, z1, u2, z3, u4, u5)
    s2, t2 = presentation["base_pair"]
    for coefficient, (basis_s, basis_t) in zip(
        parameters, presentation["kernel_pairs"]
    ):
        s2 = add(s2, basis_s, coefficient)
        t2 = add(t2, basis_t, coefficient)
    return s2, t2


def order_seven_values_at(
    point,
    presentation,
    field,
    relations,
    functionals,
    monomials,
):
    S = presentation["S"]
    T = presentation["T"]
    s2, t2 = lower_lift_at(point, presentation, field, relations)
    defect_five = poisson(s2, t2)
    defect_five = add(
        defect_five, pi_power(s2, T, 3), field.one / field(24)
    )
    defect_five = add(
        defect_five, pi_power(S, t2, 3), field.one / field(24)
    )
    defect_five = add(
        defect_five, pi_power(S, T, 5), field.one / field(1920)
    )
    fifth_vector, fifth_kernel, fifth_rank = solve_affine(
        presentation["correction_five"],
        scale(defect_five, -field.one),
        field,
    )
    assert fifth_rank == 34
    assert fifth_kernel == []
    s4, t4 = split_pair(fifth_vector, S4_SUPPORT, T4_SUPPORT)

    defect_seven = add(poisson(s2, t4), poisson(s4, t2))
    for left, right in ((s4, T), (S, t4), (s2, t2)):
        defect_seven = add(
            defect_seven,
            pi_power(left, right, 3),
            field.one / field(24),
        )
    for left, right in ((s2, T), (S, t2)):
        defect_seven = add(
            defect_seven,
            pi_power(left, right, 5),
            field.one / field(1920),
        )
    defect_seven = add(
        defect_seven,
        pi_power(S, T, 7),
        field.one / field(322560),
    )
    assert set(defect_seven).issubset(set(monomials))
    return [
        sum(
            (
                coefficient
                * defect_seven.get(monomial, field.zero)
                for monomial, coefficient in functional.items()
            ),
            field.zero,
        )
        for functional in functionals
    ]


def binomial_term(variable: str, degree: int) -> str:
    if degree == 0:
        return ""
    factors = [variable] + [f"({variable}-{offset})" for offset in range(1, degree)]
    denominator = 1
    for integer in range(2, degree + 1):
        denominator *= integer
    product_term = "*".join(factors)
    if denominator > 1:
        return f"({product_term}/{denominator})"
    return product_term


def verify_order_seven_unit(presentation, field, relations):
    S = presentation["S"]
    T = presentation["T"]
    s6_support = weight_monomials(24, 0, 16)
    assert s6_support == [
        (16, 0, 0),
        (17, 1, 0),
        (18, 2, 0),
        (19, 3, 0),
        (20, 4, 0),
    ]
    correction_seven = [
        poisson({monomial: field.one}, T) for monomial in s6_support
    ]

    # Reconstruct the complete support at one point, then quotient the five
    # current-correction columns once over the quartic field.
    S = presentation["S"]
    T = presentation["T"]
    s2, t2 = lower_lift_at((0, 0, 0, 0), presentation, field, relations)
    defect_five = poisson(s2, t2)
    defect_five = add(defect_five, pi_power(s2, T, 3), field.one / field(24))
    defect_five = add(defect_five, pi_power(S, t2, 3), field.one / field(24))
    defect_five = add(defect_five, pi_power(S, T, 5), field.one / field(1920))
    fifth_vector, fifth_kernel, fifth_rank = solve_affine(
        presentation["correction_five"],
        scale(defect_five, -field.one),
        field,
    )
    assert fifth_rank == 34 and fifth_kernel == []
    s4, t4 = split_pair(fifth_vector, S4_SUPPORT, T4_SUPPORT)
    first_defect = add(poisson(s2, t4), poisson(s4, t2))
    for left, right in ((s4, T), (S, t4), (s2, t2)):
        first_defect = add(
            first_defect, pi_power(left, right, 3), field.one / field(24)
        )
    for left, right in ((s2, T), (S, t2)):
        first_defect = add(
            first_defect, pi_power(left, right, 5), field.one / field(1920)
        )
    first_defect = add(
        first_defect, pi_power(S, T, 7), field.one / field(322560)
    )

    monomials = sorted(
        set(first_defect).union(*(set(column) for column in correction_seven))
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
    assert len(monomials) == 65
    assert len(reduced) == 5
    assert len(functionals) == 60

    indices = [
        exponent
        for total in range(4)
        for exponent in product(range(total + 1), repeat=4)
        if sum(exponent) == total
    ]
    values = {
        point: order_seven_values_at(
            point,
            presentation,
            field,
            relations,
            functionals,
            monomials,
        )
        for point in indices
    }

    # Multivariate Newton interpolation in the product binomial basis.
    newton = {}
    for alpha in indices:
        row = list(values[alpha])
        for gamma, coefficients in newton.items():
            factor = 1
            for alpha_value, gamma_value in zip(alpha, gamma):
                factor *= comb(alpha_value, gamma_value)
            if factor:
                row = [
                    value - field(factor) * coefficient
                    for value, coefficient in zip(row, coefficients)
                ]
        newton[alpha] = row

    # An extra point is not used for interpolation and certifies the degree
    # bound and reconstruction.
    test_point = (4, 1, 2, 3)
    predicted = [field.zero for _ in functionals]
    for alpha, coefficients in newton.items():
        factor = 1
        for value, degree in zip(test_point, alpha):
            factor *= comb(value, degree)
        if factor:
            predicted = [
                value + field(factor) * coefficient
                for value, coefficient in zip(predicted, coefficients)
            ]
    actual = order_seven_values_at(
        test_point,
        presentation,
        field,
        relations,
        functionals,
        monomials,
    )
    assert predicted == actual

    variable_names = ("u0", "u2", "u4", "u5")
    polynomials = []
    for equation_index in range(len(functionals)):
        terms = []
        for alpha in indices:
            coefficient = newton[alpha][equation_index]
            if not coefficient:
                continue
            factors = [
                binomial_term(variable, degree)
                for variable, degree in zip(variable_names, alpha)
                if degree
            ]
            term = singular_coefficient(field, coefficient)
            if factors:
                term += "*" + "*".join(factors)
            terms.append(term)
        polynomials.append("+".join(terms) if terms else "0")

    program = f"""ring r=(0,a),({','.join(variable_names)}),dp;
minpoly=4*a*a*a*a+66*a*a*a+561*a*a+1260*a+900;
ideal I={','.join(polynomials)};
option(redSB);
ideal G=std(I);
print("ORDER7_GB_SIZE="+string(size(G)));
print("ORDER7_DIMENSION="+string(dim(G)));
print("ORDER7_REDUCED_ONE="+string(reduce(1,G)));
"""
    singular_values, digest = run_singular(program)
    assert singular_values == {
        "ORDER7_GB_SIZE": "1",
        "ORDER7_DIMENSION": "-1",
        "ORDER7_REDUCED_ONE": "0",
    }
    return {
        "current_correction_support": [list(monomial) for monomial in s6_support],
        "current_correction_rank": 5,
        "output_dimension": 65,
        "projected_equations": 60,
        "equation_degree": 3,
        "interpolation_points": len(indices),
        "independent_interpolation_check": list(test_point),
        "groebner_basis": ["1"],
        "singular_program_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    assert sp.factor_list(MINIMAL_POLYNOMIAL.as_expr())[1] == [
        (MINIMAL_POLYNOMIAL.as_expr(), 1)
    ]
    prime_profile = component_prime_profile()

    sigma_symbol, tau_symbol = sp.symbols("sigma tau")
    tau_expression = -tau_numerator(sigma_symbol) / sp.Integer(4260)
    divisor_remainder = sp.rem(
        sp.together(
            rank_drop_quartic(sigma_symbol, tau_expression)
        ).as_numer_denom()[0],
        MINIMAL_POLYNOMIAL.as_expr().subs(
            MINIMAL_VARIABLE, sigma_symbol
        ),
        sigma_symbol,
        domain=QQ,
    )
    assert divisor_remainder == 0

    field, alpha, tau = algebraic_context()
    ranks = rank_record(field, alpha, tau)
    assert ranks == {
        "h3_rank": 77,
        "h3_kernel_dimension": 6,
        "h5_correction_rank": 34,
        "h5_strong_span_rank": 40,
        "h5_augmented_rank": 40,
        "h5_output_dimension": 110,
    }
    presentation = family_presentation(field, alpha, tau)
    order_three_columns = [
        poisson({monomial: field.one}, presentation["T"])
        for monomial in S2_SUPPORT
    ]
    order_three_columns += [
        poisson(presentation["S"], {monomial: field.one})
        for monomial in T2_SUPPORT
    ]
    order_three_rhs = scale(
        pi_power(presentation["S"], presentation["T"], 3),
        -field.one / field(24),
    )
    base_vector, rebuilt_kernel, rebuilt_rank = solve_affine(
        order_three_columns,
        order_three_rhs,
        field,
    )
    assert rebuilt_rank == 77
    assert len(rebuilt_kernel) == 6
    presentation["base_pair"] = split_pair(
        base_vector, S2_SUPPORT, T2_SUPPORT
    )
    relations = order_five_relations(field, alpha)
    order_five = verify_order_five_kuranishi(
        presentation, field, relations
    )
    order_seven = verify_order_seven_unit(presentation, field, relations)

    certificate = {
        "scope": (
            "exact quartic reconstruction and closure through hbar^7 in the "
            "parity-preserving root-weight-homogeneous degree-six filtration"
        ),
        "component": {
            "minimal_polynomial": str(MINIMAL_POLYNOMIAL.as_expr()),
            "tau_relation": (
                "4260*tau+52*sigma^3+498*sigma^2+12693*sigma+13050"
            ),
            "residue_degree": 4,
            "prime_point_profile": prime_profile,
            "contained_in_reconstructed_rank_drop_quartic": True,
        },
        "relative_ranks": ranks,
        "order_five_kuranishi": order_five,
        "order_seven": order_seven,
        "root_at_infinity_weights": {
            "variables": {"X": 1, "Q": -1, "Z": -2},
            "classical": {"S": -2, "T": -1},
            "order_two": {"S2": 4, "T2": 5},
            "order_four": {"S4": 10, "T4": 11},
            "order_six": {"S6": 16, "T6": None},
        },
        "terminal_status": {
            "ore_localized_quantization": "inconsistent at hbar^7",
            "root_at_infinity_correction_valuations": "not reached beyond the listed filtration weights",
            "conductor_gluing": "not reached",
            "weyl_relations": "not reached",
            "nonsurjectivity": "not reached",
        },
    }
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")

    print("PASS: reconstructed the residue-degree-four component over Q")
    print("PASS: exact ranks are h3=77, h5 strong=40, h5 augmented=40")
    print("PASS: the genuine order-five lift scheme is affine 4-space")
    print("PASS: 60 exact order-seven equations generate the unit ideal")
    print("CONCLUSION: the multi-prime component is obstructed at hbar^7")
    print("SCOPE: no boundary gluing or Weyl nonsurjectivity stage is reached")


if __name__ == "__main__":
    main()

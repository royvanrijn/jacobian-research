#!/usr/bin/env python3
"""Exact sparse searches for short Strong Factorial counterexamples.

The standard factorial functional is

    L(U_1^a_1 ... U_n^a_n) = a_1! ... a_n!.

Two finite support spaces are exhausted.

1. Every three-monomial support in C[x,y] whose monomials have total degree
   at most six.  After normalizing the first coefficient, L(f)=0 eliminates
   the third coefficient.  Exact univariate gcds then decide whether
   L(f^2)=L(f^3)=0 has a solution with all three coefficients nonzero.
2. Every pair of nontrivial monomial orbits of total degree at most six
   under (U_0,U_1,U_2,U_3) -> (U_2,U_3,U_0,U_1).  For

       f=(M-sigma(M)) + c(N-sigma(N)),

   all odd moments vanish by symmetry.  Exact gcds decide whether the second
   and fourth moments have a common nonzero coefficient c.
3. For every homogeneous binary form of degree d=1,2,3,4, projective
   Groebner charts decide the least initial moment cutoff.  The first d+1
   moments have no nonzero common zero, while the first d do.

A survivor in either search would violate the Strong Factorial Conjecture.
The script uses exact rational arithmetic and writes a finite-search artifact;
it makes no assertion outside the displayed support bounds.
"""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
from pathlib import Path
from typing import Iterator

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "sparse_factorial_moment_frontier.json"
)
MAX_TOTAL_DEGREE = 6
THREE_TERM_LONG_CHECK = 12
PAIRED_LONG_CHECK = 20
MAX_BINARY_FORM_DEGREE = 4

Exponent = tuple[int, ...]
RationalPolynomial = dict[int, Fraction]
IntegerPolynomial = dict[int, int]

A = sp.symbols("a")
C = sp.symbols("c")


def weak_compositions(total: int, parts: int) -> Iterator[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first,) + tail


def monomials(variable_count: int, maximum_total_degree: int) -> list[Exponent]:
    result = []
    for total_degree in range(maximum_total_degree + 1):
        result.extend(weak_compositions(total_degree, variable_count))
    return result


def multinomial(counts: tuple[int, ...]) -> int:
    total = sum(counts)
    result = factorial(total)
    for count in counts:
        result //= factorial(count)
    return result


def exponent_sum(support: tuple[Exponent, ...], counts: tuple[int, ...]) -> Exponent:
    return tuple(
        sum(count * exponent[coordinate] for count, exponent in zip(counts, support))
        for coordinate in range(len(support[0]))
    )


def factorial_weight(exponent: Exponent) -> int:
    result = 1
    for component in exponent:
        result *= factorial(component)
    return result


def add_rational_term(
    polynomial: RationalPolynomial, degree: int, coefficient: Fraction
) -> None:
    polynomial[degree] = polynomial.get(degree, Fraction(0)) + coefficient
    if polynomial[degree] == 0:
        del polynomial[degree]


def three_term_substituted_moment(support: tuple[Exponent, ...], order: int) -> RationalPolynomial:
    """Return L((M0+a*M1+b(a)*M2)^order) after eliminating L(f)=0."""
    linear_values = tuple(factorial_weight(exponent) for exponent in support)
    b_constant = Fraction(-linear_values[0], linear_values[2])
    b_linear = Fraction(-linear_values[1], linear_values[2])
    result: RationalPolynomial = {}
    for counts in weak_compositions(order, 3):
        k0, k1, k2 = counts
        scalar = Fraction(
            multinomial(counts) * factorial_weight(exponent_sum(support, counts))
        )
        for b_linear_count in range(k2 + 1):
            coefficient = (
                scalar
                * comb(k2, b_linear_count)
                * b_constant ** (k2 - b_linear_count)
                * b_linear**b_linear_count
            )
            add_rational_term(result, k1 + b_linear_count, coefficient)
    return result


def rational_poly(
    coefficients: RationalPolynomial, variable: sp.Symbol
) -> sp.Poly:
    expression = sum(
        sp.Rational(coefficient.numerator, coefficient.denominator)
        * variable**degree
        for degree, coefficient in coefficients.items()
    )
    return sp.Poly(expression, variable, domain=sp.QQ)


def integer_poly(coefficients: IntegerPolynomial, variable: sp.Symbol) -> sp.Poly:
    expression = sum(
        sp.Integer(coefficient) * variable**degree
        for degree, coefficient in coefficients.items()
    )
    return sp.Poly(expression, variable, domain=sp.QQ)


def nonexcluded_irreducible_factors(
    left: sp.Poly, right: sp.Poly, exclusion: sp.Poly
) -> list[sp.Poly]:
    common = sp.gcd(left, right)
    if common.is_zero or common.degree() <= 0:
        return []
    squarefree = common.sqf_part().monic()
    excluded = sp.gcd(squarefree, exclusion).monic()
    valid = squarefree.exquo(excluded).monic()
    if valid.degree() <= 0:
        return []
    return [
        factor.monic()
        for factor, _multiplicity in sp.factor_list(valid)[1]
        if factor.degree() > 0
    ]


def poly_as_string(polynomial: sp.Poly) -> str:
    return str(sp.factor(polynomial.as_expr()))


def first_nonzero_three_term_moment(
    support: tuple[Exponent, ...], factor: sp.Poly
) -> tuple[int | None, str | None]:
    for order in range(4, THREE_TERM_LONG_CHECK + 1):
        moment = rational_poly(three_term_substituted_moment(support, order), A)
        remainder = moment.rem(factor)
        if not remainder.is_zero:
            return order, poly_as_string(remainder)
    return None, None


def search_three_term_supports() -> dict:
    available_monomials = monomials(2, MAX_TOTAL_DEGREE)
    support_count = 0
    survivors = []
    for raw_support in combinations(available_monomials, 3):
        support = tuple(raw_support)
        support_count += 1
        second = rational_poly(three_term_substituted_moment(support, 2), A)
        third = rational_poly(three_term_substituted_moment(support, 3), A)

        l0, l1, _l2 = (
            factorial_weight(exponent) for exponent in support
        )
        exclusion = sp.Poly(A * (l0 + l1 * A), A, domain=sp.QQ)
        for factor in nonexcluded_irreducible_factors(second, third, exclusion):
            first_nonzero, remainder = first_nonzero_three_term_moment(
                support, factor
            )
            survivors.append(
                {
                    "support": [list(exponent) for exponent in support],
                    "coefficient_parameter": "a",
                    "third_coefficient": (
                        f"-({l0}+{l1}*a)/"
                        f"{factorial_weight(support[2])}"
                    ),
                    "common_factor": poly_as_string(factor),
                    "first_nonzero_checked_order": first_nonzero,
                    "first_nonzero_remainder_mod_factor": remainder,
                }
            )
    assert support_count == comb(len(available_monomials), 3)
    return {
        "variables": 2,
        "maximum_total_degree": MAX_TOTAL_DEGREE,
        "monomial_count": len(available_monomials),
        "support_count": support_count,
        "normalization": "coefficient of the first listed monomial is 1",
        "nonzero_coefficient_saturation": "a*b != 0",
        "moments_imposed": [1, 2, 3],
        "survivors": survivors,
    }


def involution(exponent: Exponent) -> Exponent:
    assert len(exponent) == 4
    return (exponent[2], exponent[3], exponent[0], exponent[1])


def nontrivial_involution_orbits() -> list[tuple[Exponent, Exponent]]:
    seen: set[Exponent] = set()
    result = []
    for exponent in monomials(4, MAX_TOTAL_DEGREE):
        partner = involution(exponent)
        if exponent == partner or exponent in seen:
            continue
        representative = min(exponent, partner)
        other = max(exponent, partner)
        result.append((representative, other))
        seen.add(representative)
        seen.add(other)
    result.sort()
    return result


def paired_moment(
    first_orbit: tuple[Exponent, Exponent],
    second_orbit: tuple[Exponent, Exponent],
    order: int,
) -> IntegerPolynomial:
    support = (
        first_orbit[0],
        first_orbit[1],
        second_orbit[0],
        second_orbit[1],
    )
    result: IntegerPolynomial = {}
    for counts in weak_compositions(order, 4):
        k0, k1, k2, k3 = counts
        coefficient = (
            (-1) ** (k1 + k3)
            * multinomial(counts)
            * factorial_weight(exponent_sum(support, counts))
        )
        c_degree = k2 + k3
        result[c_degree] = result.get(c_degree, 0) + coefficient
        if result[c_degree] == 0:
            del result[c_degree]
    return result


def symbolic_moment(
    support: tuple[Exponent, ...],
    coefficients: tuple[sp.Symbol, ...],
    order: int,
) -> sp.Expr:
    result = 0
    for counts in weak_compositions(order, len(support)):
        coefficient_monomial = sp.prod(
            coefficient**count
            for coefficient, count in zip(coefficients, counts)
        )
        result += (
            multinomial(counts)
            * factorial_weight(exponent_sum(support, counts))
            * coefficient_monomial
        )
    return sp.expand(result)


def is_unit_groebner(basis: sp.GroebnerBasis) -> bool:
    return (
        len(basis.polys) == 1
        and basis.polys[0].total_degree() == 0
        and basis.polys[0].LC() != 0
    )


def binary_homogeneous_cutoffs() -> dict:
    results = []
    for degree in range(1, MAX_BINARY_FORM_DEGREE + 1):
        support = tuple((degree - index, index) for index in range(degree + 1))
        coefficients = sp.symbols(f"q{degree}_0:{degree + 1}")
        moments = [
            symbolic_moment(support, coefficients, order)
            for order in range(1, degree + 2)
        ]

        unit_charts = []
        for chart in range(degree + 1):
            variables = tuple(
                coefficient
                for index, coefficient in enumerate(coefficients)
                if index != chart
            )
            chart_moments = [
                moment.subs(coefficients[chart], 1) for moment in moments
            ]
            basis = sp.groebner(
                chart_moments,
                *variables,
                order="grevlex",
                domain=sp.QQ,
            )
            assert is_unit_groebner(basis)
            unit_charts.append(chart)

        sharp_variables = tuple(coefficients[1:])
        sharp_basis = sp.groebner(
            [moment.subs(coefficients[0], 1) for moment in moments[:degree]],
            *sharp_variables,
            order="lex",
            domain=sp.QQ,
        )
        assert not is_unit_groebner(sharp_basis)
        assert sharp_basis.is_zero_dimensional
        eliminant = sharp_basis.polys[-1]
        assert eliminant.as_expr().free_symbols <= {sharp_variables[-1]}

        results.append(
            {
                "degree": degree,
                "coefficient_count": degree + 1,
                "sufficient_initial_moments": degree + 1,
                "projective_unit_charts": unit_charts,
                "sharp_initial_zero_count": degree,
                "sharp_chart": f"{coefficients[0]}=1",
                "sharp_chart_zero_dimensional": True,
                "sharp_eliminant_degree": eliminant.degree(sharp_variables[-1]),
                "sharp_eliminant": poly_as_string(eliminant),
            }
        )

    # Small explicit sharp witnesses used as transcription checks.
    sqrt15 = sp.sqrt(15)
    quadratic_coefficients = (
        sp.Integer(1),
        -sp.Rational(5, 2) + sp.I * sqrt15 / 2,
        sp.Rational(1, 4) - sp.I * sqrt15 / 4,
    )
    quadratic_support = ((2, 0), (1, 1), (0, 2))
    quadratic_moments = [
        sp.simplify(
            symbolic_moment(
                quadratic_support, quadratic_coefficients, order
            )
        )
        for order in range(1, 4)
    ]
    assert quadratic_moments[:2] == [0, 0]
    assert quadratic_moments[2] == 180 + 108 * sp.I * sqrt15

    sqrt21 = sp.sqrt(21)
    cubic_parameter = 6 + sp.I * sqrt21
    cubic_coefficients = (
        sp.Integer(1),
        -cubic_parameter,
        cubic_parameter,
        sp.Integer(-1),
    )
    cubic_support = ((3, 0), (2, 1), (1, 2), (0, 3))
    cubic_moments = [
        sp.simplify(
            symbolic_moment(cubic_support, cubic_coefficients, order)
        )
        for order in range(1, 5)
    ]
    assert cubic_moments[:3] == [0, 0, 0]
    assert cubic_moments[3] == 252564480 - 42301440 * sp.I * sqrt21

    return {
        "variables": 2,
        "homogeneous_degrees": [1, MAX_BINARY_FORM_DEGREE],
        "results": results,
        "explicit_sharp_witnesses": {
            "degree_2": {
                "polynomial": (
                    "x^2+(-5/2+i*sqrt(15)/2)xy"
                    "+(1/4-i*sqrt(15)/4)y^2"
                ),
                "moments": [
                    "0",
                    "0",
                    "180+108*i*sqrt(15)",
                ],
            },
            "degree_3": {
                "polynomial": (
                    "x^3-(6+i*sqrt(21))x^2y"
                    "+(6+i*sqrt(21))xy^2-y^3"
                ),
                "moments": [
                    "0",
                    "0",
                    "0",
                    "252564480-42301440*i*sqrt(21)",
                ],
            },
        },
    }


def first_nonzero_paired_moment(
    first_orbit: tuple[Exponent, Exponent],
    second_orbit: tuple[Exponent, Exponent],
    factor: sp.Poly,
) -> tuple[int | None, str | None]:
    for order in range(6, PAIRED_LONG_CHECK + 1, 2):
        moment = integer_poly(paired_moment(first_orbit, second_orbit, order), C)
        remainder = moment.rem(factor)
        if not remainder.is_zero:
            return order, poly_as_string(remainder)
    return None, None


def search_paired_supports() -> dict:
    orbits = nontrivial_involution_orbits()
    pair_count = 0
    survivors = []
    for first_orbit, second_orbit in combinations(orbits, 2):
        pair_count += 1
        assert paired_moment(first_orbit, second_orbit, 1) == {}
        assert paired_moment(first_orbit, second_orbit, 3) == {}

        second = integer_poly(paired_moment(first_orbit, second_orbit, 2), C)
        fourth = integer_poly(paired_moment(first_orbit, second_orbit, 4), C)
        exclusion = sp.Poly(C, C, domain=sp.QQ)
        for factor in nonexcluded_irreducible_factors(second, fourth, exclusion):
            first_nonzero, remainder = first_nonzero_paired_moment(
                first_orbit, second_orbit, factor
            )
            survivors.append(
                {
                    "first_orbit": [
                        list(exponent) for exponent in first_orbit
                    ],
                    "second_orbit": [
                        list(exponent) for exponent in second_orbit
                    ],
                    "coefficient_parameter": "c",
                    "common_factor": poly_as_string(factor),
                    "first_nonzero_checked_even_order": first_nonzero,
                    "first_nonzero_remainder_mod_factor": remainder,
                }
            )

    assert pair_count == comb(len(orbits), 2)
    return {
        "variables": 4,
        "involution": "(U0,U1,U2,U3) -> (U2,U3,U0,U1)",
        "maximum_total_degree": MAX_TOTAL_DEGREE,
        "nontrivial_orbit_count": len(orbits),
        "support_count": pair_count,
        "family": "(M-sigma(M))+c*(N-sigma(N))",
        "nonzero_coefficient_saturation": "c != 0",
        "automatic_zero_moments": "every odd order",
        "even_moments_imposed": [2, 4],
        "survivors": survivors,
    }


def main() -> None:
    three_term = search_three_term_supports()
    paired = search_paired_supports()
    binary_cutoffs = binary_homogeneous_cutoffs()

    artifact = {
        "format": "sparse-factorial-moment-frontier-v1",
        "field": "complex numbers",
        "factorial_functional": "L(U^alpha)=product_i alpha_i!",
        "searches": {
            "three_term_binary": three_term,
            "four_term_involution_paired": paired,
            "binary_homogeneous_cutoffs": binary_cutoffs,
        },
        "long_check_bounds": {
            "three_term": THREE_TERM_LONG_CHECK,
            "four_term_involution_paired": PAIRED_LONG_CHECK,
        },
        "scope": (
            "finite exact support searches only; absence of survivors is not "
            "a global Strong Factorial or Factorial theorem"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print(
        "PASS three-term binary factorial search: "
        f"{three_term['support_count']} supports, "
        f"{len(three_term['survivors'])} survivors"
    )
    print(
        "PASS four-term paired factorial search: "
        f"{paired['support_count']} supports, "
        f"{len(paired['survivors'])} survivors"
    )
    print("PASS binary homogeneous factorial cutoffs: degrees 1..4")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

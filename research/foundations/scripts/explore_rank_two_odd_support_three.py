#!/usr/bin/env python3
"""Exact coordinate-chart atlas for support-three odd quantization branches.

This script restricts the 41 quadratic Maurer--Cartan equations to every
coordinate projective plane in the 38-dimensional essential quotient.  The
finite-field search is retained as a profile only.  Completeness comes from
an independent characteristic-zero saturation of every chart not already
excluded by an exact monomial equation.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations

import sympy as sp
from sympy.polys.domains import GF

from explore_rank_two_odd_mixed_quantization import (
    essential_problem,
    modular_essential_pairs,
    rational,
)
from explore_rank_two_odd_quantization import second_obstruction


PRIME = 31


def quadratic_equations_mod_prime(prime: int = PRIME):
    field = GF(prime)
    S, T, pairs = modular_essential_pairs(field)
    (
        _,
        parameter_monomials,
        _,
        equation_pivots,
        _,
        _,
        _,
        sparse_equations,
    ) = second_obstruction(S, T, pairs, field)
    assert len(equation_pivots) == len(sparse_equations) == 41
    equations = []
    for equation in sparse_equations.values():
        coefficients = {
            parameter_monomials[column]: int(value) % prime
            for column, value in equation.items()
            if value
        }
        equations.append(coefficients)
    return equations


def restrict_equation(equation, triple, prime: int = PRIME):
    """Restrict to x_i=1, x_j=r, x_k=s in six-monomial coordinates."""

    left, middle, right = triple
    positions = {left: 0, middle: 1, right: 2}
    # 1, r, s, r^2, r*s, s^2
    monomial_slot = {
        (0, 0): 0,
        (0, 1): 1,
        (0, 2): 2,
        (1, 1): 3,
        (1, 2): 4,
        (2, 2): 5,
    }
    result = [0] * 6
    for (first, second), coefficient in equation.items():
        if first not in positions or second not in positions:
            continue
        pair = tuple(sorted((positions[first], positions[second])))
        result[monomial_slot[pair]] += coefficient
    return tuple(value % prime for value in result)


def evaluate(restriction, r_value, s_value, prime: int = PRIME):
    values = (
        1,
        r_value,
        s_value,
        r_value * r_value,
        r_value * s_value,
        s_value * s_value,
    )
    return sum(
        coefficient * value
        for coefficient, value in zip(restriction, values)
    ) % prime


def support_three_search():
    equations = quadratic_equations_mod_prime()
    torus_points = [
        (r_value, s_value)
        for r_value in range(1, PRIME)
        for s_value in range(1, PRIME)
    ]
    survivors = {}
    for triple in combinations(range(38), 3):
        restrictions = {
            restrict_equation(equation, triple)
            for equation in equations
        }
        restrictions.discard((0, 0, 0, 0, 0, 0))
        if not restrictions:
            survivors[triple] = None
            continue
        # A single Laurent monomial cannot vanish on the coordinate torus.
        if any(sum(value != 0 for value in restriction) == 1
               for restriction in restrictions):
            continue
        candidates = torus_points
        for restriction in sorted(
            restrictions,
            key=lambda row: sum(value != 0 for value in row),
        ):
            candidates = [
                (r_value, s_value)
                for r_value, s_value in candidates
                if evaluate(restriction, r_value, s_value) == 0
            ]
            if not candidates:
                break
        if candidates:
            survivors[triple] = candidates
    return survivors


def exact_support_three_audit():
    """Classify every coordinate chart after saturating by its coordinates."""

    _, _, _, _, parameter_monomials, sparse_equations = essential_problem()
    exact_equations = []
    for equation in sparse_equations.values():
        exact_equations.append(
            {
                parameter_monomials[column]: rational(value)
                for column, value in equation.items()
                if value
            }
        )

    r_symbol, s_symbol, inverse = sp.symbols("r s inverse")
    contained = []
    monomial_empty = []
    saturated_empty = []
    curves = {}
    points = {}
    for triple in combinations(range(38), 3):
        positions = {
            triple[0]: 0,
            triple[1]: 1,
            triple[2]: 2,
        }
        restricted_rows = []
        for equation in exact_equations:
            row = {}
            for (first, second), coefficient in equation.items():
                if first not in positions or second not in positions:
                    continue
                monomial = tuple(sorted((positions[first], positions[second])))
                row[monomial] = row.get(monomial, sp.S.Zero) + coefficient
            row = {
                monomial: coefficient
                for monomial, coefficient in row.items()
                if coefficient
            }
            if row:
                restricted_rows.append(row)
        if not restricted_rows:
            contained.append(triple)
            continue
        if any(len(row) == 1 for row in restricted_rows):
            monomial_empty.append(triple)
            continue

        coordinate_values = {
            triple[0]: sp.S.One,
            triple[1]: r_symbol,
            triple[2]: s_symbol,
        }
        restrictions = []
        for row in restricted_rows:
            expression = sp.S.Zero
            for (first_position, second_position), coefficient in row.items():
                first = triple[first_position]
                second = triple[second_position]
                expression += (
                    coefficient
                    * coordinate_values[first]
                    * coordinate_values[second]
                )
            restrictions.append(expression)
        basis = sp.groebner(
            restrictions + [inverse * r_symbol * s_symbol - 1],
            inverse,
            s_symbol,
            r_symbol,
            order="lex",
            domain=sp.QQ,
        )
        if len(basis.polys) == 1 and basis.polys[0].as_expr() == 1:
            saturated_empty.append(triple)
            continue
        elimination = tuple(
            sp.factor(polynomial.as_expr())
            for polynomial in basis.polys
            if inverse not in polynomial.as_expr().free_symbols
        )
        if basis.is_zero_dimensional:
            points[triple] = elimination
        else:
            curves[triple] = elimination
    return contained, monomial_empty, saturated_empty, curves, points


def curve_degree_profile(curves):
    r_symbol, s_symbol = sp.symbols("r s")
    profile = Counter()
    for equations in curves.values():
        common = sp.Poly(equations[0], s_symbol, r_symbol, domain=sp.QQ)
        for equation in equations[1:]:
            common = sp.gcd(
                common,
                sp.Poly(equation, s_symbol, r_symbol, domain=sp.QQ),
            )
        profile[common.total_degree()] += 1
    return sorted(profile.items())


def point_length_profile(points):
    r_symbol = sp.Symbol("r")
    profile = Counter()
    for equations in points.values():
        univariate = [
            sp.Poly(equation, r_symbol, domain=sp.QQ)
            for equation in equations
            if equation.free_symbols <= {r_symbol}
        ]
        profile[max(polynomial.degree() for polynomial in univariate)] += 1
    return sorted(profile.items())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print the coordinate triples in every nonempty exact class",
    )
    arguments = parser.parse_args()

    modular_survivors = support_three_search()
    contained_planes = [
        triple
        for triple, points_on_triple in modular_survivors.items()
        if points_on_triple is None
    ]
    modular_isolated = {
        triple: points
        for triple, points in modular_survivors.items()
        if points is not None
    }
    print(f"prime={PRIME}")
    print(f"coordinate planes contained in the quadratic locus={len(contained_planes)}")
    print(
        "other support-three triples with modular points="
        f"{len(modular_isolated)}"
    )
    point_count_profile = {}
    for points_on_triple in modular_isolated.values():
        point_count_profile[len(points_on_triple)] = (
            point_count_profile.get(len(points_on_triple), 0) + 1
        )
    print(f"modular point-count profile={sorted(point_count_profile.items())}")

    (
        exact_contained,
        monomial_empty,
        saturated_empty,
        curves,
        points,
    ) = exact_support_three_audit()
    assert exact_contained == contained_planes
    assert len(exact_contained) == 66
    assert len(monomial_empty) == 7924
    assert len(saturated_empty) == 268
    assert len(curves) == 149
    assert len(points) == 29
    assert (
        len(exact_contained)
        + len(monomial_empty)
        + len(saturated_empty)
        + len(curves)
        + len(points)
        == 8436
    )
    print(f"exact monomial-empty charts={len(monomial_empty)}")
    print(f"other exact empty charts after saturation={len(saturated_empty)}")
    print(f"positive-dimensional exact-support-three charts={len(curves)}")
    degree_profile = curve_degree_profile(curves)
    assert degree_profile == [(1, 52), (2, 97)]
    print(f"curve-degree profile={degree_profile}")
    print(f"zero-dimensional exact-support-three charts={len(points)}")
    length_profile = point_length_profile(points)
    assert length_profile == [(2, 1), (3, 19), (4, 9)]
    print(f"zero-dimensional scheme-length profile={length_profile}")
    print(
        "total zero-dimensional scheme length="
        f"{sum(length * count for length, count in length_profile)}"
    )
    if arguments.verbose:
        print(f"positive-dimensional triples={sorted(curves)}")
        print(f"zero-dimensional triples={sorted(points)}")


if __name__ == "__main__":
    main()

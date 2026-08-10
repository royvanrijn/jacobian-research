#!/usr/bin/env python3
"""Verify the chartwise Rees reduction for the quartic weighted graph.

The full graph ideal is too large for direct normalization.  This checker
uses the geometry of its source base triangle instead:

* over the function field of each open edge, it computes the complete
  two-dimensional base-point cluster by quadratic transforms;
* it derives the colength of the integrally closed ideal from that cluster;
* at the three vertices, it computes the compact Newton facets and their
  non-monomial face colors; and
* it principalizes every generic face color in a two-variable chart,
  including the degree-five algebraic color.

This is an exact normalization front end.  It does not glue the vertex
corners, compute the global normalized Rees algebra, or prove simultaneous
normalization in the weighted parameter.
"""

from __future__ import annotations

import itertools
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from verify_quartic_biprojective_graph import (  # noqa: E402
    X0,
    X1,
    X2,
    X3,
    quartic_projective_coordinates,
)


u, v, exceptional, slope = sp.symbols("u v exceptional slope")
t, alpha = sp.symbols("t alpha")
a, b, c = sp.symbols("a b c")

FUNCTION_FIELD = sp.QQ.frac_field(t)
ALGEBRAIC_COEFFICIENT_FIELD = sp.QQ.frac_field(t, alpha)


@dataclass(frozen=True)
class ClusterPoint:
    label: str
    residue_degree: int
    multiplicity: int


def surface_poly(expression: sp.Expr) -> sp.Poly:
    return sp.Poly(sp.cancel(expression), u, v, domain=FUNCTION_FIELD)


def surface_order(expression: sp.Expr) -> int:
    polynomial = surface_poly(expression)
    if polynomial.is_zero:
        return 10**9
    return min(sum(monomial) for monomial in polynomial.monoms())


def leading_form(expression: sp.Expr, order: int) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * u**monomial[0] * v**monomial[1]
            for monomial, coefficient in surface_poly(expression).terms()
            if sum(monomial) == order
        )
    )


def surface_data(
    generators: tuple[sp.Expr, ...],
) -> tuple[int, tuple[int, ...], sp.Poly, tuple[tuple[sp.Poly, int], ...]]:
    orders = tuple(surface_order(generator) for generator in generators)
    multiplicity = min(orders)
    minimal_leads = tuple(
        leading_form(generator, multiplicity)
        for generator, order in zip(generators, orders, strict=True)
        if order == multiplicity
    )
    common = surface_poly(minimal_leads[0])
    for form in minimal_leads[1:]:
        common = sp.gcd(common, surface_poly(form))
    factors = tuple(sp.factor_list(common)[1])
    return multiplicity, orders, common, factors


def assert_associate(actual: sp.Expr, expected: sp.Expr) -> None:
    quotient = sp.cancel(actual / expected)
    assert quotient != 0
    assert not ({u, v} & quotient.free_symbols), (actual, expected, quotient)


def assert_factorization(
    factors: tuple[tuple[sp.Poly, int], ...],
    expected: tuple[tuple[sp.Expr, int], ...],
) -> None:
    assert len(factors) == len(expected), (factors, expected)
    unmatched = list(factors)
    for expected_factor, expected_power in expected:
        for index, (factor, power) in enumerate(unmatched):
            if power != expected_power:
                continue
            try:
                assert_associate(factor.as_expr(), expected_factor)
            except AssertionError:
                continue
            unmatched.pop(index)
            break
        else:
            raise AssertionError((factors, expected))
    assert not unmatched


def linear_quadratic_transform(
    generators: tuple[sp.Expr, ...],
    multiplicity: int,
    direction: sp.Expr,
) -> tuple[sp.Expr, ...]:
    factor = surface_poly(direction)
    assert factor.total_degree() == 1
    u_coefficient = factor.coeff_monomial(u)
    v_coefficient = factor.coeff_monomial(v)
    if v_coefficient == 0:
        substitution = {u: exceptional * slope, v: exceptional}
    else:
        root = sp.cancel(-u_coefficient / v_coefficient)
        substitution = {
            u: exceptional,
            v: exceptional * (root + slope),
        }
    result = []
    for generator in generators:
        transformed = sp.cancel(
            sp.expand(generator.subs(substitution)) / exceptional**multiplicity
        ).subs({exceptional: u, slope: v})
        surface_poly(transformed)
        result.append(transformed)
    return tuple(result)


def algebraic_terminal_child(
    generators: tuple[sp.Expr, ...],
    multiplicity: int,
    direction: sp.Poly,
) -> int:
    """Verify that an irreducible closed direction has one terminal child."""

    minimal_polynomial = sp.Poly(
        sp.expand(direction.as_expr().subs({u: 1, v: alpha})),
        alpha,
        domain=FUNCTION_FIELD,
    ).monic()
    factor_degrees = tuple(
        (factor.degree(), power)
        for factor, power in sp.factor_list(minimal_polynomial)[1]
    )
    assert factor_degrees == ((minimal_polynomial.degree(), 1),)

    transformed = tuple(
        sp.cancel(
            sp.expand(
                generator.subs(
                    {u: exceptional, v: exceptional * (alpha + slope)}
                )
            )
            / exceptional**multiplicity
        ).subs({exceptional: u, slope: v})
        for generator in generators
    )

    def reduce_coefficient(coefficient: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.fraction(sp.cancel(coefficient))
        numerator_remainder = sp.Poly(
            numerator, alpha, domain=FUNCTION_FIELD
        ).rem(minimal_polynomial)
        denominator_remainder = sp.Poly(
            denominator, alpha, domain=FUNCTION_FIELD
        ).rem(minimal_polynomial)
        assert not denominator_remainder.is_zero
        return sp.cancel(
            numerator_remainder.as_expr() / denominator_remainder.as_expr()
        )

    orders: list[int] = []
    leads: list[sp.Expr] = []
    for generator in transformed:
        polynomial = sp.Poly(
            generator,
            u,
            v,
            domain=ALGEBRAIC_COEFFICIENT_FIELD,
        )
        reduced_terms = tuple(
            (monomial, reduce_coefficient(coefficient))
            for monomial, coefficient in polynomial.terms()
            if reduce_coefficient(coefficient) != 0
        )
        order = min(sum(monomial) for monomial, _ in reduced_terms)
        orders.append(order)
        leads.append(
            sp.expand(
                sum(
                    coefficient * u**monomial[0] * v**monomial[1]
                    for monomial, coefficient in reduced_terms
                    if sum(monomial) == order
                )
            )
        )

    child_multiplicity = min(orders)
    assert child_multiplicity == 1
    minimal_leads = tuple(
        lead for lead, order in zip(leads, orders, strict=True) if order == 1
    )
    assert len(minimal_leads) >= 2

    rows = []
    for form in minimal_leads:
        polynomial = sp.Poly(
            form, u, v, domain=ALGEBRAIC_COEFFICIENT_FIELD
        )
        rows.append(
            (polynomial.coeff_monomial(u), polynomial.coeff_monomial(v))
        )
    determinants = tuple(
        reduce_coefficient(
            rows[first][0] * rows[second][1]
            - rows[first][1] * rows[second][0]
        )
        for first, second in itertools.combinations(range(len(rows)), 2)
    )
    assert any(determinant != 0 for determinant in determinants)
    return minimal_polynomial.degree()


def monomial_colength(leading_exponents: tuple[tuple[int, int], ...]) -> int:
    bound = 1 + max(max(exponent) for exponent in leading_exponents)
    return sum(
        1
        for first in range(bound)
        for second in range(bound)
        if not any(
            first >= generator[0] and second >= generator[1]
            for generator in leading_exponents
        )
    )


def original_surface_colength(generators: tuple[sp.Expr, ...]) -> int:
    basis = sp.groebner(
        generators,
        u,
        v,
        domain=FUNCTION_FIELD,
        order="grevlex",
    )
    leading_exponents = tuple(
        polynomial.LM(order=basis.order).exponents for polynomial in basis.polys
    )
    return monomial_colength(leading_exponents)


def cluster_colength(points: tuple[ClusterPoint, ...]) -> int:
    return sum(
        point.residue_degree
        * point.multiplicity
        * (point.multiplicity + 1)
        // 2
        for point in points
    )


def verify_generic_lines(
    coordinates: tuple[sp.Expr, ...],
) -> dict[str, tuple[int, int, int]]:
    line_substitutions = {
        "L1": {X0: u, X1: v, X2: t, X3: 1},
        "L2": {X0: u, X1: t, X2: v, X3: 1},
        "L3": {X0: u, X1: t, X2: 1, X3: v},
    }
    ideals = {
        name: tuple(sp.expand(item.subs(substitution)) for item in coordinates)
        for name, substitution in line_substitutions.items()
    }

    # L1: one rational and one cubic child above the first direction, and a
    # quadratic child after the second two-point chain.
    l1 = ideals["L1"]
    multiplicity, orders, _, factors = surface_data(l1)
    assert (multiplicity, orders) == (6, (12, 6, 7, 11))
    assert_factorization(factors, ((v, 4), (5 * t * u - 3 * v, 2)))

    l1_first = linear_quadratic_transform(l1, 6, v)
    first_multiplicity, first_orders, _, first_factors = surface_data(l1_first)
    assert (first_multiplicity, first_orders) == (4, (6, 4, 5, 7))
    assert len(first_factors) == 2
    linear_factor = next(
        factor for factor, _ in first_factors if factor.total_degree() == 1
    )
    cubic_factor = next(
        factor for factor, _ in first_factors if factor.total_degree() == 3
    )
    assert_associate(linear_factor.as_expr(), t * v + u)
    terminal_multiplicity, _, terminal_gcd, terminal_factors = surface_data(
        linear_quadratic_transform(l1_first, 4, linear_factor.as_expr())
    )
    assert terminal_multiplicity == 1
    assert terminal_gcd.total_degree() == 0 and not terminal_factors
    assert algebraic_terminal_child(l1_first, 4, cubic_factor) == 3

    l1_second = linear_quadratic_transform(l1, 6, 5 * t * u - 3 * v)
    second_multiplicity, second_orders, _, second_factors = surface_data(l1_second)
    assert (second_multiplicity, second_orders) == (2, (6, 2, 3, 6))
    assert_factorization(second_factors, ((5 * t * v + 3 * u, 2),))
    l1_second_child = linear_quadratic_transform(
        l1_second, 2, second_factors[0][0].as_expr()
    )
    second_child_multiplicity, second_child_orders, _, second_child_factors = (
        surface_data(l1_second_child)
    )
    assert (second_child_multiplicity, second_child_orders) == (
        2,
        (4, 2, 3, 5),
    )
    assert len(second_child_factors) == 1
    assert second_child_factors[0][0].total_degree() == 2
    assert algebraic_terminal_child(
        l1_second_child, 2, second_child_factors[0][0]
    ) == 2

    l1_points = (
        ClusterPoint("root", 1, 6),
        ClusterPoint("first direction", 1, 4),
        ClusterPoint("first rational child", 1, 1),
        ClusterPoint("first cubic child", 3, 1),
        ClusterPoint("second direction", 1, 2),
        ClusterPoint("second rational child", 1, 2),
        ClusterPoint("second quadratic child", 2, 1),
    )

    # L2 is a single rational chain with multiplicities 4,3,1,1,1.
    l2 = ideals["L2"]
    l2_multiplicities = []
    l2_directions = (v, t * v + u, u, v)
    current = l2
    for direction in l2_directions:
        current_multiplicity, _, _, current_factors = surface_data(current)
        l2_multiplicities.append(current_multiplicity)
        assert len(current_factors) == 1
        assert_associate(current_factors[0][0].as_expr(), direction)
        current = linear_quadratic_transform(
            current, current_multiplicity, direction
        )
    final_multiplicity, _, final_gcd, final_factors = surface_data(current)
    l2_multiplicities.append(final_multiplicity)
    assert tuple(l2_multiplicities) == (4, 3, 1, 1, 1)
    assert final_gcd.total_degree() == 0 and not final_factors
    l2_points = tuple(
        ClusterPoint(f"chain {index}", 1, multiplicity)
        for index, multiplicity in enumerate(l2_multiplicities)
    )

    # L3 is a five-point rational multiplicity-two chain followed by one
    # terminal quadratic closed point.
    l3 = ideals["L3"]
    l3_directions = (
        -3 * t * v + 5 * u,
        v,
        t**2 * v + u,
        v,
    )
    current = l3
    l3_multiplicities = []
    for direction in l3_directions:
        current_multiplicity, _, _, current_factors = surface_data(current)
        l3_multiplicities.append(current_multiplicity)
        assert current_multiplicity == 2
        assert len(current_factors) == 1
        assert_associate(current_factors[0][0].as_expr(), direction)
        current = linear_quadratic_transform(
            current, current_multiplicity, direction
        )
    last_multiplicity, _, _, last_factors = surface_data(current)
    l3_multiplicities.append(last_multiplicity)
    assert tuple(l3_multiplicities) == (2, 2, 2, 2, 2)
    assert len(last_factors) == 1 and last_factors[0][0].total_degree() == 2
    assert algebraic_terminal_child(current, 2, last_factors[0][0]) == 2
    l3_points = tuple(
        ClusterPoint(f"chain {index}", 1, multiplicity)
        for index, multiplicity in enumerate(l3_multiplicities)
    ) + (ClusterPoint("quadratic child", 2, 1),)

    points = {"L1": l1_points, "L2": l2_points, "L3": l3_points}
    expected = {
        "L1": (52, 43, 9),
        "L2": (25, 19, 6),
        "L3": (18, 17, 1),
    }
    result = {}
    for name in ("L1", "L2", "L3"):
        original = original_surface_colength(ideals[name])
        normalized = cluster_colength(points[name])
        row = (original, normalized, original - normalized)
        assert row == expected[name]
        result[name] = row
    return result


def minimal_support_points(generators: tuple[sp.Expr, ...]) -> tuple[tuple[int, ...], ...]:
    support = {
        monomial
        for generator in generators
        for monomial, _ in sp.Poly(generator, a, b, c).terms()
    }
    return tuple(
        sorted(
            point
            for point in support
            if not any(
                candidate != point
                and all(candidate[index] <= point[index] for index in range(3))
                for candidate in support
            )
        )
    )


def compact_newton_facets(
    generators: tuple[sp.Expr, ...],
) -> dict[tuple[int, int, int], int]:
    points = minimal_support_points(generators)
    facets: dict[tuple[int, int, int], int] = {}
    for first, second, third in itertools.combinations(points, 3):
        vector_first = sp.Matrix(
            [second[index] - first[index] for index in range(3)]
        )
        vector_second = sp.Matrix(
            [third[index] - first[index] for index in range(3)]
        )
        normal_vector = vector_first.cross(vector_second)
        if normal_vector == sp.zeros(3, 1):
            continue
        normal = tuple(int(entry) for entry in normal_vector)
        constant = sum(normal[index] * first[index] for index in range(3))
        values = tuple(
            sum(normal[index] * point[index] for index in range(3))
            for point in points
        )
        if all(value <= constant for value in values):
            normal = tuple(-entry for entry in normal)
            constant = -constant
            values = tuple(-value for value in values)
        elif not all(value >= constant for value in values):
            continue
        if not all(entry > 0 for entry in normal):
            continue
        divisor = math.gcd(math.gcd(normal[0], normal[1]), normal[2])
        primitive = tuple(entry // divisor for entry in normal)
        primitive_constant = constant // divisor
        face = tuple(
            point
            for point in points
            if sum(primitive[index] * point[index] for index in range(3))
            == primitive_constant
        )
        affine_rank = sp.Matrix(
            [
                [point[index] - face[0][index] for index in range(3)]
                for point in face[1:]
            ]
        ).rank()
        if affine_rank == 2:
            facets[primitive] = primitive_constant
    return facets


def weight_order(expression: sp.Expr, weight: tuple[int, int, int]) -> int:
    return min(
        sum(weight[index] * monomial[index] for index in range(3))
        for monomial, _ in sp.Poly(expression, a, b, c).terms()
    )


def weight_initial(
    expression: sp.Expr,
    weight: tuple[int, int, int],
    order: int,
) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * a**monomial[0] * b**monomial[1] * c**monomial[2]
            for monomial, coefficient in sp.Poly(expression, a, b, c).terms()
            if sum(weight[index] * monomial[index] for index in range(3))
            == order
        )
    )


def torus_colors(initials: tuple[sp.Expr, ...]) -> tuple[tuple[sp.Expr, int], ...]:
    common = sp.Poly(initials[0], a, b, c, domain=sp.QQ)
    for initial in initials[1:]:
        common = sp.gcd(common, sp.Poly(initial, a, b, c, domain=sp.QQ))
    colors = []
    for factor, power in sp.factor_list(common)[1]:
        if len(factor.terms()) == 1:
            continue
        colors.append((sp.factor(factor.as_expr()), power))
    return tuple(colors)


def rational_color_cluster(
    coordinates: tuple[sp.Expr, ...],
    substitution: dict[sp.Symbol, sp.Expr],
    exceptional_order: int,
    directions: tuple[sp.Expr, ...],
    expected_multiplicities: tuple[int, ...],
) -> None:
    current = tuple(
        sp.cancel(sp.expand(item.subs(substitution)) / u**exceptional_order)
        for item in coordinates
    )
    multiplicities = []
    for direction in directions:
        multiplicity, _, _, factors = surface_data(current)
        multiplicities.append(multiplicity)
        assert len(factors) == 1
        assert_associate(factors[0][0].as_expr(), direction)
        current = linear_quadratic_transform(current, multiplicity, direction)
    final_multiplicity, _, final_gcd, final_factors = surface_data(current)
    multiplicities.append(final_multiplicity)
    assert tuple(multiplicities) == expected_multiplicities
    assert final_gcd.total_degree() == 0 and not final_factors


def verify_degree_five_color(coordinates: tuple[sp.Expr, ...]) -> None:
    local = tuple(
        sp.expand(item.subs({X0: a, X1: b, X2: c, X3: 1}))
        for item in coordinates
    )
    weight = (2, 3, 1)
    order = weight_order(local[1], weight)
    initial = weight_initial(local[1], weight, order)
    quadratic_color = a**2 + b * c
    degree_five_color = sp.Poly(
        sp.cancel(initial / quadratic_color).subs({a: t, b: alpha, c: 1}),
        alpha,
        domain=FUNCTION_FIELD,
    )
    assert tuple(
        (factor.degree(), power)
        for factor, power in sp.factor_list(degree_five_color)[1]
    ) == ((5, 1),)
    assert sp.factor(
        sp.cancel(initial / quadratic_color).subs({a: t, b: -t**2, c: 1})
    ) == t**8

    substitution = {
        X0: u**2 * t,
        X1: u**3 * (alpha + v),
        X2: u,
        X3: 1,
    }
    generators = tuple(
        sp.cancel(sp.expand(item.subs(substitution)) / u**22)
        for item in coordinates
    )

    minimal_polynomial = degree_five_color.monic()

    def reduce_coefficient(coefficient: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.fraction(sp.cancel(coefficient))
        numerator_remainder = sp.Poly(
            numerator, alpha, domain=FUNCTION_FIELD
        ).rem(minimal_polynomial)
        denominator_remainder = sp.Poly(
            denominator, alpha, domain=FUNCTION_FIELD
        ).rem(minimal_polynomial)
        assert not denominator_remainder.is_zero
        return sp.cancel(
            numerator_remainder.as_expr() / denominator_remainder.as_expr()
        )

    orders = []
    leads = []
    for generator in generators:
        polynomial = sp.Poly(
            generator,
            u,
            v,
            domain=ALGEBRAIC_COEFFICIENT_FIELD,
        )
        terms = tuple(
            (monomial, reduce_coefficient(coefficient))
            for monomial, coefficient in polynomial.terms()
            if reduce_coefficient(coefficient) != 0
        )
        generator_order = min(sum(monomial) for monomial, _ in terms)
        orders.append(generator_order)
        leads.append(
            sp.expand(
                sum(
                    coefficient * u**monomial[0] * v**monomial[1]
                    for monomial, coefficient in terms
                    if sum(monomial) == generator_order
                )
            )
        )
    assert min(orders) == 1
    minimal_leads = tuple(
        lead for lead, generator_order in zip(leads, orders, strict=True)
        if generator_order == 1
    )
    assert len(minimal_leads) == 2
    rows = []
    for lead in minimal_leads:
        polynomial = sp.Poly(
            lead, u, v, domain=ALGEBRAIC_COEFFICIENT_FIELD
        )
        rows.append(
            (polynomial.coeff_monomial(u), polynomial.coeff_monomial(v))
        )
    determinant = reduce_coefficient(
        rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
    )
    assert determinant != 0


def verify_vertex_front_end(coordinates: tuple[sp.Expr, ...]) -> None:
    vertices = {
        "V12": {X0: a, X1: b, X2: c, X3: 1},
        "V13": {X0: a, X1: b, X2: 1, X3: c},
        "V23": {X0: a, X1: 1, X2: b, X3: c},
    }
    local_ideals = {
        name: tuple(sp.expand(item.subs(substitution)) for item in coordinates)
        for name, substitution in vertices.items()
    }
    expected_facets = {
        "V12": {(1, 1, 1): 10, (2, 3, 1): 22},
        "V13": {},
        "V23": {
            (1, 1, 2): 8,
            (1, 2, 1): 9,
            (1, 2, 3): 11,
            (2, 3, 5): 21,
        },
    }
    expected_orders = {
        ("V12", (1, 1, 1)): (12, 10, 10, 11),
        ("V12", (2, 3, 1)): (24, 22, 23, 25),
        ("V23", (1, 1, 2)): (12, 8, 8, 10),
        ("V23", (1, 2, 1)): (12, 10, 9, 9),
        ("V23", (1, 2, 3)): (12, 14, 13, 11),
        ("V23", (2, 3, 5)): (24, 22, 21, 21),
    }
    expected_color_degrees = {
        ("V12", (1, 1, 1)): ((2, 3),),
        ("V12", (2, 3, 1)): ((4, 1), (18, 1)),
        ("V23", (1, 1, 2)): ((2, 2),),
        ("V23", (1, 2, 1)): (),
        ("V23", (1, 2, 3)): ((3, 1),),
        ("V23", (2, 3, 5)): ((5, 1),),
    }

    for name, generators in local_ideals.items():
        facets = compact_newton_facets(generators)
        assert facets == expected_facets[name], (name, facets)
        for weight, ideal_order in facets.items():
            coordinate_orders = tuple(
                weight_order(generator, weight) for generator in generators
            )
            assert coordinate_orders == expected_orders[(name, weight)]
            initials = tuple(
                weight_initial(generator, weight, coordinate_order)
                for generator, coordinate_order in zip(
                    generators, coordinate_orders, strict=True
                )
                if coordinate_order == ideal_order
            )
            colors = torus_colors(initials)
            signature = tuple(
                (weight_order(color, weight), power)
                for color, power in colors
            )
            assert signature == expected_color_degrees[(name, weight)], (
                name,
                weight,
                colors,
            )

    # Principalize the five rational generic colors.
    rational_color_cluster(
        coordinates,
        {X0: u * t, X1: u, X2: u * (-t**2 + v), X3: 1},
        10,
        (u, v),
        (1, 1, 1),
    )
    rational_color_cluster(
        coordinates,
        {X0: u**2 * t, X1: u**3, X2: u * (-t**2 + v), X3: 1},
        22,
        (),
        (1,),
    )
    rational_color_cluster(
        coordinates,
        {X0: u, X1: 1, X2: u * t, X3: u**2 * (sp.Rational(5, 3) * t + v)},
        8,
        (u + v,),
        (2, 2),
    )
    rational_color_cluster(
        coordinates,
        {X0: u, X1: 1, X2: u**2 * t, X3: u**3 * ((5 * t - 3) / 3 + v)},
        11,
        (),
        (1,),
    )
    rational_color_cluster(
        coordinates,
        {X0: u**2, X1: 1, X2: u**3 * t, X3: u**5 * (sp.Rational(5, 3) * t + v)},
        21,
        (u + v, v),
        (1, 1, 1),
    )
    verify_degree_five_color(coordinates)


def main() -> None:
    coordinates = quartic_projective_coordinates()
    line_rows = verify_generic_lines(coordinates)
    verify_vertex_front_end(coordinates)
    summary = ", ".join(
        f"{name}={row}" for name, row in line_rows.items()
    )
    print(
        "PASS quartic Rees stratification: "
        f"(original, normalized, defect) {summary}; "
        "six compact vertex facets and six generic colors principalized"
    )


if __name__ == "__main__":
    main()

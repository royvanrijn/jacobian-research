#!/usr/bin/env python3
"""Exact checks for universal all-degree power-shifted gauge multiplicity."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.keller_fiber import (
    compile_polynomial_to_keller_fiber,
    quadratic_gauge_map,
    stable_multiplicity_gauge_map,
)


x, y, z, S = sp.symbols("x y z S")
g1, g2, g3, g4, g5 = sp.symbols("g1 g2 g3 g4 g5", nonzero=True)

t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
P = t * q
root = x / t
Q = y + x * q
D = sp.factor(1 - root * (Q - P * root))
assert D == 1 / t


def shifted_map(coefficients: list[sp.Expr], shift: int) -> tuple[sp.Expr, ...]:
    second = y + 3 * (g3 / g1) * x * q + 2 * (g2 / g1) * t * q
    third = x * (5 - 3 * t) - (g3 / g1) * x**3 * z
    for degree, coefficient in enumerate(coefficients, start=4):
        second += (
            degree
            * (coefficient / g1)
            * t ** (shift + 2)
            * x ** (degree - 2)
            * q ** (degree + shift)
        )
        third -= (
            (degree - 2)
            * (coefficient / g1)
            * t**shift
            * x**degree
            * q ** (degree + shift)
        )
    return P, sp.expand(second), sp.expand(third)


def lifted_seed(
    coefficients: list[sp.Expr], shift: int, variable: sp.Expr
) -> sp.Expr:
    result = g1 * variable + P * (g2 * variable**2 + g3 * variable**3)
    for degree, coefficient in enumerate(coefficients, start=4):
        result += coefficient * P ** (degree + shift) * variable**degree
    return result


# Direct three-variable checks in the first two degrees.  The marked-line
# proof is uniform; these guard both a single decoration and a nontrivial sum.
for coefficients, shift in (([g4], 2), ([g4, g5], 1)):
    first, second, third = shifted_map(coefficients, shift)
    inverse = sp.factor(
        lifted_seed(coefficients, shift, root)
        - g1 * (second * root**2 + third) / 2
    )
    assert inverse == 0

    derivative = sp.diff(
        lifted_seed(coefficients, shift, S)
        - g1 * (second * S**2 + third) / 2,
        S,
    ).subs(S, root)
    assert sp.factor(derivative - g1 * D) == 0

    jacobian = sp.det(
        sp.Matrix(
            [
                [sp.diff(component, variable) for variable in (x, y, z)]
                for component in (first, second, third)
            ]
        )
    )
    assert sp.factor(jacobian) == -2


def convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    points = sorted(set(points))

    def cross(
        origin: tuple[int, int],
        first: tuple[int, int],
        second: tuple[int, int],
    ) -> int:
        return (
            (first[0] - origin[0]) * (second[1] - origin[1])
            - (first[1] - origin[1]) * (second[0] - origin[0])
        )

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)

    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def normalized_area(polygon: list[tuple[int, int]]) -> int:
    return abs(
        sum(
            polygon[index][0] * polygon[(index + 1) % len(polygon)][1]
            - polygon[index][1] * polygon[(index + 1) % len(polygon)][0]
            for index in range(len(polygon))
        )
    )


# Fitting support, hull vertices, and normalized area.
for degree in range(4, 13):
    for shift in range(7):
        support = [(0, 0), (1, 2)] + [
            (index + shift, index - 1) for index in range(4, degree + 1)
        ]
        hull = convex_hull(support)
        expected_vertices = (
            [(0, 0), (4 + shift, 3), (1, 2)]
            if degree == 4
            else [
                (0, 0),
                (4 + shift, 3),
                (degree + shift, degree - 1),
                (1, 2),
            ]
        )
        assert hull == expected_vertices
        assert normalized_area(hull) == 2 * degree - 3 + (degree - 2) * shift

        # P=0 lower Newton polygon: affine blocks of lengths two and one,
        # followed by one block of length N-3 and height N+m-1.
        vertices = [(0, 0), (2, 0), (3, 1), (degree, degree + shift)]
        slopes = [
            sp.Rational(y2 - y1, x2 - x1)
            for (x1, y1), (x2, y2) in zip(vertices, vertices[1:])
        ]
        assert slopes == [
            0,
            1,
            sp.Rational(degree + shift - 1, degree - 3),
        ]


# Symbolic shoelace formula for the Newton quadrilateral.
N, m = sp.symbols("N m", integer=True, positive=True)
vertices_symbolic = [(0, 0), (4 + m, 3), (N + m, N - 1), (1, 2)]
area_symbolic = sp.expand(
    sum(
        vertices_symbolic[index][0]
        * vertices_symbolic[(index + 1) % 4][1]
        - vertices_symbolic[index][1]
        * vertices_symbolic[(index + 1) % 4][0]
        for index in range(4)
    )
)
assert sp.expand(area_symbolic - (2 * N - 3 + (N - 2) * m)) == 0


# The public compiler realizes the symbolic family and returns its exact
# Newton-area and P=0 boundary certificate.
T = sp.symbols("T")
quartic = T**4 + T**3 + T + 2
minimal_quartic = compile_polynomial_to_keller_fiber(
    quartic,
    T,
    translation=0,
    inverse_variable=S,
    source_variables=(x, y, z),
)
shift_zero_quartic = compile_polynomial_to_keller_fiber(
    quartic,
    T,
    translation=0,
    inverse_variable=S,
    source_variables=(x, y, z),
    stable_parameter=0,
)
shift_three_quartic = compile_polynomial_to_keller_fiber(
    quartic,
    T,
    translation=0,
    inverse_variable=S,
    source_variables=(x, y, z),
    stable_parameter=3,
)
assert minimal_quartic.determinant_minus_two_map == (
    shift_zero_quartic.determinant_minus_two_map
)
assert minimal_quartic.determinant_minus_two_map == quadratic_gauge_map(
    minimal_quartic.seed, S, (x, y, z)
)
expected_shift_three = tuple(
    sp.expand(component.subs({g1: 1, g2: 0, g3: 1, g4: 1}))
    for component in shifted_map([g4], 3)
)
assert all(
    sp.expand(actual - expected) == 0
    for actual, expected in zip(
        shift_three_quartic.determinant_minus_two_map,
        expected_shift_three,
    )
)
quartic_certificate = shift_three_quartic.stable_multiplicity
assert quartic_certificate is not None
assert quartic_certificate.separation_invariant == (
    "normalized_fitting_newton_area"
)
assert quartic_certificate.separation_value == 11
assert quartic_certificate.fitting_support == (
    (0, 0),
    (1, 2),
    (7, 3),
)
assert quartic_certificate.boundary_prime_count == 1
assert quartic_certificate.boundary_ramification_index == 1
assert shift_three_quartic.inverse_polynomial == quartic.subs(T, S)
assert sp.expand(
    shift_three_quartic.lifted_seed.subs({x: 0, y: 0, z: 1})
    - shift_three_quartic.seed
) == 0

# A stable higher-degree compilation automatically avoids every vanishing
# coefficient g_4,...,g_N.  The ordinary compiler keeps its weaker historical
# translation rule.
quintic = T**5 + T**3 + T + 1
ordinary_quintic = compile_polynomial_to_keller_fiber(
    quintic,
    T,
    inverse_variable=S,
    source_variables=(x, y, z),
)
stable_quintic = compile_polynomial_to_keller_fiber(
    quintic,
    T,
    inverse_variable=S,
    source_variables=(x, y, z),
    stable_parameter=2,
)
assert ordinary_quintic.translation == 0
assert stable_quintic.translation == 1
for derivative_order in (1, 3, 4, 5):
    assert sp.diff(quintic, T, derivative_order).subs(
        T, stable_quintic.translation
    )
quintic_certificate = stable_quintic.stable_multiplicity
assert quintic_certificate is not None
assert quintic_certificate.separation_value == 13
assert quintic_certificate.fitting_support == (
    (0, 0),
    (1, 2),
    (6, 3),
    (7, 4),
)
assert quintic_certificate.boundary_prime_count == 2
assert quintic_certificate.boundary_ramification_index == 1
try:
    compile_polynomial_to_keller_fiber(
        quintic,
        T,
        translation=0,
        stable_parameter=2,
    )
except ValueError as error:
    assert "orders 1,3,4,5" in str(error)
else:
    raise AssertionError("stable compiler accepted a vanishing higher coefficient")
try:
    stable_multiplicity_gauge_map(
        S**5 + S**3 + S,
        S,
        (x, y, z),
        1,
    )
except ValueError as error:
    assert "g_4,...,g_N" in str(error)
else:
    raise AssertionError("raw stable gauge accepted a missing higher coefficient")
for invalid_parameter in (-1, sp.Rational(1, 2), True):
    try:
        compile_polynomial_to_keller_fiber(
            quartic,
            T,
            translation=0,
            stable_parameter=invalid_parameter,
        )
    except ValueError as error:
        assert "nonnegative integer" in str(error)
    else:
        raise AssertionError(f"accepted invalid stable parameter {invalid_parameter!r}")


print("PASS: common power shifts preserve the determinant-minus-two identity")
print("PASS: the inverse and derivative reconstruction identities are exact")
print("PASS: the selected inverse polynomial is independent of the shift")
print("PASS: the P=0 Newton ledger has the stated three blocks")
print("PASS: Fitting Newton area is 2*N-3+(N-2)*m")
print("PASS: the public compiler returns exact power-shift certificates")

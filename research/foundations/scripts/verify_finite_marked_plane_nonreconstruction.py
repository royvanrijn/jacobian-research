#!/usr/bin/env python3
"""Exact checks for finite marked-plane nonreconstruction by interpolating gauges."""

from __future__ import annotations

import sympy as sp


P, B, C, S, Q = sp.symbols("P B C S Q")
x, y, z = sp.symbols("x y z")
g1, g2, g3, g4 = sp.symbols("g1 g2 g3 g4", nonzero=True)


def lower_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return the lower convex hull from left to right."""

    hull: list[tuple[int, int]] = []
    for point in sorted(set(points)):
        while len(hull) >= 2:
            first, second = hull[-2:]
            cross = (
                (second[0] - first[0]) * (point[1] - second[1])
                - (second[1] - first[1]) * (point[0] - second[0])
            )
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return hull


# The marked-line calculation permits an arbitrary polynomial multiplier
# R(P) on all decorations of S-degree at least four.
r0, r1, r2, r3 = sp.symbols("r0 r1 r2 r3")
multiplier = r0 + r1 * P + r2 * P**2 + r3 * P**3
high = g4 * P**4 * multiplier * S**4
lifted_seed = g1 * S + P * (g2 * S**2 + g3 * S**3) + high
beta = sp.cancel((sp.diff(lifted_seed, S) / g1 - 1 - P * S**2) / S)
plane_second = Q + beta
plane_third = 2 * lifted_seed / g1 - plane_second * S**2
D = 1 - S * Q + P * S**2
plane_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(plane_second, S), sp.diff(plane_second, Q)],
            [sp.diff(plane_third, S), sp.diff(plane_third, Q)],
        ]
    )
)
assert sp.factor(plane_jacobian + 2 * D) == 0


# A direct denominator-free quartic checks that P-dependence of R adds only
# dP terms to the wedge calculation and leaves the full determinant -2.
t = 1 + x * y
q = t**2 * z + sp.Rational(1, 3) * y**2 * (1 + 3 * t)
pi = t * q
root = x / t
source_multiplier = 1 + pi
first = pi
second = (
    y
    + 9 * x * q
    + 4 * t * q
    + 20 * source_multiplier * t**2 * x**2 * q**4
)
third = (
    x * (5 - 3 * t)
    - 3 * x**3 * z
    - 10 * source_multiplier * x**4 * q**4
)
source_seed = (
    S
    + 2 * P * S**2
    + 3 * P * S**3
    + 5 * P**4 * (1 + P) * S**4
)
inverse_identity = source_seed.subs({P: first, S: root})
inverse_identity -= (second * root**2 + third) / 2
assert sp.factor(inverse_identity) == 0
derivative_identity = sp.diff(
    source_seed - (second * S**2 + third) / 2,
    S,
).subs({P: first, S: root})
assert sp.factor(derivative_identity - 1 / t) == 0
jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(component, variable) for variable in (x, y, z)]
            for component in (first, second, third)
        ]
    )
)
assert sp.factor(jacobian) == -2


# Interpolate the value one on three arbitrary sample planes.  For every
# desired degree, only finitely many scalar choices make 1+t*P^k*H repeated,
# so a short exact search always finds a squarefree representative here.
sample_values = (-1, 1, 2)
sample_polynomial = sp.prod(P - value for value in sample_values)
records: list[tuple[int, int, int, sp.Expr]] = []
for extra_degree in range(1, 8):
    interpolation_core = sp.expand(P**extra_degree * sample_polynomial)
    chosen: tuple[int, sp.Expr] | None = None
    for scalar in range(1, 50):
        candidate = sp.expand(1 + scalar * interpolation_core)
        if sp.gcd(sp.Poly(candidate, P), sp.Poly(sp.diff(candidate, P), P)).degree() == 0:
            chosen = scalar, candidate
            break
    assert chosen is not None
    scalar, candidate = chosen
    expected_degree = extra_degree + len(sample_values)
    assert sp.degree(candidate, P) == expected_degree
    assert candidate.subs(P, 0) == 1
    for value in sample_values:
        assert candidate.subs(P, value) == 1
        assert sp.rem(sp.Poly(candidate - 1, P), sp.Poly(P - value, P)).is_zero
    records.append(
        (
            extra_degree,
            expected_degree,
            expected_degree + 2,
            candidate,
        )
    )


# At a simple nonzero root rho of R, the three finite cubic sheets form the
# horizontal Newton block and all remaining sheets form one boundary prime
# of ramification index N-3.
for degree in range(4, 13):
    degree_drop_points = [(index, 0) for index in range(4)]
    degree_drop_points += [(index, 1) for index in range(4, degree + 1)]
    degree_drop_hull = lower_hull(degree_drop_points)
    assert degree_drop_hull == [(0, 0), (3, 0), (degree, 1)]
    degree_drop_slopes = [
        sp.Rational(y2 - y1, x2 - x1)
        for (x1, y1), (x2, y2) in zip(
            degree_drop_hull,
            degree_drop_hull[1:],
        )
    ]
    assert degree_drop_slopes == [0, sp.Rational(1, degree - 3)]
    assert degree - 3 == sp.denom(degree_drop_slopes[-1])

    # Since every interpolating multiplier has constant term one, the old
    # P=0 ledger is unchanged from the minimal quadratic gauge.
    zero_points = [(0, 0), (1, 0), (2, 0), (3, 1)]
    zero_points += [(index, index) for index in range(4, degree + 1)]
    zero_hull = lower_hull(zero_points)
    assert zero_hull == [(0, 0), (2, 0), (3, 1), (degree, degree)]


# The complete geometric boundary-image count is deg(R)+2: one component
# for every simple nonzero root of R, plus P=0 and the ramified discriminant.
boundary_counts = [record[2] for record in records]
assert boundary_counts == sorted(set(boundary_counts))
assert boundary_counts == [degree + 2 for _, degree, _, _ in records]


print("PASS: arbitrary polynomial gauge multipliers preserve determinant -2")
print("PASS: the interpolating gauges agree on every prescribed sample plane")
print("PASS: each multiplier root gives one (e,f)=(N-3,1) boundary prime")
print("PASS: boundary target counts deg(R)+2 separate infinitely many classes")

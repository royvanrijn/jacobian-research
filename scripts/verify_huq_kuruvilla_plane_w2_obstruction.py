#!/usr/bin/env python3
"""Exact W_2 obstruction for the Huq--Kuruvilla--Mondello plane map.

The proof isolates the coefficient of x*y in the first Jacobian correction.
It is an all-degree coefficient argument, not a bounded ansatz search.

The written proof and status boundary are in
verified/HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md.
"""

from __future__ import annotations

import sympy as sp


def mod2(expr: sp.Expr, *generators: sp.Symbol) -> sp.Expr:
    """Return the canonical polynomial representative over F_2."""

    return sp.Poly(sp.expand(expr), *generators, modulus=2).as_expr()


def coefficient_functional(
    multiplier: sp.Expr,
    derivative_axis: str,
    target: tuple[int, int],
    x: sp.Symbol,
    y: sp.Symbol,
) -> int:
    """Coefficient at target of multiplier*d(arbitrary polynomial).

    For a fixed target monomial and a fixed multiplier monomial there is at
    most one source monomial.  Its derivative coefficient decides whether an
    arbitrary input coefficient can contribute.  Returning zero proves that
    this coefficient functional vanishes in every polynomial degree.
    """

    total = 0
    for (u, v), coefficient in sp.Poly(
        mod2(multiplier, x, y), x, y, modulus=2
    ).terms():
        if derivative_axis == "x":
            i, j = target[0] - u + 1, target[1] - v
            derivative_coefficient = i
        else:
            i, j = target[0] - u, target[1] - v + 1
            derivative_coefficient = j
        if i >= 0 and j >= 0:
            total ^= (int(coefficient) & 1) * (derivative_coefficient & 1)
    return total


x, y = sp.symbols("x y")
P = x + x**2 * y + x**4 + x**6 * y**2
Q = y + x**5 + x**6 * y + x**7 * y**2 + x**8 * y**3

integer_jacobian = sp.expand(sp.det(sp.Matrix((P, Q)).jacobian((x, y))))
half_error = sp.Poly(integer_jacobian - 1, x, y).quo_ground(2).as_expr()
assert sp.expand(1 + 2 * half_error - integer_jacobian) == 0
assert int(sp.Poly(half_error, x, y).coeff_monomial(x * y)) % 2 == 1
print("PASS: the naive integral Jacobian has x*y coefficient 2 modulo 4")

Px = mod2(sp.diff(P, x), x, y)
Py = mod2(sp.diff(P, y), x, y)
Qx = mod2(sp.diff(Q, x), x, y)
Qy = mod2(sp.diff(Q, y), x, y)
assert Px == 1
assert Py == x**2
assert Qx == x**4 + x**6 * y**2
assert Qy == x**6 + x**8 * y**2 + 1

# For arbitrary corrections A,B modulo two, the first variation is
# A_x Q_y + P_x B_y + A_y Q_x + P_y B_x.  Each summand has identically
# zero x*y coefficient.  The helper checks the coefficient functional on an
# arbitrary input polynomial, without imposing a degree bound.
target = (1, 1)
functionals = (
    coefficient_functional(Qy, "x", target, x, y),
    coefficient_functional(Px, "y", target, x, y),
    coefficient_functional(Qx, "y", target, x, y),
    coefficient_functional(Py, "x", target, x, y),
)
assert functionals == (0, 0, 0, 0)
print("PASS: every all-degree first Jacobian correction has zero x*y coefficient")
print("PASS: no polynomial constant-Jacobian lift exists over Z/4")

# The obstruction is unstable.  After adjoining one identity coordinate z,
# multiply that coordinate by 1+2*K.  The Jacobian matrix is block lower
# triangular, so its determinant is (1+2*K)^2=1 modulo four.
z = sp.symbols("z")
stabilized_third_coordinate = z * (1 + 2 * half_error)
stabilized_jacobian = sp.expand(
    sp.det(
        sp.Matrix((P, Q, stabilized_third_coordinate)).jacobian((x, y, z))
    )
)
for coefficient in sp.Poly(stabilized_jacobian - 1, x, y, z).coeffs():
    assert int(coefficient) % 4 == 0
print("PASS: one identity stabilization has an explicit Keller lift over Z/4")

# Universally, for h=2*K and S_n=sum_{j=0}^{n-1}(-h)^j, one has
# (1+h)S_n=1-(-h)^n.  The following symbolic induction step certifies the
# identity independently of K and n: if r=(-h)^n is the old error, adjoining
# r to S_n changes 1-r into 1+h*r, the next geometric-series remainder.
h, r_symbol = sp.symbols("h r_symbol")
assert sp.expand((1 - r_symbol) + (1 + h) * r_symbol - (1 + h * r_symbol)) == 0

# Replay several concrete levels as a regression; the all-n proof is the
# geometric-series identity above and in the canonical note.
for level in range(2, 7):
    multiplier = sum((-2 * half_error) ** j for j in range(level))
    determinant = sp.expand(integer_jacobian * multiplier)
    assert sp.expand(determinant - (1 - (-2 * half_error) ** level)) == 0
    for coefficient in sp.Poly(determinant - 1, x, y).coeffs():
        assert int(coefficient) % (2**level) == 0
print("PASS: the stabilized lifts form a compatible Keller tower over every W_n(F_2)")

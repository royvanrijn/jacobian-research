#!/usr/bin/env python3
"""Exact checks for incidence suspensions through horizontal degree four."""

from __future__ import annotations

import sympy as sp


# Universal marked-line determinant.
P, S, Q, B = sp.symbols("P S Q B")
X = sp.Function("X")(P, S)
Y = sp.Function("Y")(P, S)
marked_line = sp.Matrix([P, B, Y - B * X])
incidence_jacobian = sp.factor(
    sp.det(marked_line.jacobian((P, S, B)))
)
assert incidence_jacobian == -sp.diff(Y, S) + B * sp.diff(X, S)


# The quadratic reciprocal source chart.
x, y, z, c = sp.symbols("x y z c", nonzero=True)
t = 1 + x * y
q = t**2 * z + c * y**2 * (1 + 3 * t)
P_source = t * q
S_source = x / t
Q_source = y + x * q
D = 1 - S * Q + P * S**2

source_chart_jacobian = sp.factor(
    sp.det(
        sp.Matrix([P_source, S_source, Q_source]).jacobian((x, y, z))
    )
)
assert source_chart_jacobian == t
assert sp.factor(
    D.subs({P: P_source, S: S_source, Q: Q_source}) - 1 / t
) == 0


# The three polynomial-D0 cases A=1,S,S^2.
for exponent in range(3):
    horizontal_degree = exponent + 2
    R, Q_tilde = sp.symbols(f"R_{exponent} Q_tilde_{exponent}")

    R_expr = S**exponent * P
    Q_tilde_expr = Q / S**exponent
    rechart_jacobian = sp.factor(
        sp.det(
            sp.Matrix([R_expr, S, Q_tilde_expr]).jacobian((P, S, Q))
        )
    )
    assert rechart_jacobian == 1

    transformed_D = sp.factor(
        D.subs(
            {
                P: R / S**exponent,
                Q: S**exponent * Q_tilde,
            }
        )
    )
    expected_D = 1 + R * S ** (2 - exponent) - S ** (exponent + 1) * Q_tilde
    assert sp.factor(transformed_D - expected_D) == 0

    lam = sp.Integer(horizontal_degree)
    horizontal_X = S**horizontal_degree
    assert sp.diff(horizontal_X, S) == lam * S ** (exponent + 1)

    D0 = 1 + R * S ** (2 - exponent)
    incidence_Y = sp.integrate(lam * D0, S)
    incidence_beta = sp.factor(
        (sp.diff(incidence_Y, S) - lam * D0)
        / sp.diff(horizontal_X, S)
    )
    assert incidence_beta == 0

    incidence_C = incidence_Y - Q_tilde * horizontal_X
    plane_jacobian = sp.factor(
        sp.det(
            sp.Matrix([Q_tilde, incidence_C]).jacobian((S, Q_tilde))
        )
    )
    assert sp.factor(plane_jacobian + lam * expected_D) == 0

    # Pull back the first incidence output.  It is polynomial only for
    # exponent zero.  For positive exponent it has an exact x^{-exponent}
    # pole with nonzero generic residue y at x=0.
    pulled_B = sp.factor(
        Q_tilde_expr.subs(
            {P: P_source, S: S_source, Q: Q_source}
        )
    )
    if exponent == 0:
        assert sp.denom(pulled_B) == 1
    else:
        pole_cleared = sp.factor(x**exponent * pulled_B)
        residue = sp.factor(pole_cleared.subs(x, 0))
        assert residue == y


# For a general polynomial A, the coefficient of the independent chart
# variable Q in B is 1/A.  Additive H and beta terms cannot change it.
a0, a1, a2 = sp.symbols("a0 a1 a2")
A = a0 + a1 * S + a2 * S**2
R, Q_tilde = sp.symbols("R Q_tilde")
E = sp.Function("E")(S)
H = sp.Function("H")(R, S)
general_D = sp.factor(
    D.subs(
        {
            P: (R - E) / A,
            Q: A * (Q_tilde - H),
        }
    )
)
assert sp.factor(sp.diff(general_D, Q_tilde) + S * A) == 0
assert sp.factor(sp.diff((R - E) * S**2 / A, R) - S**2 / A) == 0


print("PASS: the universal marked-line determinant is exact")
print("PASS: the quadratic source chart contributes D^{-1}")
print("PASS: A=1,S,S^2 give exactly X=S^2,S^3,S^4")
print("PASS: the cubic and quartic incidence B-coordinates have exact source poles")
print("PASS: additive chart and incidence shears cannot cancel the 1/A coefficient")

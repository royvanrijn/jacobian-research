"""Exact regressions for the one-additional-weight classification."""

from __future__ import annotations

import sympy as sp


x, y, z, Avar, W = sp.symbols("x y z A W")


def direct_jacobian(a: int, b: int, c: int, theta: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    """Differentiate a representative nonlinear seed and compare with (9)."""

    f = y**2 + y + 1
    g = y * Avar + Avar**2
    A = 1 + x * f
    B = A**b * z + g.subs(Avar, A)
    P = A**a * B
    Q = y + x * A**c * B
    h = c - a
    s = x * A**h
    source_w = sp.cancel(s * f)

    # R itself is unnecessary: replacing its third row by R_s ds at fixed
    # (P,Q) multiplies the exact (s,P,Q) Jacobian by Theta(W(s)).
    jac_spq = sp.det(sp.Matrix([s, P, Q]).jacobian([x, y, z]))
    lhs = sp.cancel(jac_spq * theta.subs(W, source_w))
    n = b + c - 1
    rhs = sp.cancel(-A**n * ((h + 1) * A - h) * theta.subs(W, (A - 1) * A**h))
    return sp.cancel(lhs), sp.cancel(rhs)


def bounded_weight_classification() -> int:
    """Search monic Theta of degree at most six in a bounded weight box."""

    solutions = 0
    for a in range(1, 6):
        for b in range(1, 6):
            for c in range(0, 6):
                h = c - a
                n = b + c - 1
                for degree in range(0, 7):
                    # Exhaust all monic degrees, rather than testing only the
                    # degree predicted by the theorem.
                    coefficients = sp.symbols(f"u0:{degree}")
                    theta = W**degree + sum(
                        coefficients[j] * W**j for j in range(degree)
                    )
                    expression = sp.together(
                        Avar**n
                        * ((h + 1) * Avar - h)
                        * theta.subs(W, (Avar - 1) * Avar**h)
                    )
                    numerator, denominator = sp.fraction(expression)
                    target = sp.symbols("target")
                    equation = sp.Poly(
                        sp.expand(numerator - target * denominator), Avar
                    )
                    solved = sp.solve(
                        equation.all_coeffs(), (*coefficients, target), dict=True
                    )

                    if h == -1 and degree == n:
                        expected_theta = sp.expand((W - 1) ** n)
                        expected_target = sp.Integer((-1) ** n)
                        assert solved == [
                            {
                                **{
                                    coefficients[j]: expected_theta.coeff(W, j)
                                    for j in range(degree)
                                },
                                target: expected_target,
                            }
                        ]
                        solutions += 1
                    else:
                        assert solved == []
    return solutions


for weights, theta in [
    ((1, 2, 0), 1 - 2 * W + W**2),
    ((2, 3, 1), (1 - W) ** 4),
    ((2, 1, 0), 1 + W + W**2),
    ((1, 1, 1), 2 - W),
]:
    direct, formula = direct_jacobian(*weights, theta)
    assert sp.cancel(direct - formula) == 0

count = bounded_weight_classification()
assert count > 0
print(f"PASS: three-weight Jacobian formula and {count} bounded rigid solutions")

#!/usr/bin/env python3
"""Exact regression for the finite differential system behind GMC(2).

For

    G(t,U) = Q(t,U)^(-1/2),   Q = (1-t*C)^2 - 4*t^2*D,
    M_k(t) = L(U^k G),        L(U^j) = j!,

integration by parts and 2*Q*G_U + Q_U*G = 0 give

    2 sum_j q_j M_(k+j)
      - sum_j (2k+j) q_j M_(k+j-1)
      = 2 delta_(k,0) Q(t,0)G(t,0),

where Q=sum_j q_j U^j.  Its leading term is 2*q_s*M_(k+s), so
all radial moments reduce to M_0,...,M_(s-1).  The identity

    2*Q*G_t + Q_t*G = 0

then gives an s-dimensional rational differential system.

This checker verifies the universal integrand identities, constructs the
four-dimensional system for an exact centered degree-(2,3) example, and
checks the system against its factorial moment series.
"""

from __future__ import annotations

import math

import sympy as sp


def factorial_functional(expression: sp.Expr, radial: sp.Symbol) -> sp.Expr:
    """Apply L(U^j)=j! to a polynomial in U."""

    return sp.expand(
        sum(
            coefficient * math.factorial(exponent[0])
            for exponent, coefficient in sp.Poly(
                sp.expand(expression), radial
            ).terms()
        )
    )


def radial_reductions(
    coefficients: list[sp.Expr],
    base_moments: tuple[sp.Symbol, ...],
    boundary: sp.Expr,
) -> list[sp.Expr]:
    """Reduce M_s,...,M_(2s-1) to the first s radial moments."""

    degree = len(coefficients) - 1
    assert len(base_moments) == degree
    assert coefficients[-1] != 0
    reductions: list[sp.Expr] = list(base_moments)
    for k in range(degree):
        right_side = 2 * boundary if k == 0 else sp.Integer(0)
        lower_shift = sum(
            (2 * k + j) * coefficients[j] * reductions[k + j - 1]
            for j in range(degree + 1)
            if k + j - 1 >= 0
        )
        same_shift = 2 * sum(
            coefficients[j] * reductions[k + j]
            for j in range(degree)
        )
        reductions.append(
            sp.cancel(
                (right_side + lower_shift - same_shift)
                / (2 * coefficients[-1])
            )
        )
    return reductions


def pfaffian_system(
    coefficients: list[sp.Expr],
    reductions: list[sp.Expr],
    base_moments: tuple[sp.Symbol, ...],
    parameter: sp.Symbol,
) -> tuple[sp.Matrix, sp.Matrix, tuple[sp.Symbol, ...]]:
    """Return A,b,x for the linear system A*x=b satisfied by M'."""

    degree = len(coefficients) - 1
    derivatives = sp.symbols(f"x0:{degree}")

    def total_derivative(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, parameter) + sum(
            sp.diff(expression, base_moments[index])
            * derivatives[index]
            for index in range(degree)
        )

    equations = []
    for k in range(degree):
        equations.append(
            sp.together(
                2
                * sum(
                    coefficients[j]
                    * total_derivative(reductions[k + j])
                    for j in range(degree + 1)
                )
                + sum(
                    sp.diff(coefficients[j], parameter)
                    * reductions[k + j]
                    for j in range(degree + 1)
                )
            )
        )
    matrix, right_side = sp.linear_eq_to_matrix(equations, derivatives)
    return matrix, right_side, derivatives


def main() -> None:
    radial, parameter = sp.symbols("U t")

    # Universal algebraic resolvent identities.
    q0, q1, q2, q3, q4 = sp.symbols("q0:5")
    generic_q = sum(
        coefficient * radial**power
        for power, coefficient in enumerate((q0, q1, q2, q3, q4))
    )
    generic_g = generic_q ** sp.Rational(-1, 2)
    assert sp.simplify(
        2 * generic_q * sp.diff(generic_g, radial)
        + sp.diff(generic_q, radial) * generic_g
    ) == 0
    time_coefficients = [
        sp.Function(f"q{power}")(parameter) for power in range(5)
    ]
    time_q = sum(
        coefficient * radial**power
        for power, coefficient in enumerate(time_coefficients)
    )
    time_g = time_q ** sp.Rational(-1, 2)
    assert sp.simplify(
        2 * time_q * sp.diff(time_g, parameter)
        + sp.diff(time_q, parameter) * time_g
    ) == 0

    # A centered degree-(2,3) pair.  Its resolvent has radial degree four.
    centered = radial**2 + radial - 3
    circuit = radial + 2 * radial**2 + radial**3
    q_expression = sp.expand(
        (1 - parameter * centered) ** 2
        - 4 * parameter**2 * circuit
    )
    q_polynomial = sp.Poly(q_expression, radial)
    degree = q_polynomial.degree()
    assert degree == 4
    coefficients = [
        q_polynomial.coeff_monomial(radial**power)
        for power in range(degree + 1)
    ]

    # The formal square-root branch has G(t,0)=1/(1+3t).
    boundary = 1 + 3 * parameter
    base_moments = sp.symbols(f"m0:{degree}")
    reductions = radial_reductions(
        coefficients, base_moments, boundary
    )
    assert len(reductions) == 2 * degree

    matrix, right_side, derivatives = pfaffian_system(
        coefficients, reductions, base_moments, parameter
    )
    assert matrix.shape == (degree, degree)
    assert sp.factor(matrix.det()) != 0
    solution = [
        sp.cancel(value)
        for value in matrix.inv() * right_side
    ]
    assert len(solution) == degree

    # Check both the radial recurrence and the differential system against
    # the exact factorial expansion.  Multiplying by each rational
    # denominator avoids making any analyticity assumption at t=0.
    cutoff = 10
    g_series = sp.series(
        q_expression ** sp.Rational(-1, 2),
        parameter,
        0,
        cutoff,
    ).removeO()
    moment_series = [
        factorial_functional(radial**k * g_series, radial)
        for k in range(2 * degree)
    ]

    for k in range(degree):
        left_side = (
            2
            * sum(
                coefficients[j] * moment_series[k + j]
                for j in range(degree + 1)
            )
            - sum(
                (2 * k + j)
                * coefficients[j]
                * moment_series[k + j - 1]
                for j in range(degree + 1)
                if k + j - 1 >= 0
            )
        )
        expected = 2 * boundary if k == 0 else 0
        assert (
            sp.series(
                left_side - expected, parameter, 0, cutoff
            ).removeO().expand()
            == 0
        )

    substitution = {
        base_moments[index]: moment_series[index]
        for index in range(degree)
    }
    for index, rational_derivative in enumerate(solution):
        numerator, denominator = sp.fraction(rational_derivative)
        residual = sp.expand(
            denominator
            * sp.diff(moment_series[index], parameter)
            - numerator.subs(substitution)
        )
        assert (
            sp.series(
                residual, parameter, 0, cutoff - 3
            ).removeO().expand()
            == 0
        )

    assert tuple(map(str, derivatives)) == ("x0", "x1", "x2", "x3")
    print("PASS resolvent: universal U- and t-differential identities")
    print("PASS resolvent: radial recurrence closes in dimension 4")
    print("PASS resolvent: exact Pfaffian system matches factorial series")


if __name__ == "__main__":
    main()

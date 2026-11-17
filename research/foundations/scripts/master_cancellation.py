#!/usr/bin/env python3
"""Exact algebra for the master cancellation construction.

The routines in this file deliberately work over ``QQ`` and finite etale
``QQ``-algebras.  They are shared by the universal verifier, the displayed
instance verifier, and the command-line regression generator.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import sympy as sp


def integrate_unit_interval(polynomial: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Integrate a polynomial from zero to one without analytic machinery."""
    poly = sp.Poly(sp.expand(polynomial), variable)
    return sp.factor(
        sum(coefficient / sp.Rational(exponent + 1)
            for (exponent,), coefficient in poly.terms())
    )


def raw_parameter_polynomial(m: int, r: int, q: sp.Symbol) -> sp.Expr:
    """Return Phi(0,q), with the normalization inherited from the integral."""
    n = m * r
    return sp.factor(sum(
        (-1) ** k
        * sp.binomial(n, k)
        * sp.factorial(r)
        * sp.factorial(k)
        / sp.factorial(r + k + 1)
        * q**k
        for k in range(n + 1)
    ))


def parameter_polynomial(m: int, r: int, q: sp.Symbol) -> sp.Expr:
    """Return the monic truncated-binomial parameter polynomial M_{m,r}."""
    n = m * r
    total = n + r + 1
    return sp.expand(sum(
        (-1) ** j * sp.binomial(total, j) * q ** (n - j)
        for j in range(n + 1)
    ))


def parameter_discriminant(m: int, r: int) -> sp.Integer:
    """Return the closed-form discriminant of M_{m,r}."""
    n = m * r
    total = n + r + 1
    leading_reciprocal = sp.binomial(n + r, r)
    return sp.Integer(
        (-1) ** (n * (n - 1) // 2)
        * (r + 1)
        * total ** (n - 1)
        * leading_reciprocal ** (n - 2)
    )


def parameter_discriminant_is_square(m: int, r: int) -> bool:
    """Decide whether the parameter discriminant is a square in QQ."""
    n = m * r
    if n % 4 not in (0, 1):
        return False
    if n % 2 == 0:
        square_part = (r + 1) * (n + r + 1)
    else:
        square_part = (r + 1) * sp.binomial(n + r, r)
    return bool(sp.integer_nthroot(int(square_part), 2)[1])


def even_square_discriminant_family(r: int, k: int) -> int:
    """Return an m in an infinite square-discriminant family for fixed r."""
    if r < 1 or k < 1:
        raise ValueError("r and k must be positive")
    factors = sp.factorint(r + 1)
    squarefree = math.prod(
        int(prime) for prime, exponent in factors.items() if exponent % 2
    )
    square_root = sp.integer_nthroot((r + 1) // squarefree, 2)
    assert square_root[1]
    a = int(square_root[0])
    return 4 * squarefree * k * (a + r * k)


def phi(m: int, r: int, A: sp.Symbol, H: sp.Expr) -> sp.Expr:
    """The polynomial numerator whose A-adic vanishing is cancellation."""
    u = sp.Dummy("u")
    bracket = (
        A
        + (1 - A) * u
        * (1 - (1 - A) * H * (1 - u)) ** m
    )
    return sp.expand(integrate_unit_interval(bracket**r, u))


def cancellation_operator(m: int, r: int, A: sp.Symbol, H: sp.Expr) -> sp.Expr:
    """Return L_{m,r}(H): the remainder of Phi(A,H) modulo A^(r+1)."""
    expression = phi(m, r, A, H)
    return sp.expand(sum(
        expression.coeff(A, degree) * A**degree for degree in range(r + 1)
    ))


def reduce_q(expression: sp.Expr, q: sp.Symbol, modulus: sp.Expr) -> sp.Expr:
    """Reduce a rational-coefficient expression in q modulo a monic modulus."""
    return sp.rem(sp.Poly(sp.cancel(expression), q), sp.Poly(modulus, q)).as_expr()


def hensel_jet(m: int, r: int, A: sp.Symbol, q: sp.Symbol) -> sp.Expr:
    """Construct the unique degree-r cancellation jet with constant term q.

    The answer lives in QQ[q]/(M_{m,r}).  Squarefreeness of M makes the
    derivative of Phi(0,q) invertible in this algebra.
    """
    modulus = parameter_polynomial(m, r, q)
    raw = raw_parameter_polynomial(m, r, q)
    inverse_derivative = sp.invert(
        sp.Poly(sp.diff(raw, q), q), sp.Poly(modulus, q)
    ).as_expr()
    H = q
    for degree in range(1, r + 1):
        residual = sp.expand(phi(m, r, A, H)).coeff(A, degree)
        coefficient = reduce_q(-residual * inverse_derivative, q, modulus)
        H += coefficient * A**degree
    return sp.collect(sp.expand(H), A)


def fiber_antiderivative(
    m: int,
    r: int,
    T: sp.Symbol,
    P: sp.Expr,
    Q: sp.Expr,
) -> sp.Expr:
    """Return integral_0^T (1-t(Q-Pt)^m)^r dt exactly."""
    t = sp.Dummy("t")
    integrand = sp.expand((1 - t * (Q - P * t) ** m) ** r)
    poly = sp.Poly(integrand, t)
    return sp.expand(sum(
        coefficient * T ** (exponent + 1) / sp.Rational(exponent + 1)
        for (exponent,), coefficient in poly.terms()
    ))


@dataclass(frozen=True)
class DisplayedInstance:
    label: str
    m: int
    r: int


DISPLAYED_INSTANCES = (
    DisplayedInstance("N_1", 1, 1),
    DisplayedInstance("N_2", 2, 1),
    DisplayedInstance("H_2", 1, 2),
    DisplayedInstance("M_2_2", 2, 2),
)

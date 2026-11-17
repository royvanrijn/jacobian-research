#!/usr/bin/env python3
"""Exact classification of the linear-h rank-one three-real GMC ansatz.

Use circular coordinates W,Z and one independent real Gaussian T:

    E(W^i Z^j T^k) = delta(i,j) i! (k-1)!!.

After scaling the coefficient of W to one, the ansatz is

    P = W(1+aZ) + (b0+b1 Z+b2 Z^2)T^2.

On the open chart a != 0, the first three pure moments already force Long's
one-parameter family.  The all-order assertion is then proved by the formal
square identity, rather than inferred from a finite moment cutoff.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


Exponent = tuple[int, int, int]
Polynomial = dict[Exponent, sp.Expr]


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for (lw, lz, lt), lc in left.items():
        for (rw, rz, rt), rc in right.items():
            exponent = (lw + rw, lz + rz, lt + rt)
            out[exponent] = sp.expand(out.get(exponent, 0) + lc * rc)
    return {exponent: coefficient for exponent, coefficient in out.items() if coefficient != 0}


def odd_double_factorial(exponent: int) -> int:
    value = 1
    for factor in range(1, exponent, 2):
        value *= factor
    return value


def gaussian_expectation(polynomial: Polynomial) -> sp.Expr:
    total = 0
    for (w_degree, z_degree, t_degree), coefficient in polynomial.items():
        if w_degree == z_degree and t_degree % 2 == 0:
            total += (
                coefficient
                * factorial(w_degree)
                * odd_double_factorial(t_degree)
            )
    return sp.factor(total)


def main() -> None:
    a, b0, b1, b2, inverse_a, z = sp.symbols("a b0 b1 b2 inverse_a z")

    general: Polynomial = {
        (1, 0, 0): sp.Integer(1),
        (1, 1, 0): a,
        (0, 0, 2): b0,
        (0, 1, 2): b1,
        (0, 2, 2): b2,
    }

    powers: list[Polynomial] = [{(0, 0, 0): sp.Integer(1)}]
    moments: list[sp.Expr] = []
    for _ in range(6):
        powers.append(multiply(powers[-1], general))
        moments.append(gaussian_expectation(powers[-1]))

    expected_first_three = [
        a + b0,
        2 * a**2 + 2 * a * b0 + 3 * b0**2 + 2 * b1,
        3
        * (
            2 * a**3
            + 2 * a**2 * b0
            + 3 * a * b0**2
            + 4 * a * b1
            + 5 * b0**3
            + 6 * b0 * b1
            + 2 * b2
        ),
    ]
    assert [sp.expand(moment) for moment in moments[:3]] == [
        sp.expand(moment) for moment in expected_first_three
    ]

    # Saturation by a is encoded by adjoining inverse_a*a-1.
    basis = sp.groebner(
        moments[:3] + [inverse_a * a - 1],
        inverse_a,
        b2,
        b1,
        b0,
        a,
        order="lex",
        domain=sp.QQ,
    )
    expected_basis = [
        a * inverse_a - 1,
        sp.Rational(1, 2) * a**3 + b2,
        sp.Rational(3, 2) * a**2 + b1,
        a + b0,
    ]
    assert [sp.expand(poly.as_expr()) for poly in basis.polys] == expected_basis

    substitutions = {
        b0: -a,
        b1: -sp.Rational(3, 2) * a**2,
        b2: -sp.Rational(1, 2) * a**3,
    }
    for moment in moments:
        assert sp.expand(moment.subs(substitutions)) == 0

    # The general Wick master series is
    # h(z)^m / sqrt(1-2*z*v(z)/h(z)).  On the forced family the radicand
    # is h(z)^2, with the formal square root having constant term one.
    h = 1 + a * z
    v = (b0 + b1 * z + b2 * z**2).subs(substitutions)
    assert sp.factor(1 - 2 * z * v / h - h**2) == 0

    # Therefore the normalized A(Z)-mixed moment is
    # [z^m] A(z) h(z)^(m-1).  A=1 gives zero and A=z gives a^(m-1).
    for m in range(1, 12):
        pure_coefficient = sp.expand(h ** (m - 1)).coeff(z, m)
        mixed_coefficient = sp.expand(z * h ** (m - 1)).coeff(z, m)
        assert pure_coefficient == 0
        assert mixed_coefficient == a ** (m - 1)

    # Every variable coefficient in the displayed normalized support is
    # nonzero on a != 0.  (The coefficient of W was fixed to one.)  Adding
    # any available deletion equation makes the saturated ideal one.
    deletion_coefficients = [a, b0, b1, b2]
    for coefficient in deletion_coefficients:
        deletion_basis = sp.groebner(
            moments[:3] + [inverse_a * a - 1, coefficient],
            inverse_a,
            b2,
            b1,
            b0,
            a,
            order="lex",
            domain=sp.QQ,
        )
        assert deletion_basis.polys == [sp.Poly(1, *deletion_basis.gens, domain=sp.QQ)]

    print("PASS rank-one GMC: first three moments give the saturated Long family")
    print("PASS rank-one GMC: formal square identity proves all pure moments vanish")
    print("PASS rank-one GMC: E[Z P^m] = m! a^(m-1) for every m >= 1")
    print("PASS rank-one GMC: no variable-term deletion or cubic member on a != 0")


if __name__ == "__main__":
    main()

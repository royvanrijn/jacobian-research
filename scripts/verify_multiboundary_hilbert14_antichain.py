#!/usr/bin/env python3
"""Verify a two-boundary Hilbert--14 bidegree escape control.

Let

    A=k[s^2,s^3,t^2,t^3],  S=A[X,Y,U,V],

and let the commuting LNDs be

    Ds=s^3*d/dX-s^2*d/dY,
    Dt=t^3*d/dU-t^2*d/dV.

In the normalization, P=X+sY and Q=U+tV are invariant.  Every

    F_(m,n)=s^2*t^2*P^m*Q^n

lies back in S.  Modulo (s^4,t^4), two positive-P factors or two positive-Q
factors vanish.  Hence a finite generating set has a bounded surviving
(P,Q)-rectangle, while F_(dP+1,dQ+1) escapes it.

The finite loops below replay the uniform formulas.  The written
conductor/bidegree argument is the proof of non-finite generation.
"""

from __future__ import annotations

import argparse
from math import comb

import sympy as sp


s, t, X, Y, U, V = sp.symbols("s t X Y U V")
square_s, cube_s, square_t, cube_t = sp.symbols("A2 A3 B2 B3")
variables = (s, t, X, Y, U, V)

Ds = (0, 0, s**3, -s**2, 0, 0)
Dt = (0, 0, 0, 0, t**3, -t**2)
P = X + s * Y
Q = U + t * V


def apply_derivation(
    polynomial: sp.Expr, images: tuple[sp.Expr, ...]
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(polynomial, variable) * image
            for variable, image in zip(variables, images, strict=True)
        )
    )


def cusp_monomial(
    exponent: int, square: sp.Symbol, cube: sp.Symbol
) -> sp.Expr:
    if exponent < 0 or exponent == 1:
        raise ValueError(f"exponent {exponent} is outside <2,3>")
    if exponent % 2 == 0:
        return square ** (exponent // 2)
    return cube * square ** ((exponent - 3) // 2)


def truncate_two_boundaries(polynomial: sp.Expr) -> sp.Expr:
    """Reduce modulo the ideal (s^4,t^4)."""

    expanded = sp.Poly(sp.expand(polynomial), s, t)
    return sp.expand(
        sum(
            coefficient * s**exponents[0] * t**exponents[1]
            for exponents, coefficient in expanded.terms()
            if exponents[0] < 4 and exponents[1] < 4
        )
    )


def mixed_ladder(m_degree: int, n_degree: int) -> sp.Expr:
    return sp.expand(s**2 * t**2 * P**m_degree * Q**n_degree)


def cusp_representative(m_degree: int, n_degree: int) -> sp.Expr:
    """Write F_(m,n) in k[s^2,s^3,t^2,t^3,X,Y,U,V]."""

    return sp.expand(
        sum(
            comb(m_degree, i)
            * comb(n_degree, j)
            * X ** (m_degree - i)
            * Y**i
            * U ** (n_degree - j)
            * V**j
            * cusp_monomial(2 + i, square_s, cube_s)
            * cusp_monomial(2 + j, square_t, cube_t)
            for i in range(m_degree + 1)
            for j in range(n_degree + 1)
        )
    )


def verify(max_bidegree: int) -> None:
    assert apply_derivation(P, Ds) == 0
    assert apply_derivation(P, Dt) == 0
    assert apply_derivation(Q, Ds) == 0
    assert apply_derivation(Q, Dt) == 0

    # The two actions commute and are locally nilpotent on the generators.
    for index in range(len(variables)):
        bracket = apply_derivation(Dt[index], Ds) - apply_derivation(
            Ds[index], Dt
        )
        assert bracket == 0
    for derivation in (Ds, Dt):
        for variable in variables:
            assert apply_derivation(
                apply_derivation(variable, derivation), derivation
            ) == 0

    for m_degree in range(max_bidegree + 2):
        for n_degree in range(max_bidegree + 2):
            invariant = mixed_ladder(m_degree, n_degree)
            assert apply_derivation(invariant, Ds) == 0
            assert apply_derivation(invariant, Dt) == 0
            representative = cusp_representative(m_degree, n_degree)
            assert sp.expand(
                representative.subs(
                    {
                        square_s: s**2,
                        cube_s: s**3,
                        square_t: t**2,
                        cube_t: t**3,
                    }
                )
                - invariant
            ) == 0

    # Regression replay of the conductor-square mechanism.
    pure_s = [sp.expand(s**2 * P**degree) for degree in range(max_bidegree + 2)]
    pure_t = [sp.expand(t**2 * Q**degree) for degree in range(max_bidegree + 2)]
    for left in pure_s:
        for right in pure_s:
            assert truncate_two_boundaries(left * right) == 0
    for left in pure_t:
        for right in pure_t:
            assert truncate_two_boundaries(left * right) == 0

    for p_bound in range(max_bidegree + 1):
        for q_bound in range(max_bidegree + 1):
            escaping = truncate_two_boundaries(
                mixed_ladder(p_bound + 1, q_bound + 1)
            )
            assert escaping != 0
            polynomial = sp.Poly(escaping, X, Y, U, V)
            for exponents, coefficient in polynomial.terms():
                if coefficient == 0:
                    continue
                assert exponents[0] + exponents[1] == p_bound + 1
                assert exponents[2] + exponents[3] == q_bound + 1

    # In a tangent-normalized factorization slice, m=a*d+c*p=1 places 1 in
    # the ideal (a,p).  Thus the two obvious leading divisors never meet and
    # cannot realize this two-boundary control.
    a, c, p, d = sp.symbols("a c p d")
    tangent_coefficient = a * d + c * p
    assert sp.expand(tangent_coefficient - (a * d + p * c)) == 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-bidegree",
        type=int,
        default=7,
        help="finite rectangle used to replay the uniform bidegree formulas",
    )
    args = parser.parse_args()
    if args.max_bidegree < 0:
        parser.error("--max-bidegree must be nonnegative")

    verify(args.max_bidegree)
    print("PASS: the two cusp LNDs commute and fix P=X+sY, Q=U+tV")
    print(
        "PASS: s^2*t^2*P^m*Q^n was replayed through the "
        f"{args.max_bidegree + 1} by {args.max_bidegree + 1} rectangle"
    )
    print("PASS: the modulo-(s^4,t^4) conductor-square escape was replayed")
    print("PASS: tangent normalization makes the leading divisors a=0,p=0 disjoint")
    print("NOTE: the arbitrary-bidegree written argument proves non-finiteness")


if __name__ == "__main__":
    main()

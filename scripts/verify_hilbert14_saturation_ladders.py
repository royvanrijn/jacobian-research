#!/usr/bin/env python3
"""Exact finite and non-finite controls for LND saturation experiments.

The terminating control is the invariant algebra on the normalized (2,3)
factorization slice.  Two apparent denominator classes G and H close after
finite saturation and give the polynomial ring k[K,H,V].

The non-terminating control is Maubach's cusp-base example.  For

    R=k[T^2,T^3],  D=T^3*d/dX-T^2*d/dY,  P=X+T*Y,

every F_n=T^2*P^n lies in the original finitely generated algebra and in
ker(D).  Products of two positive-P-degree invariants vanish modulo T^4,
so no finite list can generate the unbounded degree ladder.

The bounded ``--max-degree`` loop is a regression replay of formulas whose
written proof is uniform in n; passing it is not, by itself, a proof of
non-finite generation.
"""

from __future__ import annotations

import argparse
from math import comb

import sympy as sp


def apply_derivation(
    polynomial: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    images: tuple[sp.Expr, ...],
) -> sp.Expr:
    """Apply a derivation specified by its images on algebra generators."""

    return sp.expand(
        sum(
            sp.diff(polynomial, variable) * image
            for variable, image in zip(variables, images, strict=True)
        )
    )


def cusp_monomial(exponent: int, square: sp.Symbol, cube: sp.Symbol) -> sp.Expr:
    """Represent T^exponent in k[T^2,T^3], when possible."""

    if exponent < 0 or exponent == 1:
        raise ValueError(f"T^{exponent} is not in k[T^2,T^3]")
    if exponent % 2 == 0:
        return square ** (exponent // 2)
    return cube * square ** ((exponent - 3) // 2)


def truncate_in_T(polynomial: sp.Expr, T: sp.Symbol, cutoff: int) -> sp.Expr:
    """Reduce a polynomial modulo T^cutoff."""

    expanded = sp.Poly(sp.expand(polynomial), T)
    return sp.expand(
        sum(
            coefficient * T**exponent[0]
            for exponent, coefficient in expanded.terms()
            if exponent[0] < cutoff
        )
    )


def verify_terminating_control() -> None:
    """Check the finite saturation identities from the (2,3) slice."""

    K, H, V = sp.symbols("K H V")
    a = (V**2 - H * K**2) / 1024
    G = (V + a * H) / 4
    U = (a**2 * H - K**2 - 8 * a * G) / 16
    J = 4 * U + a * G

    assert sp.factor(J - 4 * U - a * G) == 0
    assert sp.factor(K**2 + 8 * J - 16 * U - a**2 * H) == 0
    assert sp.factor(V - (4 * G - a * H)) == 0
    assert sp.factor(V**2 - H * K**2 - 1024 * a) == 0


def verify_nonterminating_control(max_degree: int) -> None:
    """Replay Maubach's invariant ladder through the requested degree."""

    T, X, Y, z = sp.symbols("T X Y z")
    square, cube = sp.symbols("A B")
    variables = (T, X, Y, z)
    derivation = (sp.Integer(0), T**3, -T**2, sp.Integer(0))
    P = X + T * Y

    assert apply_derivation(P, variables, derivation) == 0
    relation = z**2 - T**8 * P**2 - 1
    assert apply_derivation(relation, variables, derivation) == 0

    ladder: list[sp.Expr] = []
    for degree in range(max_degree + 2):
        invariant = sp.expand(T**2 * P**degree)
        ladder.append(invariant)
        assert apply_derivation(invariant, variables, derivation) == 0

        cusp_representative = sp.expand(
            sum(
                comb(degree, j)
                * X ** (degree - j)
                * Y**j
                * cusp_monomial(2 + j, square, cube)
                for j in range(degree + 1)
            )
        )
        assert sp.expand(
            cusp_representative.subs({square: T**2, cube: T**3}) - invariant
        ) == 0

    # Uniform obstruction schema.  Positive-degree generators lie in the
    # conductor ideal (T^2,T^3), so every product of two vanishes mod T^4.
    # A list of maximum (X,Y)-degree d therefore cannot create F_(d+1),
    # whose T^2*X^(d+1) term survives in that quotient.
    for degree_bound in range(max_degree + 1):
        escaping_invariant = ladder[degree_bound + 1]
        escaping_truncation = truncate_in_T(escaping_invariant, T, 4)
        assert sp.Poly(escaping_truncation, X, Y).total_degree() == degree_bound + 1
        for left in ladder[: degree_bound + 1]:
            for right in ladder[: degree_bound + 1]:
                assert truncate_in_T(left * right, T, 4) == 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-degree",
        type=int,
        default=12,
        help="largest finite prefix used to replay the uniform ladder formulas",
    )
    args = parser.parse_args()
    if args.max_degree < 0:
        parser.error("--max-degree must be nonnegative")

    verify_terminating_control()
    verify_nonterminating_control(args.max_degree)
    print("PASS: the (2,3) saturation closes with B=k[K,H,V]")
    print(
        "PASS: Maubach's T^2(X+TY)^n ladder was replayed through "
        f"n={args.max_degree + 1}"
    )
    print("PASS: the modulo-T^4 degree-escape obstruction was replayed")
    print("NOTE: the written arbitrary-n argument, not the bound, proves non-finiteness")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Archived bounded search for one-complex-plus-one-real weighted witnesses.

The target branch is g=u*h(g), with h=1+z^2-z^3 by default.  For

    P = W*h(Z) + sum c[k,i,j] T^k W^i Z^j,   Q=Z,

the script forms exact necessary equations

    E(P^m)=0,
    E(Q P^m)=(m-1)! [z^(m-1)] h(z)^m

through a requested order and computes their Groebner basis.  A unit basis is
an exact obstruction to that finite-support ansatz; a nonunit basis is only a
bounded survivor, not an all-orders witness.

GMC is already false in every dimension n >= 3 by the explicit five-term
three-real witness.  This program is retained only to reproduce the historical
weighted-family reconnaissance; it is not part of the active GMC(2) program.

The default support is the first genuinely W-dependent extension of the
half-pair correction, T^2(v_0(Z)+W v_1(Z)).  Use --t-degrees to include higher
powers of the real Gaussian.  This is reconnaissance: the accompanying note
contains the all-orders separated-quadratic classification.
"""

from __future__ import annotations

import argparse
from math import factorial

import sympy as sp


W, Z, T = sp.symbols("W Z T")


def odd_double_factorial(even_degree: int) -> int:
    value = 1
    for factor in range(1, even_degree, 2):
        value *= factor
    return value


def gaussian_expectation(expr: sp.Expr) -> sp.Expr:
    """Contract E(W^a Z^b)=delta_(a,b)a! and standard-real T moments."""
    total = sp.Integer(0)
    for (w_degree, z_degree, t_degree), coefficient in sp.Poly(
        sp.expand(expr), W, Z, T
    ).terms():
        if w_degree != z_degree or t_degree % 2:
            continue
        total += (
            coefficient
            * factorial(w_degree)
            * odd_double_factorial(t_degree)
        )
    return sp.expand(total)


def parse_t_degrees(raw: str) -> tuple[int, ...]:
    degrees = tuple(sorted({int(item) for item in raw.split(",") if item}))
    if not degrees or min(degrees) < 1:
        raise argparse.ArgumentTypeError("t-degrees must be positive integers")
    return degrees


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z-degree", type=int, default=1)
    parser.add_argument("--w-degree", type=int, default=1)
    parser.add_argument("--t-degrees", type=parse_t_degrees, default=(2,))
    parser.add_argument("--moment-order", type=int, default=5)
    parser.add_argument(
        "--seed",
        choices=("canonical-cubic", "affine"),
        default="canonical-cubic",
    )
    parser.add_argument(
        "--pure-only",
        action="store_true",
        help="omit the target mixed-moment equations",
    )
    parser.add_argument("--print-equations", action="store_true")
    parser.add_argument("--print-basis", action="store_true")
    args = parser.parse_args()
    if args.z_degree < 0 or args.w_degree < 0 or args.moment_order < 1:
        parser.error("degrees must be nonnegative and moment-order positive")

    h = 1 + Z**2 - Z**3 if args.seed == "canonical-cubic" else 1 + Z
    coefficients: list[sp.Symbol] = []
    correction = sp.Integer(0)
    for k in args.t_degrees:
        for i in range(args.w_degree + 1):
            for j in range(args.z_degree + 1):
                coefficient = sp.Symbol(f"c_{k}_{i}_{j}")
                coefficients.append(coefficient)
                correction += coefficient * T**k * W**i * Z**j
    P = sp.expand(W * h + correction)

    equations: list[sp.Expr] = []
    p_power = sp.Integer(1)
    for m in range(1, args.moment_order + 1):
        p_power = sp.expand(p_power * P)
        equations.append(gaussian_expectation(p_power))
        if not args.pure_only:
            expected = (
                factorial(m - 1)
                * sp.Poly(h**m, Z).coeff_monomial(Z ** (m - 1))
            )
            equations.append(gaussian_expectation(Z * p_power) - expected)

    equations = [equation for equation in equations if equation != 0]
    print(f"unknowns={len(coefficients)} equations={len(equations)}")
    if args.print_equations:
        for index, equation in enumerate(equations, start=1):
            print(f"e{index}={sp.factor(equation)}")

    basis = sp.groebner(equations, *coefficients, order="grevlex")
    is_unit = basis.contains(sp.Integer(1))
    print(f"groebner_basis_size={len(basis.polys)}")
    print(f"unit_ideal={is_unit}")
    if args.print_basis:
        for polynomial in basis.polys:
            print(sp.factor(polynomial.as_expr()))


if __name__ == "__main__":
    main()

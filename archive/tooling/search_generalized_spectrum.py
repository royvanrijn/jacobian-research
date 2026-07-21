#!/usr/bin/env python3
"""Exact stress test for generalized-cancellation spectral intersections.

The all-degree coprimality theorem is proved by disjoint root annuli; this
bounded search remains an independent regression of its conclusion.
"""
from __future__ import annotations

import argparse

import sympy as sp


def spectral_polynomial(n: int, exponent: int, q: sp.Symbol) -> sp.Poly:
    expression = sum(
        (-q) ** j
        * sp.binomial(n, j)
        / ((exponent + j + 1) * sp.binomial(exponent + j, j))
        for j in range(n + 1)
    )
    return sp.Poly(expression, q, domain=sp.QQ)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-e", type=int, default=12)
    parser.add_argument("--max-n", type=int, default=80)
    args = parser.parse_args()
    if args.max_e < 1 or args.max_n < 1:
        parser.error("bounds must be positive")

    q = sp.symbols("q")
    intersections: list[tuple[int, int, int, sp.Expr]] = []
    for exponent in range(1, args.max_e + 1):
        polynomials: list[sp.Poly] = []
        for n in range(1, args.max_n + 1):
            current = spectral_polynomial(n, exponent, q)
            for previous_n, previous in enumerate(polynomials, start=1):
                gcd = sp.gcd(current, previous)
                if gcd.degree() > 0:
                    intersections.append(
                        (exponent, previous_n, n, sp.factor(gcd.as_expr()))
                    )
            polynomials.append(current)
        print(f"e={exponent}: checked N<= {args.max_n}")

    if intersections:
        for exponent, left, right, gcd in intersections:
            print(f"COMMON e={exponent}, N={left},{right}: {gcd}")
        raise SystemExit(1)
    print(
        "PASS: no pairwise common factor for "
        f"e<= {args.max_e} and N<= {args.max_n}"
    )


if __name__ == "__main__":
    main()

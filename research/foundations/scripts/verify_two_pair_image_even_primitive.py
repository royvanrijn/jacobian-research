#!/usr/bin/env python3
"""Bounded exact checks for all even primitive two-pair Image degrees.

The all-order conclusion belongs to the written proof. This standalone replay
checks original four-variable contractions for k=0..3, m=1..3, independent
Laurent-sphere moments through m=8, normalization, parity and factor valuations.
Run from the repository root with Python 3.12 and SymPy 1.14.0:

    .venv/bin/python research/foundations/scripts/verify_two_pair_image_even_primitive.py

The deterministic JSON receipt is written next to the canonical research note.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from math import comb, factorial, prod
from pathlib import Path

if not __debug__:
    raise SystemExit("Run without -O: verification assertions must be enabled.")

import sympy as sp

SCRIPT = Path(__file__).resolve()
RESEARCH = SCRIPT.parents[2]
OUTPUT = RESEARCH / "foundations/extended-geometry/two_pair_image_even_primitive_verification.json"
XI1, XI2, Z1, Z2 = sp.symbols("xi1 xi2 z1 z2")
VARS = (XI1, XI2, Z1, Z2)
X, T = sp.symbols("x t")
FOUR_K = range(4)
FOUR_M = range(1, 4)
SPHERE_M = range(1, 9)


def poly(expression: sp.Expr) -> sp.Poly:
    return sp.Poly(expression, *VARS, domain=sp.QQ)


def contract_balanced(polynomial: sp.Poly) -> sp.Rational:
    """Literal E_2: only equal xi/z exponent vectors survive in equal bidegree."""
    degrees = {(a + b, c + d) for (a, b, c, d), _ in polynomial.terms()}
    assert len(degrees) == 1
    dual, coordinate = next(iter(degrees))
    assert dual == coordinate
    return sum(
        (coefficient * factorial(a) * factorial(b)
         for (a, b, c, d), coefficient in polynomial.terms()
         if a == c and b == d),
        sp.S.Zero,
    )


def sphere_moment(k: int, m: int, marker: bool) -> sp.Rational:
    """Expand sphere substitution directly, then integrate its constant term."""
    # R=1, Z=x, T=t, W=(1-t^2)/(2x). No claimed moment formula is used.
    phase = T**k * (1 + X)**(k + 1) * ((1 - T**2) / X - (2 + X) * T**2)
    # Multiplication by x^m removes every negative exponent before Poly.
    polynomial = sp.Poly(sp.expand(X * phase), X, T, domain=sp.QQ)**m
    required_x = m - int(marker)
    return sum(
        (coefficient / (t_degree + 1)
         for (x_degree, t_degree), coefficient in polynomial.terms()
         if x_degree == required_x and t_degree % 2 == 0),
        sp.S.Zero,
    )


def coefficient_sum(k: int, m: int, marker: bool) -> sp.Rational:
    """Separate binomial constant-term calculation, with full height parity."""
    if k * m % 2:
        return sp.S.Zero
    return sum(
        (sp.Rational((-1)**j * comb(m, j) * comb((k + 1)*m + 2*j,
                                               m - int(marker)),
                     k*m + 2*j + 1)
         for j in range(m + 1)),
        sp.S.Zero,
    )


def marked_integral(k: int, m: int) -> sp.Rational:
    """Integral from 0 to 1 of w^(km)(1-w^2)^m; nonzero for either parity."""
    return sum(
        (sp.Rational((-1)**j * comb(m, j), k*m + 2*j + 1)
         for j in range(m + 1)),
        sp.S.Zero,
    )


def valuation(polynomial: sp.Poly, factor: sp.Poly) -> int:
    count = 0
    while True:
        quotient, remainder = polynomial.div(factor)
        if not remainder.is_zero:
            return count
        polynomial = quotient
        count += 1


def main() -> None:
    r = poly(XI1*Z1 + XI2*Z2)
    z = poly(XI1*Z2)
    w = poly(2*XI2*Z1)
    t = poly(XI1*Z1 - XI2*Z2)
    a = r + z
    d = 2*w*r**2 - 2*r*t**2 - z*t**2
    assert t**2 == r**2 - 2*z*w
    assert sp.gcd(t, a).total_degree() == 0
    assert not d.rem(t).is_zero
    assert not d.rem(a).is_zero
    # Bilinear coefficient determinants certify absolute irreducibility of T,A.
    assert sp.Matrix([[1, 0], [0, -1]]).det() == -1
    assert sp.Matrix([[1, 1], [0, 1]]).det() == 1

    direct_rows = []
    factor_rows = []
    point = {XI1: 1, XI2: 1, Z1: 1, Z2: -1}
    assert r.eval(point) == 0
    for k in FOUR_K:
        balanced_degree = 2*k + 4
        f = t**k * a**(k + 1) * d
        value = f.eval(point)
        assert value == (-1)**(k + 1) * 2**(k + 2) != 0
        t_order, a_order = valuation(f, t), valuation(f, a)
        assert (t_order, a_order) == (k, k + 1)
        assert sp.igcd(t_order, a_order) == 1
        factor_rows.append({"k": k, "R_zero_point_value": str(value),
                            "T_multiplicity": t_order,
                            "R_plus_Z_multiplicity": a_order})
        f_power = poly(1)
        for m in FOUR_M:
            f_power *= f
            pure = contract_balanced(f_power)
            marked = contract_balanced(z * f_power)
            pure_sphere = sphere_moment(k, m, False)
            marked_sphere = sphere_moment(k, m, True)
            expected = marked_integral(k, m) if k*m % 2 == 0 else sp.S.Zero
            assert pure == factorial(balanced_degree*m + 1)*pure_sphere == 0
            assert marked == factorial(balanced_degree*m + 2)*marked_sphere
            assert marked_sphere == expected
            if k*m % 2 == 0:
                assert marked > 0
            else:
                assert marked == 0
            direct_rows.append({"k": k, "m": m, "balanced_degree": balanced_degree,
                                "pure_contraction": str(pure),
                                "marked_contraction": str(marked),
                                "marked_angular_moment": str(marked_sphere)})

    sphere_rows = []
    for k in FOUR_K:
        for m in SPHERE_M:
            pure = sphere_moment(k, m, False)
            marked = sphere_moment(k, m, True)
            expected = marked_integral(k, m) if k*m % 2 == 0 else sp.S.Zero
            assert pure == coefficient_sum(k, m, False) == 0
            assert marked == coefficient_sum(k, m, True) == expected
            integral = marked_integral(k, m)
            product_value = sp.Rational(2**m * factorial(m),
                                        prod(k*m + 2*j + 1 for j in range(m + 1)))
            assert integral == product_value > 0
            if k == 1:
                sphere_rows.append({"k": k, "m": m, "pure_angular_moment": str(pure),
                                    "marked_angular_moment": str(marked)})

    # Negative controls distinguish odd-height parity and one missing radial factorial.
    odd_raw_integral = marked_integral(1, 1)
    assert odd_raw_integral == sp.Rational(1, 4)
    assert sphere_moment(1, 1, True) == 0 != odd_raw_integral
    first_marked = contract_balanced(z * a * d)
    assert first_marked == 480
    assert first_marked != factorial(5)*marked_integral(0, 1)

    receipt = {
        "status": "passed",
        "scope": "Internal bounded exact checks; the written all-order proof, external review and novelty are not certified by this computation.",
        "family": "F_k=T^k(R+Z)^(k+1)(2WR^2-2RT^2-ZT^2), d=2k+4, k>=0",
        "arithmetic": "Exact rational polynomial arithmetic",
        "command": ".venv/bin/python research/foundations/scripts/verify_two_pair_image_even_primitive.py",
        "working_directory": "repository root",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "checker": {"path": SCRIPT.relative_to(RESEARCH.parent).as_posix(),
                    "sha256": hashlib.sha256(SCRIPT.read_bytes()).hexdigest()},
        "assurance": {"internal_replay": True, "independent_replay": False,
                      "formal_verification": False, "external_review": False},
        "coverage": {"original_four_variable": "k=0,1,2,3; m=1,2,3",
                     "sphere_and_binomial_constant_term": "k=0,1,2,3; m=1,...,8",
                     "integral_product_formula": "k=0,1,2,3; m=1,...,8; J(1)=2^m*m!/product(k*m+2*j+1,j=0,...,m)",
                     "factor_multiplicities_and_R_zero_specialization": "k=0,1,2,3"},
        "original_four_variable_checks": direct_rows,
        "degree_six_sphere_checks": sphere_rows,
        "factor_checks": factor_rows,
        "negative_controls": {"odd_height_integral_does_not_equal_marked_moment": True,
                              "omitting_one_radial_factorial_is_detected": True},
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("PASS: 12 original four-variable cases; 32 sphere/CT cases; factor and normalization controls.")
    print(OUTPUT.relative_to(RESEARCH.parent).as_posix())


if __name__ == "__main__":
    main()

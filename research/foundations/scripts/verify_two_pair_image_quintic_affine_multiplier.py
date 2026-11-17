#!/usr/bin/env python3
"""Exact finite certificate excluding affine multiples of the quartic witness.

For F=(aR+bT+cZ+dW)(R+Z)D of balanced degree five, compare the first four
original contractions with independently assembled sphere moments. Verify the
elimination and an explicit integer Bezout identity. These identities force
b=c=d=0 when those four pure moments vanish over a characteristic-zero field.
This covers the displayed affine-multiplier ansatz, not all primitive quintics.

Run from the repository root with Python 3.12 and SymPy 1.14.0:

    .venv/bin/python research/foundations/scripts/verify_two_pair_image_quintic_affine_multiplier.py
"""

from __future__ import annotations

import hashlib
import json
import platform
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

if not __debug__:
    raise SystemExit("Run without -O: verification assertions must be enabled.")

import sympy as sp

SCRIPT = Path(__file__).resolve()
RESEARCH = SCRIPT.parents[2]
OUTPUT = RESEARCH / "foundations/extended-geometry/two_pair_image_quintic_affine_multiplier_verification.json"
X, Y, T = sp.symbols("x y t")
A, B, C, D = sp.symbols("a b c d")
PARAMETERS = (A, B, C, D)
XI1, XI2, Z1, Z2 = sp.symbols("xi1 xi2 z1 z2")
VARS = (XI1, XI2, Z1, Z2)
CUTOFF = 4


@lru_cache(maxsize=None)
def sphere_monomial(i: int, j: int, k: int) -> sp.Rational:
    """Uniform phase and t in [-1,1], with xy=(1-t^2)/2."""
    if i != j or k % 2:
        return sp.S.Zero
    return sum(
        (sp.Rational((-1)**h * comb(i, h), 2**i * (k + 2*h + 1))
         for h in range(i + 1)),
        sp.S.Zero,
    )


def sphere_polynomial(order: int) -> sp.Expr:
    p0 = sp.Poly((1 + X)*(2*Y - (2 + X)*T**2), X, Y, T, domain=sp.QQ)
    base = p0**order
    result = sp.S.Zero
    for ia in range(order + 1):
        for ib in range(order - ia + 1):
            for ic in range(order - ia - ib + 1):
                id_ = order - ia - ib - ic
                value = sum(
                    (coefficient * sphere_monomial(i + ic, j + id_, k + ib)
                     for (i, j, k), coefficient in base.terms()),
                    sp.S.Zero,
                )
                if value:
                    multinomial = sp.Rational(
                        factorial(order),
                        factorial(ia)*factorial(ib)*factorial(ic)*factorial(id_),
                    )
                    result += multinomial*value*A**ia*B**ib*C**ic*D**id_
    return sp.factor(result)


def original_contraction(polynomial: sp.Poly) -> sp.Expr:
    """Literal E_2 on a balanced polynomial, retaining symbolic parameters."""
    degrees = {(i + j, k + l) for (i, j, k, l), _ in polynomial.terms()}
    assert len(degrees) == 1
    dual, coordinate = next(iter(degrees))
    assert dual == coordinate
    return sp.expand(sum(
        (coefficient*factorial(i)*factorial(j)
         for (i, j, k, l), coefficient in polynomial.terms()
         if i == k and j == l),
        sp.S.Zero,
    ))


def eliminate_even_b(expression: sp.Expr) -> sp.Expr:
    """Substitute c=3d/10 and b^2=-(130ad+27d^2)/25 exactly."""
    b_squared = -(130*A*D + 27*D**2)/25
    polynomial = sp.Poly(expression.subs(C, 3*D/10), B)
    assert all(exponent % 2 == 0 for (exponent,), _ in polynomial.terms())
    return sp.factor(sum(
        (coefficient*b_squared**(exponent//2)
         for (exponent,), coefficient in polynomial.terms()),
        sp.S.Zero,
    ))


def main() -> None:
    r = XI1*Z1 + XI2*Z2
    z = XI1*Z2
    w = 2*XI2*Z1
    t = XI1*Z1 - XI2*Z2
    cubic = 2*w*r**2 - 2*r*t**2 - z*t**2
    assert sp.expand(t**2 - r**2 + 2*z*w) == 0
    quartic = (r + z)*cubic
    expression = (A*r + B*t + C*z + D*w)*quartic
    original = sp.Poly(expression, *VARS, domain=sp.QQ.poly_ring(*PARAMETERS))

    moments = []
    rows = []
    for order in range(1, CUTOFF + 1):
        moment = sphere_polynomial(order)
        contraction = original_contraction(original**order)
        assert sp.expand(contraction - factorial(5*order + 1)*moment) == 0
        moments.append(moment)
        rows.append({
            "order": order,
            "radial_factor": f"({5*order + 1})!",
            "angular_moment_polynomial": str(moment),
            "original_four_variable_identity_passed": True,
        })

    assert sp.expand(moments[0] - (10*C - 3*D)/15) == 0
    second = 2*(56*A*C + 4*A*D + 4*B**2 + 28*C**2 - 4*C*D + 3*D**2)/105
    assert sp.expand(moments[1] - second) == 0
    reduced_second = 8*(130*A*D + 25*B**2 + 27*D**2)/2625
    assert sp.expand(moments[1].subs(C, 3*D/10) - reduced_second) == 0
    reduced_third = 8*D*(32175*A**2 + 56680*A*D + 14254*D**2)/625625
    reduced_fourth = 256*D*(60775*A**3 + 79220*A**2*D + 1802*A*D**2 - 643*D**3)/31906875
    assert sp.expand(eliminate_even_b(moments[2]) - reduced_third) == 0
    assert sp.expand(eliminate_even_b(moments[3]) - reduced_fourth) == 0

    u = sp.Symbol("u")
    p = 32175*u**2 + 56680*u + 14254
    q = 60775*u**3 + 79220*u**2 + 1802*u - 643
    bezout_a = -1531111624900*u**2 - 1247587548025*u + 290136844604
    bezout_b = 810588507300*u + 1031834055265
    bezout_constant = 3472141285450021
    assert sp.expand(bezout_a*p + bezout_b*q) == bezout_constant
    # Controls: radial padding passes; the naive t-multiple fails already at order 2.
    assert all(moment.subs({B: 0, C: 0, D: 0}) == 0 for moment in moments)
    assert moments[1].subs({A: 0, B: 1, C: 0, D: 0}) == sp.Rational(8, 105)

    receipt = {
        "status": "passed",
        "scope": "Exact symbolic first-four-moment and Bezout certificate for the specified four-parameter affine-multiplier ansatz. These finite identities force b=c=d=0 over characteristic zero. They do not classify arbitrary primitive quintics or establish independent review.",
        "input": {
            "R": str(r), "Z": str(z), "W": str(w), "T": str(t),
            "cubic_D": str(cubic),
            "family": "F=(aR+bT+cZ+dW)(R+Z)D",
            "sphere_quartic": "p0=(1+x)(2y-(2+x)t^2)",
            "sphere_relation": "t^2+2xy=1",
            "parameters": "a,b,c,d arbitrary over a characteristic-zero field",
            "cutoff": CUTOFF,
        },
        "arithmetic": "Exact rational polynomial arithmetic; symbolic parameters retained",
        "command": ".venv/bin/python research/foundations/scripts/verify_two_pair_image_quintic_affine_multiplier.py",
        "working_directory": "repository root",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "checker": {"path": SCRIPT.relative_to(RESEARCH.parent).as_posix(),
                    "sha256": hashlib.sha256(SCRIPT.read_bytes()).hexdigest()},
        "assurance": {"internal_replay": True, "independent_replay": False,
                      "formal_verification": False, "external_review": False},
        "moment_identity_checks": rows,
        "elimination": {
            "first_moment": "c=3d/10",
            "second_moment": "b^2=-(130ad+27d^2)/25",
            "reduced_third": str(reduced_third),
            "reduced_fourth": str(reduced_fourth),
            "d_nonzero_coordinate": "u=a/d",
            "P": str(p), "Q": str(q),
            "bezout_A": str(bezout_a), "bezout_B": str(bezout_b),
            "bezout_constant": bezout_constant,
            "bezout_identity_passed": True,
            "conclusion": "d=0, c=0, b=0; only F=aR(R+Z)D remains",
        },
        "controls": {"radial_padding_passes_first_four": True,
                     "naive_T_multiplier_second_angular_moment": "8/105"},
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("PASS: four symbolic original/sphere moment identities; elimination; integer Bezout certificate.")
    print(OUTPUT.relative_to(RESEARCH.parent).as_posix())


if __name__ == "__main__":
    main()

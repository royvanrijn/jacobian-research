#!/usr/bin/env python3
"""Exact degree-five bilinear-multiplier obstruction for the SIC2 seed."""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_degree_five_multiplier_obstruction.json"
)


def contraction(polynomial: sp.Poly) -> sp.Expr:
    result = 0
    for exponent, coefficient in polynomial.terms():
        xi1, xi2, zz1, zz2 = exponent
        if xi1 == zz1 and xi2 == zz2:
            result += coefficient * factorial(xi1) * factorial(xi2)
    return sp.factor(result)


def main() -> None:
    x1, x2, z1, z2 = sp.symbols("xi1 xi2 z1 z2")
    a, b, c, e, u = sp.symbols("a b c e u")
    variables = (x1, x2, z1, z2)

    r = x1 * z1 + x2 * z2
    z = x1 * z2
    w = 2 * x2 * z1
    t = x1 * z1 - x2 * z2
    f = sp.expand(
        (r + z) * (r**2 * w - sp.Rational(1, 2) * (2 * r + z) * t**2)
    )
    linear = a * r + b * z + c * w + e * t
    seed = sp.expand(linear * f)

    seed_polynomial = sp.Poly(seed, *variables)
    seed_power = sp.Poly(1, *variables)
    moments: dict[int, sp.Expr] = {}
    for order in range(1, 4):
        seed_power *= seed_polynomial
        moments[order] = contraction(seed_power)
    normalized = {
        1: 10 * b - 3 * c,
        2: (
            56 * a * b
            + 4 * a * c
            + 28 * b**2
            - 4 * b * c
            + 3 * c**2
            + 4 * e**2
        ),
        3: (
            3432 * a**2 * b
            + 6864 * a * b**2
            - 156 * a * c**2
            + 1144 * b**3
            + 1716 * b**2 * c
            + 130 * b * c**2
            - 49 * c**3
            - 260 * c * e**2
        ),
    }
    scales = {1: 24, 2: 190080, 3: 1045094400}
    for order in range(1, 4):
        assert sp.expand(moments[order] - scales[order] * normalized[order]) == 0

    centered = {
        order: sp.factor(
            normalized[order].subs(c, sp.Rational(10, 3) * b)
        )
        for order in (2, 3)
    }
    assert sp.expand(
        centered[2]
        - sp.Rational(4, 3)
        * (52 * a * b + 36 * b**2 + 3 * e**2)
    ) == 0
    assert sp.expand(
        centered[3]
        - sp.Rational(8, 27)
        * b
        * (
            11583 * a**2
            + 17316 * a * b
            + 21916 * b**2
            - 2925 * e**2
        )
    ) == 0

    nonradial_substitution = {
        b: 1,
        c: sp.Rational(10, 3),
        e: u,
        a: -(36 + 3 * u**2) / 52,
    }
    q = 8019 * u**4 - 623736 * u**2 + 3219760
    assert sp.expand(
        centered[3].subs(nonradial_substitution) - q / 702
    ) == 0

    specialized_seed = sp.Poly(
        seed.subs(nonradial_substitution),
        *variables,
    )
    specialized_power = sp.Poly(1, *variables)
    for _ in range(4):
        specialized_power *= specialized_seed
    fourth = contraction(specialized_power)
    p = (
        136323 * u**6
        - 5359284 * u**4
        - 174020976 * u**2
        - 802761152
    )
    fourth_scale = -sp.Rational(1235533824000, 169)
    assert sp.expand(fourth - fourth_scale * p) == 0

    gcd = sp.gcd(sp.Poly(p, u, domain=sp.QQ), sp.Poly(q, u, domain=sp.QQ))
    resultant = sp.resultant(p, q, u)
    expected_resultant = (
        97842802725670657299880334741299150646484103418606547061702656
    )
    assert gcd.as_expr() == 1
    assert resultant == expected_resultant

    artifact = {
        "format": "two-pair-degree-five-multiplier-obstruction-v1",
        "field": "characteristic zero",
        "family": "L*F with L=a*R+b*Z+c*W+e*T",
        "moments_used": [1, 2, 3, 4],
        "solution_set": "b=c=e=0; equivalently L=a*R",
        "nonradial_branch": {
            "q": str(q),
            "p": str(p),
            "gcd": "1",
            "resultant": str(resultant),
        },
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_DEGREE_FIVE_MULTIPLIER_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS degree-five multiplier obstruction: moments 1--4 force L=aR")
    print("PASS nonradial branch: exact gcd 1 and nonzero resultant")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

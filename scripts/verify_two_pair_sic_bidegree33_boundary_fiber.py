#!/usr/bin/env python3
"""Exact characteristic-zero finite-fiber certificate on the s0 boundary.

After the mu_2, A, and B constant pivots, the common boundary has six base
variables and two fiber variables s5,t4.  At one rational mu_3-zero base
point, mu_4 and mu_5 cut out a quotient of length six in the fiber.
"""

from __future__ import annotations

from math import factorial
import json
from pathlib import Path

import sympy as sp

from verify_two_pair_sic_bidegree33_anchor_jacobians import (
    COEFFICIENT_MAP,
    NORMALIZED_QUADRATIC,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_boundary_fiber.json"
)
BASE_POINT = {
    "s1": -3,
    "s2": -3,
    "s3": sp.Rational(45973, 1026),
    "t0": -3,
    "t1": 3,
    "t2": -3,
}


def main() -> None:
    s5, t4, x, y = sp.symbols("s5 t4 x y")
    values: list[sp.Expr] = [sp.Integer(0)] * 12
    values[0] = 1
    values[1] = BASE_POINT["s1"]
    values[2] = BASE_POINT["s2"]
    values[3] = BASE_POINT["s3"]
    values[5] = s5
    values[7] = BASE_POINT["t0"]
    values[8] = BASE_POINT["t1"]
    values[9] = BASE_POINT["t2"]
    values[11] = t4

    s1, s2, s3 = values[1], values[2], values[3]
    t0, t1, t2 = values[7], values[8], values[9]
    # A=0 and B=0, using their constant t3 and s4 pivots.
    values[10] = (
        -6 * s1**2 * t1
        + 3 * s1 * s2 * t0
        + 3 * s1 * t2
        + 3 * s2 * t1
        - 2 * s3 * t0
        + 3 * t0
    )
    values[4] = sp.Rational(1, 3) * (
        12 * s1 * s3
        + 28 * s1 * t0 * t1
        - 18 * s1
        - 9 * s2**2
        - 14 * s2 * t0**2
        - 2 * t0 * t2
        - 12 * t1**2
    )
    # mu_2=0, using the constant s6 pivot on s0=1.
    values[6] = sp.Rational(1, 3) * (
        18 * s1 * s5
        - 45 * s2 * values[4]
        + 30 * s3**2
        + 14 * t0 * t4
        - 56 * t1 * values[10]
        + 42 * t2**2
        + 70
    )

    polynomial = sp.Integer(0)
    for i, j, parameter, coefficient in COEFFICIENT_MAP:
        polynomial += coefficient * values[parameter] * x**i * y**j
    for i, j, coefficient in NORMALIZED_QUADRATIC:
        polynomial += coefficient * x**i * y**j

    moments: dict[int, sp.Poly] = {}
    for order in range(2, 6):
        power = sp.Poly(sp.expand(polynomial**order), x, y)
        moment = sum(
            factorial(3 * order - diagonal)
            * factorial(diagonal)
            * power.coeff_monomial(x**diagonal * y**diagonal)
            for diagonal in range(3 * order + 1)
        )
        moments[order] = sp.Poly(moment, s5, t4)

    assert moments[2].is_zero
    assert moments[3].is_zero
    assert moments[4].total_degree() == 2 and len(moments[4].terms()) == 6
    assert moments[5].total_degree() == 3 and len(moments[5].terms()) == 8
    assert sp.gcd(moments[4], moments[5]).total_degree() == 0

    basis = sp.groebner(
        [moments[4].as_expr(), moments[5].as_expr()],
        s5,
        t4,
        order="grevlex",
    )
    leading_exponents = [
        tuple(polynomial.LM(order=basis.order).exponents)
        for polynomial in basis.polys
    ]
    assert leading_exponents == [(0, 4), (1, 2), (2, 0)]
    standard_monomials = [
        (left, right)
        for left in range(5)
        for right in range(5)
        if not any(
            left >= lead_left and right >= lead_right
            for lead_left, lead_right in leading_exponents
        )
    ]
    assert standard_monomials == [
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
    ]

    payload = {
        "field": "QQ",
        "base_point": {key: str(value) for key, value in BASE_POINT.items()},
        "vanishing_moments": [2, 3],
        "fiber_variables": ["s5", "t4"],
        "mu4": str(moments[4].as_expr()),
        "mu5": str(moments[5].as_expr()),
        "groebner_basis": [str(polynomial.as_expr()) for polynomial in basis.polys],
        "leading_exponents": leading_exponents,
        "standard_monomials": standard_monomials,
        "quotient_length": len(standard_monomials),
        "scope": (
            "exact characteristic-zero rank-six fiber certificate; "
            "not a full common-boundary zero-fiber exclusion"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("PASS s0 common boundary: mu_2=mu_3=0 at the recorded base point")
    print("PASS fiber quotient: (mu_4,mu_5) has length 6 over QQ")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

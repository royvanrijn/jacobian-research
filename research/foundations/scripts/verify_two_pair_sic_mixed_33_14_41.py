#!/usr/bin/env python3
"""Exact nonzero-(1,4) theorem on the (3,3)+(1,4)+(4,1) stratum.

Let F=A+B+C have the displayed bidegree blocks.  On B != 0, the
dual-linear theorem normalizes B=xi2*z1^4.  This checker verifies:

1. the maximal positive-weight contractions and their linear solution;
2. two further polynomial-valued coefficients forcing A upper triangular;
3. the four-variable central core through moment five;
4. exact power memberships proving that core has radical the origin; and
5. the residual two-step exponent cone giving eventual mixed vanishing.

The B=0 boundary reduces to the open balanced bidegree-(3,3) problem and
is deliberately not claimed here.
"""

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
    / "two_pair_sic_mixed_33_14_41.json"
)

XI1, XI2, Z1, Z2 = sp.symbols("xi1 xi2 z1 z2")
VARIABLES = (XI1, XI2, Z1, Z2)


def contraction(expression: sp.Expr) -> sp.Expr:
    result = sp.Integer(0)
    for (xi1, xi2, z1, z2), coefficient in sp.Poly(
        sp.expand(expression), *VARIABLES
    ).terms():
        if xi1 > z1 or xi2 > z2:
            continue
        result += (
            coefficient
            * sp.Rational(factorial(z1), factorial(z1 - xi1))
            * sp.Rational(factorial(z2), factorial(z2 - xi2))
            * Z1 ** (z1 - xi1)
            * Z2 ** (z2 - xi2)
        )
    return sp.expand(result)


def output_coefficient(
    expression: sp.Expr, z1_degree: int, z2_degree: int
) -> sp.Expr:
    return sp.Poly(sp.expand(expression), Z1, Z2).coeff_monomial(
        Z1**z1_degree * Z2**z2_degree
    )


def exponent_data(
    exponent: tuple[int, int, int, int],
) -> dict[str, int | list[int]]:
    xi1, xi2, z1, z2 = exponent
    weight = z1 + z2 - xi1 - xi2
    delta = z1 - xi1
    epsilon = xi2 - z2
    assert epsilon == delta - weight
    return {
        "exponent": list(exponent),
        "weight": weight,
        "delta": delta,
        "epsilon": epsilon,
    }


def main() -> None:
    a = sp.symbols("a0:16")
    c = sp.symbols("c0:10")

    A = sp.Integer(0)
    for i in range(4):
        for j in range(4):
            A += (
                a[4 * i + j]
                * XI1**i
                * XI2 ** (3 - i)
                * Z1**j
                * Z2 ** (3 - j)
            )

    B = XI2 * Z1**4

    C = sp.Integer(0)
    for i in range(5):
        for j in range(2):
            C += (
                c[2 * i + j]
                * XI1**i
                * XI2 ** (4 - i)
                * Z1**j
                * Z2 ** (1 - j)
            )

    highest = {
        order: contraction(A * B**order)
        for order in range(4)
    }
    assert highest[0] == 2 * (3 * a[0] + a[5] + a[10] + 3 * a[15])
    assert sp.expand(
        highest[1]
        - (24 * a[4] + 40 * a[9] + 120 * a[14]) * Z1**3
        - (72 * a[8] + 120 * a[13]) * Z1**2 * Z2
        - 72 * a[12] * Z1 * Z2**2
    ) == 0
    assert sp.expand(
        highest[2]
        - (336 * a[8] + 1008 * a[13]) * Z1**6
        - 2016 * a[12] * Z1**5 * Z2
    ) == 0
    assert highest[3] == 7920 * a[12] * Z1**9

    highest_substitution = {
        a[12]: 0,
        a[13]: 0,
        a[8]: 0,
        a[4]: -5 * a[14] - sp.Rational(5, 3) * a[9],
        a[0]: (
            -sp.Rational(1, 3) * a[10]
            - a[15]
            - sp.Rational(1, 3) * a[5]
        ),
    }
    reduced_A = sp.expand(A.subs(highest_substitution))

    # Weight three in F^3 and weight six in F^4.
    third_weight_three = contraction(
        3 * reduced_A**2 * B + 3 * B**2 * C
    )
    fourth_weight_six = contraction(
        6 * reduced_A**2 * B**2 + 4 * B**3 * C
    )
    r, s = a[14], a[9]
    q1 = 118 * r**2 + 51 * r * s + 6 * s**2
    q2 = 72 * r**2 + 21 * r * s + 2 * s**2
    assert sp.expand(
        output_coefficient(third_weight_three, 2, 1) - 4320 * q1
    ) == 0
    assert sp.expand(
        output_coefficient(fourth_weight_six, 6, 0) - 94080 * q2
    ) == 0

    triangular_basis = sp.groebner((q1, q2), r, s, order="lex")
    expected_triangular_basis = (
        199 * r**2 - 4 * s**2,
        597 * r * s + 98 * s**2,
        s**3,
    )
    assert tuple(triangular_basis.polys) == tuple(
        sp.Poly(item, r, s) for item in expected_triangular_basis
    )

    b, d, e, h = sp.symbols("b d e h")
    leading_diagonal = -(b + d + 3 * e) / 3
    diagonal_A = (
        leading_diagonal * (XI2 * Z2) ** 3
        + b * (XI1 * Z1) * (XI2 * Z2) ** 2
        + d * (XI1 * Z1) ** 2 * (XI2 * Z2)
        + e * (XI1 * Z1) ** 3
    )
    core = diagonal_A + B + h * XI1**4 * Z2

    normalized_moments: dict[int, sp.Expr] = {}
    contents: dict[int, sp.Rational] = {}
    for order in range(1, 6):
        value = contraction(core**order)
        assert not (value.free_symbols & {Z1, Z2})
        if order == 1:
            assert value == 0
            continue
        content, primitive = sp.Poly(value, b, d, e, h).primitive()
        contents[order] = content
        normalized_moments[order] = primitive.as_expr()

    assert contents == {
        2: 24,
        3: 288,
        4: 23040,
        5: 1612800,
    }
    f2 = normalized_moments[2]
    f3 = normalized_moments[3]
    assert f2 == (
        2 * b**2
        + 5 * b * d
        + 13 * b * e
        + 4 * d**2
        + 25 * d * e
        + 57 * e**2
        + 2 * h
    )
    assert f3 == (
        -20 * b**3
        - 70 * b**2 * d
        - 230 * b**2 * e
        - 80 * b * d**2
        - 520 * b * d * e
        - 780 * b * e**2
        + 11 * b * h
        - 30 * d**3
        - 290 * d**2 * e
        - 780 * d * e**2
        + 26 * d * h
        + 93 * e * h
    )

    core_basis = sp.groebner(
        tuple(normalized_moments.values()),
        b,
        d,
        e,
        h,
        order="grevlex",
    )
    assert core_basis.is_zero_dimensional
    power_certificates = {b: 10, d: 10, e: 10, h: 5}
    for variable, power in power_certificates.items():
        assert core_basis.reduce(variable**power)[1] == 0

    # Audit the residual cone after the diagonal core and c_(4,0) vanish.
    residual_exponents: dict[str, tuple[int, int, int, int]] = {
        "B": (0, 1, 4, 0),
    }
    for i in range(4):
        for j in range(i + 1, 4):
            residual_exponents[f"a{i}{j}"] = (
                i,
                3 - i,
                j,
                3 - j,
            )
    for i in range(5):
        for j in range(2):
            if (i, j) == (4, 0):
                continue
            residual_exponents[f"c{i}{j}"] = (
                i,
                4 - i,
                j,
                1 - j,
            )

    cone = {
        name: exponent_data(exponent)
        for name, exponent in residual_exponents.items()
    }
    for item in cone.values():
        epsilon = int(item["epsilon"])
        delta = int(item["delta"])
        assert epsilon >= 0
        if epsilon == 0:
            assert delta == -3
        else:
            assert delta <= 4

    payload = {
        "claim": (
            "On the nonzero V_(1,4) branch, SIC(2) holds on "
            "V_(3,3)+V_(1,4)+V_(4,1); the first five pure "
            "contractions suffice"
        ),
        "claim_boundary": (
            "The V_(1,4)=0 boundary reduces to the open balanced "
            "bidegree-(3,3) problem"
        ),
        "normalized_positive_piece": "B=xi2*z1^4",
        "highest_weight_contractions": {
            str(order): str(sp.factor(value))
            for order, value in highest.items()
        },
        "triangular_equations": {
            "q1": str(q1),
            "q2": str(q2),
            "groebner_basis": [
                str(item.as_expr()) for item in triangular_basis.polys
            ],
            "conclusion": "a14=a9=a4=0",
        },
        "central_core": {
            "variables": ["b", "d", "e", "h"],
            "normalized_moments": {
                str(order): str(polynomial)
                for order, polynomial in normalized_moments.items()
            },
            "groebner_basis_size": len(core_basis),
            "groebner_basis": [
                str(item.as_expr()) for item in core_basis.polys
            ],
            "power_certificates": {
                str(variable): power
                for variable, power in power_certificates.items()
            },
            "radical": "(b,d,e,h)",
        },
        "residual_cone": cone,
        "multiplier_cutoff_rule": (
            "For a multiplier monomial Q with epsilon_Q=xi2-z2 and "
            "delta_Q=z1-xi1, put R=max(0,-epsilon_Q).  No monomial of "
            "E_2(Q*F^m) survives when m>(delta_Q+7*R)/3."
        ),
        "status": "exact characteristic-zero theorem on B!=0",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    print(
        "PASS mixed (3,3)+(1,4)+(4,1), nonzero (1,4) branch: "
        "triangular reduction, central radical, and multiplier cone"
    )


if __name__ == "__main__":
    main()

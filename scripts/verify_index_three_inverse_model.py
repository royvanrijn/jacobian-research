#!/usr/bin/env python3
"""Exact lower-degree calibration for the cubic index-three attack.

The triangular cubic map

    F = (x0 + x1^3, x1 + x2^3, x2)

has ``(JH)^3=0`` but ``(JH)^2!=0``.  Its inverse has a nonzero degree-nine
term, so no general index-three inverse theorem can truncate at degree seven.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "index_three_inverse_degree_model.json"
)


def encode(expressions: list[sp.Expr], variables: tuple[sp.Symbol, ...]):
    encoded = []
    for expression in expressions:
        poly = sp.Poly(expression, *variables, domain=sp.QQ)
        encoded.append(
            [
                {
                    "monomial": [
                        [index, power]
                        for index, power in enumerate(exponents)
                        if power
                    ],
                    "coefficient": str(coefficient),
                }
                for exponents, coefficient in poly.terms()
            ]
        )
    return encoded


def homogeneous_part(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    degree: int,
) -> sp.Expr:
    poly = sp.Poly(expression, *variables, domain=sp.QQ)
    return sp.Add(
        *[
            coefficient
            * sp.prod(
                variable**power
                for variable, power in zip(variables, exponents)
            )
            for exponents, coefficient in poly.terms()
            if sum(exponents) == degree
        ]
    )


def main() -> None:
    x = sp.symbols("x0:3")
    y = sp.symbols("y0:3")
    h = [x[1] ** 3, x[2] ** 3, sp.Integer(0)]
    f = [x[index] + h[index] for index in range(3)]
    jacobian = sp.Matrix(h).jacobian(x)
    assert jacobian**2 != sp.zeros(3)
    assert jacobian**3 == sp.zeros(3)
    assert (sp.eye(3) + jacobian).det() == 1

    inverse = [
        y[0] - (y[1] - y[2] ** 3) ** 3,
        y[1] - y[2] ** 3,
        y[2],
    ]
    inverse_substitution = dict(zip(x, inverse))
    forward_substitution = dict(zip(y, f))
    assert [
        sp.expand(expression.subs(inverse_substitution, simultaneous=True))
        for expression in f
    ] == list(y)
    assert [
        sp.expand(expression.subs(forward_substitution, simultaneous=True))
        for expression in inverse
    ] == list(x)

    correction = [
        sp.expand(inverse[index] - y[index]) for index in range(3)
    ]
    degrees = (3, 5, 7, 9)
    parts = {
        degree: [
            sp.expand(homogeneous_part(expression, y, degree))
            for expression in correction
        ]
        for degree in degrees
    }
    assert all(any(component != 0 for component in parts[d]) for d in degrees)
    assert parts[9][0] == y[2] ** 9
    assert max(
        sp.Poly(expression, *y, domain=sp.QQ).total_degree()
        for expression in inverse
    ) == 9

    payload = {
        "format": "cubic-index-three-inverse-degree-model-v1",
        "field": "QQ",
        "dimension": 3,
        "H": encode(h, x),
        "nilpotency_index_JH": 3,
        "jacobian_determinant": "1",
        "inverse": encode(inverse, y),
        "inverse_degree": 9,
        "nonzero_inverse_correction_degrees": list(degrees),
        "degree_nine_term": "y2^9 in the first inverse component",
        "consequence": (
            "a uniform cubic index-three inverse-degree bound, if true, "
            "cannot be smaller than 9"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS index-three model: (JH)^2!=0 and (JH)^3=0")
    print("PASS index-three model: exact two-sided polynomial inverse")
    print("LOWER CALIBRATION index-three model: inverse degree is 9")
    print(f"PASS index-three model: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

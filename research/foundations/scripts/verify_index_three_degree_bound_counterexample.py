#!/usr/bin/env python3
"""Verify van den Essen's cubic index-three inverse-degree counterexample.

The cubic homogeneous map ``F=X+H`` in dimension five with

    H = (
        3*x4^2*x2 - 2*x4*x5*x3,
        x4^2*x5,
        x4^3,
        x5^3,
        0,
    )

has weak Jacobian nilpotency index three and generic Jacobian rank three.
It is nevertheless a polynomial automorphism whose inverse has degree
thirteen.  Thus the Druzkowski--Rusek bound ``deg(F^-1)<=3^(nu-1)``, and in
particular the proposed degree-nine full-class index-three bound, is false.

The example is from A. van den Essen, "A counterexample to a conjecture of
Druzkowski and Rusek", Ann. Polon. Math. 62 (1995), 173--176.
"""

from __future__ import annotations

from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "index_three_degree_bound_counterexample.json"
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


def derivative_tensor(
    h: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    children: list[sp.Matrix],
) -> sp.Matrix:
    answer = []
    for component in h:
        value = sp.Integer(0)
        for indices in product(range(len(variables)), repeat=len(children)):
            derivative = component
            for index in indices:
                derivative = sp.diff(derivative, variables[index])
            value += derivative * sp.prod(
                children[slot][index]
                for slot, index in enumerate(indices)
            )
        answer.append(sp.expand(value))
    return sp.Matrix(answer)


def main() -> None:
    x = sp.symbols("x1:6")
    x1, x2, x3, x4, x5 = x
    h = sp.Matrix(
        [
            3 * x4**2 * x2 - 2 * x4 * x5 * x3,
            x4**2 * x5,
            x4**3,
            x5**3,
            0,
        ]
    )
    jacobian = h.jacobian(x)
    assert jacobian**2 != sp.zeros(5)
    assert jacobian**3 == sp.zeros(5)
    assert (sp.eye(5) + jacobian).det() == 1

    nonzero_three_minor = sp.factor(
        jacobian.extract((0, 1, 2), (1, 3, 4)).det()
    )
    assert nonzero_three_minor == -9 * x4**6
    assert all(
        sp.expand(jacobian.extract(rows, columns).det()) == 0
        for rows in combinations(range(5), 4)
        for columns in combinations(range(5), 4)
    )

    # Weak index three is strictly smaller than the strong index four.
    points = (
        (0, 0, 0, 1, 0),
        (0, 1, 0, 1, 1),
        (0, 0, 0, 0, 1),
    )
    mixed_cube = sp.eye(5)
    for point in points:
        mixed_cube *= jacobian.subs(dict(zip(x, point)), simultaneous=True)
    assert mixed_cube[0, 4] == 18
    assert sum(bool(entry) for entry in mixed_cube) == 1
    a = sp.symbols("a1:6")
    b = sp.symbols("b1:6")
    c = sp.symbols("c1:6")
    d = sp.symbols("d1:6")
    independent_jacobians = [
        jacobian.subs(dict(zip(x, point)), simultaneous=True)
        for point in (a, b, c, d)
    ]
    assert sp.simplify(sp.prod(independent_jacobians)) == sp.zeros(5)

    y = sp.symbols("y1:6")
    y1, y2, y3, y4, y5 = y
    u = y4 - y5**3
    inverse = sp.Matrix(
        [
            y1
            - 3 * u**2 * (y2 - u**2 * y5)
            + 2 * u * y5 * (y3 - u**3),
            y2 - u**2 * y5,
            y3 - u**3,
            y4 - y5**3,
            y5,
        ]
    )
    forward = sp.Matrix(x) + h
    assert [
        sp.expand(expression.subs(dict(zip(x, inverse)), simultaneous=True))
        for expression in forward
    ] == list(y)
    assert [
        sp.expand(expression.subs(dict(zip(y, forward)), simultaneous=True))
        for expression in inverse
    ] == list(x)

    correction = [
        sp.expand(inverse[index] - y[index]) for index in range(5)
    ]
    degrees = (3, 5, 7, 9, 11, 13)
    parts = {
        degree: [
            sp.expand(homogeneous_part(expression, y, degree))
            for expression in correction
        ]
        for degree in degrees
    }
    assert all(any(component != 0 for component in parts[d]) for d in degrees)
    assert parts[11] == [-4 * y4 * y5**10, 0, 0, 0, 0]
    assert parts[13] == [y5**13, 0, 0, 0, 0]
    inverse_degree = max(sp.Poly(expression, *y).total_degree() for expression in inverse)
    assert inverse_degree == 13

    # Evaluate the repository's degree-eleven quotient representative
    # Omega_11=-C(NH,H,H)/2-B(B(H,H),H)/2-N(C(H,H,H))/6.
    nh = jacobian * h
    bhh = derivative_tensor(h, x, [h, h])
    chhh = derivative_tensor(h, x, [h, h, h])
    omega_11 = (
        -sp.Rational(1, 2) * derivative_tensor(h, x, [nh, h, h])
        - sp.Rational(1, 2) * derivative_tensor(h, x, [bhh, h])
        - sp.Rational(1, 6) * jacobian * chhh
    )
    omega_11 = sp.Matrix([sp.expand(component) for component in omega_11])
    assert omega_11 == sp.Matrix([-4 * x4 * x5**10, 0, 0, 0, 0])

    payload = {
        "format": "cubic-index-three-degree-bound-counterexample-v1",
        "field": "QQ",
        "source": {
            "author": "Arno van den Essen",
            "title": "A counterexample to a conjecture of Druzkowski and Rusek",
            "journal": "Annales Polonici Mathematici 62.2 (1995), 173-176",
            "url": "https://eudml.org/doc/262620",
        },
        "dimension": 5,
        "H": encode(list(h), x),
        "generic_rank_JH": 3,
        "rank_three_minor": "-9*x4^6",
        "weak_nilpotency_index_JH": 3,
        "strong_nilpotency_index_JH": 4,
        "jacobian_determinant": "1",
        "inverse": encode(list(inverse), y),
        "inverse_degree": inverse_degree,
        "nonzero_inverse_correction_degrees": list(degrees),
        "degree_eleven_term": "-4*y4*y5^10 in the first inverse component",
        "degree_thirteen_term": "y5^13 in the first inverse component",
        "omega_11_evaluation": "(-4*x4*x5^10,0,0,0,0)",
        "consequences": [
            "the proposed uniform inverse-degree bound 9 is false",
            "the degree-11 inverse covariant survives the full coefficient variety (JH)^3=0",
            "the surviving tensor is an automorphism and does not refute index-three invertibility",
            "the first degree-bound failure occurs in dimension 5 at generic rank 3",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS van den Essen index-three map: generic rank 3 and weak index 3")
    print("PASS van den Essen index-three map: exact two-sided inverse has degree 13")
    print("PASS van den Essen index-three map: Omega_11 is nonzero")
    print(f"PASS van den Essen index-three map: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

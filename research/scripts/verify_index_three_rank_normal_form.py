#!/usr/bin/env python3
"""Eliminate index three on a nontriangular cubic rank-three normal form.

Consider

    H=(0, lambda*x1^3, x2*Delta+p(x1,x2), x1*Delta+q(x1,x2)),
    Delta=x1*x3-x2*x4,

with arbitrary binary cubics ``p,q``.  This family is a standard
four-variable homogeneous nilpotent-Jacobian normal-form ansatz.  Exact
coefficient extraction shows that ``(JH)^3=0`` forces ``lambda=0``.

On that locus only two components remain, ``rank(JH)<=2``, and the map has
an explicit inverse of degree at most five.  Thus the rank-three branch of
this entire family has nilpotency index at least four and cannot meet either
open endpoint ``rank=3,index=3``.
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
    / "index_three_rank_normal_form_exclusion.json"
)


def main() -> None:
    x1, x2, x3, x4 = x = sp.symbols("x1:5")
    lam = sp.symbols("lambda")
    a = sp.symbols("a0:4")
    b = sp.symbols("b0:4")
    p = sum(a[index] * x1 ** (3 - index) * x2**index for index in range(4))
    q = sum(b[index] * x1 ** (3 - index) * x2**index for index in range(4))
    delta = x1 * x3 - x2 * x4
    h = sp.Matrix([0, lam * x1**3, x2 * delta + p, x1 * delta + q])
    jacobian = h.jacobian(x)
    cube = jacobian**3

    coefficient_equations = sorted(
        {
            sp.factor(coefficient)
            for entry in cube
            for coefficient in sp.Poly(sp.expand(entry), *x).coeffs()
            if coefficient
        },
        key=str,
    )
    assert sp.Integer(3) * lam in coefficient_equations
    assert -sp.Integer(3) * lam in coefficient_equations
    assert all(sp.rem(equation, lam, lam) == 0 for equation in coefficient_equations)

    h0 = sp.expand(h.subs(lam, 0))
    jacobian0 = h0.jacobian(x)
    assert jacobian0**3 == sp.zeros(4)
    assert all(
        minor == 0
        for minor in (
            sp.Matrix(jacobian0).extract(rows, columns).det()
            for rows in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
            for columns in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
        )
    )

    # With lambda=0, the last two coordinates have the affine-linear block
    # z -> (I+B)z+c over QQ[x1,x2], and B^2=0.
    block = sp.Matrix([[x1 * x2, -x2**2], [x1**2, -x1 * x2]])
    constant = sp.Matrix([p, q])
    assert block**2 == sp.zeros(2)
    y1, y2, y3, y4 = y = sp.symbols("y1:5")
    block_y = block.subs({x1: y1, x2: y2})
    constant_y = constant.subs({x1: y1, x2: y2})
    inverse_tail = (sp.eye(2) - block_y) * (
        sp.Matrix([y3, y4]) - constant_y
    )
    inverse = sp.Matrix([y1, y2, *map(sp.expand, inverse_tail)])

    forward0 = sp.Matrix(x) + h0
    assert [
        sp.expand(expression.subs(dict(zip(x, inverse)), simultaneous=True))
        for expression in forward0
    ] == list(y)
    assert [
        sp.expand(expression.subs(dict(zip(y, forward0)), simultaneous=True))
        for expression in inverse
    ] == list(x)
    inverse_degree = max(
        sp.Poly(expression, *y, domain=sp.QQ.frac_field(*a, *b)).total_degree()
        for expression in inverse
    )
    assert inverse_degree <= 5

    payload = {
        "format": "cubic-index-three-rank-normal-form-exclusion-v1",
        "field": "characteristic zero",
        "family": (
            "H=(0,lambda*x1^3,x2*(x1*x3-x2*x4)+p,"
            "x1*(x1*x3-x2*x4)+q), p,q arbitrary binary cubics"
        ),
        "coefficient_equations_of_JH_cube": [
            str(equation) for equation in coefficient_equations
        ],
        "index_three_condition": "lambda=0",
        "index_three_locus_rank_bound": 2,
        "index_three_locus_inverse_degree_bound": 5,
        "rank_three_consequence": (
            "the lambda!=0 branch has (JH)^3!=0, hence index at least 4"
        ),
        "frontier_consequence": (
            "the entire displayed normal-form family contains no map simultaneously "
            "on the rank-three and index-three open strata"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS rank/index normal form: (JH)^3=0 forces lambda=0")
    print("PASS rank/index normal form: index-three locus has rank at most 2")
    print("PASS rank/index normal form: explicit inverse has degree at most 5")
    print(f"PASS rank/index normal form: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

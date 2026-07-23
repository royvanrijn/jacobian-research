#!/usr/bin/env python3
"""Verify the cubic remainder synchronization identity for degree thirty.

After centering the inner and outer cubics, the open pair ``{3,10}``
amounts to

    G(y)^3 + u G(y) = C(D(y)) + delta*y + constant,
    D(y) = y^3 + p*y + q,

with ``G`` monic of degree ten.  Write the unique remainder of ``G`` over
``QQ[D]`` as

    G = g0(D) + y*g1(D) + y^2*g2(D).

The Hessian equations say that the ``y^2`` remainder of ``G^3+uG`` is zero
and that every positive ``D``-coefficient of its ``y`` remainder is zero.
Synchronization is the assertion that the remaining constant coefficient
of the ``y`` remainder is already in that ideal.

This checker asks Singular for that exact characteristic-zero ideal
membership in thirteen variables.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


y, z = sp.symbols("y z")
p, q, u = sp.symbols("p q u")
a = sp.symbols("a0:4")
b = sp.symbols("b0:3")
c = sp.symbols("c0:3")
VARIABLES = a + b + c + (p, q, u)


def serialize_singular(expression: sp.Expr) -> str:
    """Serialize a rational polynomial for Singular."""

    return str(sp.expand(expression)).replace("**", "^")


def remainder_data() -> tuple[list[sp.Expr], sp.Expr]:
    """Return the Hessian remainder ideal and its linear defect."""

    g0 = sum(a[index] * z**index for index in range(4))
    g1 = z**3 + sum(b[index] * z**index for index in range(3))
    g2 = sum(c[index] * z**index for index in range(3))
    g = g0 + y * g1 + y**2 * g2
    cubic = y**3 + p * y + q - z
    remainder = sp.rem(
        sp.Poly(sp.expand(g**3 + u * g), y),
        sp.Poly(cubic, y),
    )
    coefficient_y = sp.Poly(remainder, y).nth(1)
    coefficient_y2 = sp.Poly(remainder, y).nth(2)
    y_polynomial = sp.Poly(coefficient_y, z)
    y2_polynomial = sp.Poly(coefficient_y2, z)
    equations = [
        y2_polynomial.nth(power)
        for power in range(y2_polynomial.degree() + 1)
    ]
    equations.extend(
        y_polynomial.nth(power)
        for power in range(1, y_polynomial.degree() + 1)
    )
    equations = [sp.primitive(sp.Poly(item, *VARIABLES))[1].as_expr()
                 for item in equations if item != 0]
    defect = y_polynomial.nth(0)
    return equations, defect


def main() -> None:
    equations, defect = remainder_data()
    singular = shutil.which("Singular")
    assert singular is not None
    program = (
        f'ring r=0,({",".join(map(str, VARIABLES))}),dp;\n'
        f'ideal I={",".join(serialize_singular(item) for item in equations)};\n'
        "ideal G=slimgb(I);\n"
        'print("CUBIC_REMAINDER_SYNC");\n'
        f"print(reduce({serialize_singular(defect)},G)==0);\n"
        "print(size(G));\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=300,
    )
    compact = result.stdout.split()
    marker = compact.index("CUBIC_REMAINDER_SYNC")
    assert compact[marker + 1] == "1", result.stdout + result.stderr
    print(
        "PASS: the degree-30 cubic remainder defect is in the exact "
        f"Hessian remainder ideal; basis size {compact[marker + 2]}"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify the conic-divisible gate in the clean double-conic packet.

Normalize the smooth conic to q=x*z-y^2.  Its projective automorphism group
has two orbits on residual lines, represented by z (tangent) and y (secant).
For a general quintic h=q*G3, exact rational saturation proves that

    det(Hess(h)) = k*q^4*line,  k != 0,

has no solution in either orbit.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


x, y, z = sp.symbols("x y z")
q = x * z - y**2
cubic_monomials = [
    x**i * y**j * z ** (3 - i - j)
    for i in range(4)
    for j in range(4 - i)
]
coefficients = sp.symbols("g0:10")
k = sp.symbols("k")
G3 = sum(coefficient * monomial for coefficient, monomial in zip(coefficients, cubic_monomials))
h5 = sp.expand(q * G3)
determinant = sp.expand(sp.hessian(h5, (x, y, z)).det())
ring_variables = coefficients + (k,)


def singular_polynomial(expression: sp.Expr) -> str:
    _, polynomial = sp.Poly(
        expression, *ring_variables, domain=sp.QQ
    ).clear_denoms(convert=True)
    return str(sp.expand(polynomial.as_expr())).replace("**", "^")


def coefficient_equations(expression: sp.Expr) -> list[sp.Expr]:
    polynomial = sp.Poly(sp.expand(expression), x, y, z)
    return [
        polynomial.coeff_monomial(x**i * y**j * z ** (9 - i - j))
        for i in range(10)
        for j in range(10 - i)
        if polynomial.coeff_monomial(x**i * y**j * z ** (9 - i - j)) != 0
    ]


def verify_orbit(name: str, residual_line: sp.Expr) -> None:
    equations = coefficient_equations(determinant - k * q**4 * residual_line)
    source = "\n".join(
        [
            'LIB "elim.lib";',
            f"ring r=0,({','.join(map(str, ring_variables))}),dp;",
            "option(redSB);",
            f"ideal I={','.join(singular_polynomial(equation) for equation in equations)};",
            "ideal K=k;",
            "ideal J=sat(I,K);",
            'print("NONZERO_DETERMINANT_SATURATION");',
            "J;",
            "exit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=source,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "NONZERO_DETERMINANT_SATURATION\nJ[1]=1" in completed.stdout
    print(f"PASS: {name} residual-line orbit has unit nonzero-determinant saturation")


def main() -> None:
    assert shutil.which("Singular") is not None, "Singular is required"
    verify_orbit("tangent", z)
    verify_orbit("secant", y)

    # Calibration: the deepest obvious family h=q^2*L has only q^3 in its
    # Hessian, and its residual quadratic is never another q.
    a, b, c = sp.symbols("a b c")
    linear = a * x + b * y + c * z
    deep_determinant = sp.factor(sp.hessian(q**2 * linear, (x, y, z)).det())
    residual_quadratic = (
        6 * a**2 * x**2
        + 12 * a * b * x * y
        + (8 * a * c + b**2) * x * z
        + (4 * a * c + 5 * b**2) * y**2
        + 12 * b * c * y * z
        + 6 * c**2 * z**2
    )
    assert sp.expand(deep_determinant - 16 * q**3 * linear * residual_quadratic) == 0
    mu = sp.symbols("mu")
    divisibility_equations = [
        coefficient
        for _, coefficient in sp.Poly(
            residual_quadratic - mu * q, x, y, z
        ).terms()
    ]
    divisibility_basis = sp.groebner(
        divisibility_equations, a, b, c, mu, order="grevlex"
    )
    assert divisibility_basis.reduce(a**2)[1] == 0
    assert divisibility_basis.reduce(b**3)[1] == 0
    assert divisibility_basis.reduce(c**2)[1] == 0

    print("PASS: the q^2*L calibration has exact conic multiplicity three")
    print("THEOREM: a clean double-conic survivor cannot have h5 divisible by q")


if __name__ == "__main__":
    main()

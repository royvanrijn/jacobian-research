#!/usr/bin/env python3
"""Test the degree-24 all-factor-order intersection on a prime-word chart.

The chart starts with a generic complete decomposition of degree word
``3 o 2 o 2 o 2``.  It therefore lies on C_(3,8), C_(6,4), and C_(12,2).
We pull back the canonical residuals for the three opposite cuts C_(2,12),
C_(4,6), and C_(8,3), then ask Singular for the exact minimal primes.

This is a compatibility experiment for the four-vertex path obtained by
moving the cubic through three quadratic factors.  It is deliberately run on
one complete-decomposition chart rather than in the full degree-24
coefficient space.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


W, Z = sp.symbols("W Z")
a1, a2, b1, b2, b3 = sp.symbols("a1 a2 b1 b2 b3")
PARAMETERS = (a1, a2, b1, b2, b3)


def dickson(degree: int, variable: sp.Expr, parameter: sp.Expr) -> sp.Expr:
    """Return the monic Dickson polynomial D_degree(variable, parameter)."""

    previous, current = sp.Integer(2), variable
    for _ in range(2, degree + 1):
        previous, current = current, sp.expand(variable * current - parameter * previous)
    return current


def canonical_residuals(polynomial: sp.Expr, a: int, b: int) -> list[sp.Expr]:
    """Pull the canonical C_(a,b) residuals back to ``polynomial``."""

    degree = a * b
    source = sp.Poly(sp.expand(polynomial), W)
    inner = {j: sp.symbols(f"b{a}{b}_{j}") for j in range(1, b)}
    outer = {j: sp.symbols(f"a{a}{b}_{j}") for j in range(1, a)}
    B = W**b + sum(inner[j] * W**j for j in inner)
    A = Z**a + sum(outer[j] * Z**j for j in outer)
    composition = sp.Poly(sp.expand(A.subs(Z, B)), W)

    reconstruction: dict[sp.Symbol, sp.Expr] = {}
    used_degrees = set()
    for j in range(b - 1, 0, -1):
        coefficient_degree = degree - b + j
        equation = sp.expand(
            composition.nth(coefficient_degree).subs(reconstruction)
            - source.nth(coefficient_degree)
        )
        reconstruction[inner[j]] = sp.factor(sp.solve(equation, inner[j])[0])
        used_degrees.add(coefficient_degree)
    for j in range(a - 1, 0, -1):
        coefficient_degree = j * b
        equation = sp.expand(
            composition.nth(coefficient_degree).subs(reconstruction)
            - source.nth(coefficient_degree)
        )
        reconstruction[outer[j]] = sp.factor(sp.solve(equation, outer[j])[0])
        used_degrees.add(coefficient_degree)

    residuals = []
    for coefficient_degree in range(2, degree):
        if coefficient_degree in used_degrees:
            continue
        residual = sp.factor(
            sp.together(
                composition.nth(coefficient_degree).subs(reconstruction)
                - source.nth(coefficient_degree)
            ).as_numer_denom()[0]
        )
        if residual == 0:
            continue
        primitive = sp.primitive(sp.Poly(residual, *PARAMETERS))[1].as_expr()
        residuals.append(sp.factor(primitive))
    return residuals


def singular_minimal_primes(equations: list[sp.Expr], timeout: int = 1800) -> str:
    singular = shutil.which("Singular")
    assert singular is not None
    serialized = [str(equation).replace("**", "^") for equation in equations]
    program = (
        'LIB "primdec.lib";\n'
        f'ring q=0,({",".join(map(str, PARAMETERS))}),dp;\n'
        "option(redSB);\n"
        f'ideal I={",".join(serialized)};\n'
        'print("IDEAL_DIM"); dim(std(I));\n'
        "list L=minAssGTZ(I);\n"
        'print("COMPONENTS"); size(L);\n'
        "for (int i=1; i<=size(L); i++) {\n"
        '  print("COMPONENT"); print(i); print(dim(std(L[i]))); print(L[i]);\n'
        "}\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "COMPONENTS" in result.stdout, result.stdout + result.stderr
    return result.stdout


def main() -> None:
    cubic = Z**3 + a2 * Z**2 + a1 * Z
    quadratic_1 = Z**2 + b1 * Z
    quadratic_2 = Z**2 + b2 * Z
    quadratic_3 = W**2 + b3 * W
    polynomial = sp.expand(
        cubic.subs(
            Z,
            quadratic_1.subs(
                Z,
                quadratic_2.subs(Z, quadratic_3),
            ),
        )
    )
    assert sp.degree(polynomial, W) == 24

    equations_by_cut = {
        (2, 12): canonical_residuals(polynomial, 2, 12),
        (4, 6): canonical_residuals(polynomial, 4, 6),
        (8, 3): canonical_residuals(polynomial, 8, 3),
    }
    all_equations = []
    for cut, equations in equations_by_cut.items():
        print(f"C_{cut}: {len(equations)} nonzero pulled residuals")
        all_equations.extend(equations)

    t, parameter = sp.symbols("t parameter")
    c2 = dickson(2, t, parameter)
    c4 = dickson(4, t, parameter)
    c8 = dickson(8, t, parameter)
    dickson_substitution = {
        b3: 2 * t,
        b2: 2 * c2,
        b1: 2 * c4,
        a2: 3 * c8,
        a1: 3 * (c8**2 - parameter**8),
    }
    assert all(
        sp.factor(equation.subs(dickson_substitution)) == 0
        for equation in all_equations
    )
    expected = sp.expand(
        dickson(24, W + t, parameter) - dickson(24, t, parameter)
    )
    assert sp.expand(polynomial.subs(dickson_substitution) - expected) == 0

    output = singular_minimal_primes(all_equations)
    compact = " ".join(output.split())
    assert "IDEAL_DIM 2" in compact, output
    assert "COMPONENTS 1" in compact, output
    assert "COMPONENT 1 2" in compact, output
    print(output)
    print("PASS: the degree-24 all-order path has one Dickson/Chebyshev component")


if __name__ == "__main__":
    main()

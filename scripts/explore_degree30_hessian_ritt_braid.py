#!/usr/bin/env python3
"""Test the degree-30 all-factor-order braid on a prime-word chart.

The generic chart has complete degree word ``2 o 3 o 5`` and hence lies on
C_(2,15) and C_(6,5).  The other four proper cuts are imposed by pulling back
their canonical composition residuals.  Exact minimal primes distinguish a
genuine braid-compatible family from a coincidence of pairwise Ritt moves.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


W, Z = sp.symbols("W Z")
p1, q1, q2, r1, r2, r3, r4 = sp.symbols("p1 q1 q2 r1 r2 r3 r4")
PARAMETERS = (p1, q1, q2, r1, r2, r3, r4)


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


def singular_minimal_primes(equations: list[sp.Expr], timeout: int = 3600) -> str:
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
    factor_2 = Z**2 + p1 * Z
    factor_3 = Z**3 + q2 * Z**2 + q1 * Z
    factor_5 = W**5 + r4 * W**4 + r3 * W**3 + r2 * W**2 + r1 * W
    polynomial = sp.expand(factor_2.subs(Z, factor_3.subs(Z, factor_5)))
    assert sp.degree(polynomial, W) == 30

    equations_by_cut = {
        (3, 10): canonical_residuals(polynomial, 3, 10),
        (5, 6): canonical_residuals(polynomial, 5, 6),
        (10, 3): canonical_residuals(polynomial, 10, 3),
        (15, 2): canonical_residuals(polynomial, 15, 2),
    }
    all_equations = []
    for cut, equations in equations_by_cut.items():
        print(f"C_{cut}: {len(equations)} nonzero pulled residuals")
        all_equations.extend(equations)

    t, parameter = sp.symbols("t parameter")
    c5 = dickson(5, t, parameter)
    c15 = dickson(15, t, parameter)
    dickson_substitution = {
        r4: 5 * t,
        r3: 10 * t**2 - 5 * parameter,
        r2: 10 * t**3 - 15 * parameter * t,
        r1: 5 * t**4 - 15 * parameter * t**2 + 5 * parameter**2,
        q2: 3 * c5,
        q1: 3 * (c5**2 - parameter**5),
        p1: 2 * c15,
    }
    assert all(
        sp.factor(equation.subs(dickson_substitution)) == 0
        for equation in all_equations
    )
    expected = sp.expand(
        dickson(30, W + t, parameter) - dickson(30, t, parameter)
    )
    assert sp.expand(polynomial.subs(dickson_substitution) - expected) == 0

    output = singular_minimal_primes(all_equations)
    compact = " ".join(output.split())
    assert "IDEAL_DIM 2" in compact, output
    assert "COMPONENTS 1" in compact, output
    assert "COMPONENT 1 2" in compact, output
    print(output)
    print("PASS: the degree-30 all-order braid has one Dickson/Chebyshev component")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Restrict the degree-18 Hessian--Ritt loci to the 2-versus-9 power core.

Ritt's second theorem puts the monic source-translated power collision in the
form

    f = (y R(y^2))^2,  y = W-s,  deg(R)=4.

It has both 2-o-9 and 9-o-2 decompositions.  This script reconstructs the
canonical target factors for 3-o-6 and 6-o-3, pulls back their residual
equations to QQ[s,r0,r1,r2,r3], and asks Singular for exact minimal primes.

The calculation certifies the unique minimal prime, hence the reduced support.
It does not test whether the unreduced pullback ideal has nilpotents or
embedded components.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


W, Z = sp.symbols("W Z")
s, r0, r1, r2, r3 = sp.symbols("s r0 r1 r2 r3")
parameters = (s, r0, r1, r2, r3)


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
        residuals.append(sp.primitive(sp.Poly(residual, *parameters))[1].as_expr())
    assert len(residuals) == degree - a - b
    return residuals


def singular_minimal_primes(equations: list[sp.Expr], timeout: int = 600) -> str:
    singular = shutil.which("Singular")
    assert singular is not None
    serialized = [str(equation).replace("**", "^") for equation in equations]
    program = (
        'LIB "primdec.lib";\n'
        f'ring q=0,({",".join(map(str, parameters))}),dp;\n'
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


def assert_one_threefold(output: str) -> None:
    compact = " ".join(output.split())
    assert "IDEAL_DIM 3" in compact, output
    assert "COMPONENTS 1" in compact, output
    assert "COMPONENT 1 3" in compact, output


def normalized_clean_seed(polynomial: sp.Expr) -> sp.Expr | None:
    tangent = sp.diff(polynomial, W).subs(W, 0)
    endpoint = sp.factor(
        polynomial.subs(W, 1) - polynomial.subs(W, 0) - tangent
    )
    if endpoint != 0:
        return None
    derivative_gap = sp.factor(sp.diff(polynomial, W).subs(W, 1) - tangent)
    if derivative_gap == 0:
        return None
    H = sp.factor(
        -(polynomial - polynomial.subs(W, 0) - tangent * W) / derivative_gap
    )
    K = sp.diff(H, W, 2)
    Q = sp.cancel(H / W**2)
    clean_quantities = (
        sp.discriminant(K, W),
        K.subs(W, 0),
        sp.discriminant(Q, W),
        K.subs(W, 1) + 2,
    )
    return None if any(quantity == 0 for quantity in clean_quantities) else H


def main() -> None:
    y = W - s
    R = y**8 + r3 * y**6 + r2 * y**4 + r1 * y**2 + r0
    h = sp.expand(y * R)
    power_collision = sp.expand(h**2)
    assert sp.degree(power_collision, W) == 18

    residuals_36 = canonical_residuals(power_collision, 3, 6)
    residuals_63 = canonical_residuals(power_collision, 6, 3)
    nested_substitution = {
        r2: r3**2 / 3,
        r0: r1 * r3 / 3 - r3**4 / 81,
    }
    assert all(
        sp.factor(equation.subs(nested_substitution)) == 0
        for equation in residuals_36 + residuals_63
    )

    u = r3 / 3
    v = r1 - r3**3 / 27
    B = y**3 + u * y
    A_of_B = B**3 + v * B
    assert sp.factor((h - A_of_B).subs(nested_substitution)) == 0

    clean_y = W + 3
    clean_B = clean_y**3 + 5 * clean_y
    clean_A = clean_B**3 - sp.Rational(721476, 31) * clean_B
    assert normalized_clean_seed(sp.expand(clean_A**2)) is not None

    print(f"C_36 residuals: {len(residuals_36)}")
    output_36 = singular_minimal_primes(residuals_36)
    assert_one_threefold(output_36)
    print(output_36)
    print(f"C_63 residuals: {len(residuals_63)}")
    output_63 = singular_minimal_primes(residuals_63)
    assert_one_threefold(output_63)
    print(output_63)
    print("C_36 and C_63 residuals")
    output_all = singular_minimal_primes(residuals_36 + residuals_63)
    assert_one_threefold(output_all)
    print(output_all)
    print("PASS: both middle factor orders have the same unique nested-power minimal prime")
    print("PASS: the nested-power threefold meets the marked clean open")


if __name__ == "__main__":
    main()

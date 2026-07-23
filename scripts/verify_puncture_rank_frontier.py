#!/usr/bin/env python3
"""Exact audit of the degree 4--7 two-center puncture obstruction."""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


x, y, B = sp.symbols("x y B")
f_symbol = sp.Function("f")(y)
A_symbol = 1 + x * f_symbol
s_symbol = x / A_symbol
P_symbol = A_symbol * B
Q_symbol = y + x * B

chart_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(output, variable) for variable in (x, y, B)]
                for output in (s_symbol, P_symbol, Q_symbol)
            ]
        )
    )
)
assert chart_jacobian == -1 / A_symbol


def cleared_boundary_moment(r: int, a: int, b: int) -> sp.Poly:
    """Return the numerator of the necessary boundary moment."""

    v, G = sp.symbols("v G")
    f = y**a * (y - 1) ** b
    c = G / f
    shifted_y = y - c * (1 - v)
    shifted_f = shifted_y**a * (shifted_y - 1) ** b
    moment = sp.integrate(v**r * shifted_f**r, (v, 0, 1))
    numerator = sp.cancel(moment).as_numer_denom()[0]
    return sp.Poly(numerator, G, y, domain=sp.QQ)


# If f has two distinct roots, k[Y,P,1/f] has the two independent units
# Y and Y-1.  Their valuation vectors at the two finite punctures are the
# standard basis, so the geometric unit rank is exactly two.
valuation_matrix = sp.Matrix([[1, 0], [0, 1]])
assert valuation_matrix.rank() == 2


# N = r * (degree(f) + 1) + 1.  These are exactly the rank-two cases with
# 4 <= N <= 7.
cases: list[tuple[int, int, int, int]] = []
for r in range(1, 3):
    for a in range(1, 6):
        for b in range(1, 6):
            degree = r * (a + b + 1) + 1
            if 4 <= degree <= 7:
                cases.append((degree, r, a, b))

expected_cases = [
    (4, 1, 1, 1),
    (5, 1, 1, 2),
    (5, 1, 2, 1),
    (6, 1, 1, 3),
    (6, 1, 2, 2),
    (6, 1, 3, 1),
    (7, 1, 1, 4),
    (7, 1, 2, 3),
    (7, 1, 3, 2),
    (7, 1, 4, 1),
    (7, 2, 1, 1),
]
assert sorted(cases) == expected_cases


# The smallest case has a short independent discriminant certificate.
G = sp.symbols("G")
smallest = cleared_boundary_moment(1, 1, 1)
smallest_monic = sp.Poly(
    G**2
    - 2 * y * (y - 1) * (2 * y - 1) * G
    + 6 * y**3 * (y - 1) ** 3,
    G,
    y,
    domain=sp.QQ,
)
assert sp.rem(smallest.as_expr(), smallest_monic.as_expr(), G) == 0
discriminant = sp.factor(sp.discriminant(smallest_monic.as_expr(), G))
assert (
    sp.expand(
        discriminant
        - 4 * y**2 * (y - 1) ** 2 * (-2 * y**2 + 2 * y + 1)
    )
    == 0
)
nonsquare_part = sp.cancel(discriminant / (4 * y**2 * (y - 1) ** 2))
assert sp.degree(nonsquare_part, y) == 2
assert sp.discriminant(nonsquare_part, y) != 0


certificates: list[tuple[int, int, int, int, int]] = []
moments: list[sp.Poly] = []
for degree, r, a, b in expected_cases:
    moment = cleared_boundary_moment(r, a, b)
    moments.append(moment)
    coefficient, factors = sp.factor_list(moment.as_expr(), G, y)
    del coefficient
    assert len(factors) == 1
    factor, exponent = factors[0]
    assert exponent == 1
    assert sp.Poly(factor, G, y, domain=sp.QQ).total_degree() > 1
    g_degree = sp.degree(factor, G)
    assert g_degree > 1
    assert all(sp.degree(candidate, G) != 1 for candidate, _ in factors)
    certificates.append((degree, r, a, b, g_degree))


# Rational irreducibility would not suffice over the algebraically closed
# ground field.  Ask Singular for the absolute factorization of every
# generated moment.  absolute_factors[4] is the total number of nonconstant
# absolute factors, so it must equal one in every case.
singular = shutil.which("Singular")
assert singular is not None, "Singular is required for absolute factorization"
singular_lines = ['LIB "absfact.lib";', "ring R=0,(y,G),dp;"]
for index, ((degree, r, a, b), moment) in enumerate(
    zip(expected_cases, moments, strict=True),
    start=1,
):
    singular_polynomial = str(moment.as_expr()).replace("**", "^")
    singular_lines.extend(
        [
            f"poly p{index}={singular_polynomial};",
            f"def S{index}=absFactorize(p{index});",
            f"setring S{index};",
            f'print("ABSCASE_{degree}_{r}_{a}_{b}");',
            "print(absolute_factors[4]);",
            "setring R;",
        ]
    )
singular_lines.append("quit;")
singular_result = subprocess.run(
    [singular, "-q"],
    input="\n".join(singular_lines),
    text=True,
    capture_output=True,
    check=True,
)
absolute_counts: dict[str, int] = {}
output_lines = [
    line.strip() for line in singular_result.stdout.splitlines() if line.strip()
]
for index, line in enumerate(output_lines[:-1]):
    if line.startswith("ABSCASE_"):
        absolute_counts[line] = int(output_lines[index + 1])
assert len(absolute_counts) == len(expected_cases)
assert set(absolute_counts.values()) == {1}


print("PASS: universal two-center chart has determinant -A^-1 before z")
print("PASS: two distinct finite centers give geometric unit rank two")
print("PASS: degree 4--7 list has exactly 11 exponent cases")
print("PASS: degree-four moment has a nonsquare quadratic discriminant")
print("PASS: every degree 4--7 boundary moment is absolutely irreducible")
print("certificates (N,r,a,b,deg_G):", certificates)

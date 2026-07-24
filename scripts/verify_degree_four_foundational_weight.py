#!/usr/bin/env python3
"""Exact exclusion of degree-four Keller counterexamples at weights (1,-1,-2).

The symbolic part enumerates the complete equivariant support, derives the
Keller coefficient ideal, and verifies the three explicit automorphism
families.  Singular independently computes the radical of the coefficient
ideal and checks it against the intersection of the three displayed primes.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


x, y, z = sp.symbols("x y z")
variables = (x, y, z)
weights = (1, -1, -2)


def equivariant_exponents(target_weight: int, degree_bound: int = 4):
    """Return every positive-degree exponent of bounded total degree."""
    result = []
    for total_degree in range(1, degree_bound + 1):
        for x_degree in range(total_degree + 1):
            for y_degree in range(total_degree - x_degree + 1):
                z_degree = total_degree - x_degree - y_degree
                exponent = (x_degree, y_degree, z_degree)
                if sum(entry * weight for entry, weight in zip(exponent, weights)) == target_weight:
                    result.append(exponent)
    return result


expected_supports = {
    1: [(1, 0, 0), (2, 1, 0), (3, 0, 1)],
    -1: [(0, 1, 0), (1, 0, 1), (1, 2, 0), (2, 1, 1)],
    -2: [(0, 0, 1), (0, 2, 0), (1, 1, 1), (1, 3, 0), (2, 0, 2)],
}
assert {
    target_weight: equivariant_exponents(target_weight)
    for target_weight in weights
} == expected_supports

a, b, c, d, e, f, g, h, i = sp.symbols("a b c d e f g h i")
parameters = (a, b, c, d, e, f, g, h, i)

F = sp.Matrix(
    [
        x + a * x**2 * y + b * x**3 * z,
        y + c * x * z + d * x * y**2 + e * x**2 * y * z,
        z + f * y**2 + g * x * y * z + h * x * y**3 + i * x**2 * z**2,
    ]
)

determinant_remainder = sp.Poly(
    sp.expand(F.jacobian(variables).det() - 1), *variables
)
equations_by_monomial = dict(determinant_remainder.terms())

expected_equations = {
    (6, 0, 3): 4 * b * e * i,
    (5, 3, 1): -4 * b * e * h,
    (5, 1, 2): 2 * a * e * i + 8 * b * d * i + b * e * g,
    (4, 4, 0): -h * (5 * a * e - b * d),
    (4, 2, 1): 6 * a * d * i - a * e * g - 6 * b * c * h + 5 * b * d * g - 2 * b * e * f,
    (4, 0, 2): -2 * b * c * g + 3 * b * e + 4 * b * i + 2 * e * i,
    (3, 3, 0): -5 * a * c * h + 3 * a * d * g - 4 * a * e * f + 2 * b * d * f - b * h - 3 * e * h,
    (3, 1, 1): -2 * (a * c * g - 2 * a * i + 2 * b * c * f - 3 * b * d - b * g - 2 * d * i),
    (2, 2, 0): -4 * a * c * f + 3 * a * d + 2 * a * g - 3 * c * h + 2 * d * g - 2 * e * f,
    (2, 0, 1): -a * c + 3 * b - c * g + e + 2 * i,
    (1, 1, 0): 2 * a - 2 * c * f + 2 * d + g,
}
assert equations_by_monomial.keys() == expected_equations.keys()
assert all(
    sp.expand(equations_by_monomial[monomial] - expected_equations[monomial]) == 0
    for monomial in expected_equations
)


def singular_polynomial(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
assert singular is not None, "Singular is required for the radical certificate"
ideal_generators = ",".join(
    singular_polynomial(expression) for expression in expected_equations.values()
)
singular_program = f"""
ring r=0,(a,b,c,d,e,f,g,h,i),dp;
ideal I={ideal_generators};
ideal P1=i,h,g,e,b,a,c*f-d;
ideal P2=i,g,e,d,c,b,a;
ideal P3=h,e,d,b,a,g^2-4*f*i,c*g-2*i,2*c*f-g;
ideal J=intersect(P1,intersect(P2,P3));
LIB "primdec.lib";
ideal R=radical(I);
ideal left=reduce(std(J),std(R));
ideal right=reduce(std(R),std(J));
if ((size(left)==0) && (size(right)==0))
{{
  print("PASS_RADICAL");
}}
else
{{
  print("FAIL_RADICAL");
}}
quit;
"""
singular_result = subprocess.run(
    [singular, "-q"],
    input=singular_program,
    text=True,
    capture_output=True,
    check=True,
)
assert "PASS_RADICAL" in singular_result.stdout, singular_result.stdout
assert "FAIL_RADICAL" not in singular_result.stdout, singular_result.stdout

# The radical has the compact generating set recorded in the theorem note.
radical_generators = (
    e,
    b,
    a,
    h * i,
    d * i,
    g * h,
    d * h,
    c * h,
    g**2 - 4 * f * i,
    d * g,
    c * g - 2 * i,
    2 * c * f - 2 * d - g,
)
for prime_substitution in (
    {a: 0, b: 0, e: 0, g: 0, h: 0, i: 0, d: c * f},
    {a: 0, b: 0, c: 0, d: 0, e: 0, g: 0, i: 0},
    {a: 0, b: 0, d: 0, e: 0, h: 0, g: 2 * c * f, i: c**2 * f},
):
    assert all(sp.expand(entry.subs(prime_substitution)) == 0 for entry in radical_generators)
    assert sp.expand(determinant_remainder.as_expr().subs(prime_substitution)) == 0

# Family P1:
#   (x, y+c*x*(z+f*y^2), z+f*y^2).
X, Y, Z = sp.symbols("X Y Z")
P1 = sp.Matrix([x, y + c * x * (z + f * y**2), z + f * y**2])
P1_inverse = sp.Matrix(
    [
        X,
        Y - c * X * Z,
        Z - f * (Y - c * X * Z) ** 2,
    ]
)

# Family P2:
#   (x, y, z+f*y^2+h*x*y^3).
P2 = sp.Matrix([x, y, z + f * y**2 + h * x * y**3])
P2_inverse = sp.Matrix([X, Y, Z - f * Y**2 - h * X * Y**3])

# Family P3:
#   let Y=y+c*x*z, then (x,Y,z+f*Y^2).
P3 = sp.Matrix([x, y + c * x * z, z + f * (y + c * x * z) ** 2])
P3_inverse = sp.Matrix(
    [
        X,
        Y - c * X * (Z - f * Y**2),
        Z - f * Y**2,
    ]
)


def verify_inverse(mapping: sp.Matrix, inverse: sp.Matrix) -> None:
    forward_then_inverse = inverse.subs(
        {X: mapping[0], Y: mapping[1], Z: mapping[2]},
        simultaneous=True,
    )
    inverse_then_forward = mapping.subs(
        {x: inverse[0], y: inverse[1], z: inverse[2]},
        simultaneous=True,
    )
    assert all(
        sp.expand(got - expected) == 0
        for got, expected in zip(forward_then_inverse, (x, y, z))
    )
    assert all(
        sp.expand(got - expected) == 0
        for got, expected in zip(inverse_then_forward, (X, Y, Z))
    )


for mapping, inverse in (
    (P1, P1_inverse),
    (P2, P2_inverse),
    (P3, P3_inverse),
):
    assert sp.expand(mapping.jacobian(variables).det() - 1) == 0
    verify_inverse(mapping, inverse)

print("PASS: the complete degree-four equivariant support has nine nonlinear coefficients")
print("PASS: the Keller identity gives exactly eleven coefficient equations")
print("PASS: Singular certifies radical(I) = P1 intersect P2 intersect P3")
print("PASS: every component is an explicit tame-automorphism family")

#!/usr/bin/env python3
"""Verify generic Schur rigidity on a two-parameter sextic surface.

The sextic is

  (x^6+y^6+z^6)/30
  + mu*x^2*y^2*z^2
  + nu*sum_{i != j} x_i^4*x_j^2.

On nu != 0, six boundary coefficients solve the quadratic Schur quotient.
After clearing the resulting nu^2 denominator, 114 intrinsic equations
remain in the fifteen quartic coefficients.  Over Q(mu,nu), their exact
Groebner basis makes every quartic coefficient nilpotent of exponent three.

This proves rigidity at the generic point of the parameter surface.  It
does not claim that the specialization locus inside nu != 0 is empty.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import sympy as sp


x, y, z, mu, nu = sp.symbols("x y z mu nu")
quartic_coefficients = sp.symbols("s0:15")

quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quartic_coefficients, quartic_monomials
    )
)

mixed_42 = sum(
    left**4 * right**2
    for left in (x, y, z)
    for right in (x, y, z)
    if left != right
)
h6 = (
    (x**6 + y**6 + z**6) / 30
    + mu * x**2 * y**2 * z**2
    + nu * mixed_42
)
hessian = sp.hessian(h6, (x, y, z))
hessian_determinant = sp.expand(hessian.det())
gradient = sp.Matrix(
    [sp.diff(quartic, variable) for variable in (x, y, z)]
)
schur_numerator = sp.expand(
    (gradient.T * hessian.adjugate() * gradient)[0]
)

quadratic_monomials = (
    x**2,
    y**2,
    z**2,
    x * y,
    x * z,
    y * z,
)
degree_14_monomials = [
    x**i * y**j * z ** (14 - i - j)
    for i in range(15)
    for j in range(15 - i)
]
quotient_matrix = sp.Matrix(
    [
        [
            sp.Poly(
                hessian_determinant * quadratic_monomial,
                x,
                y,
                z,
            ).coeff_monomial(degree_14_monomial)
            for quadratic_monomial in quadratic_monomials
        ]
        for degree_14_monomial in degree_14_monomials
    ]
)
numerator_vector = sp.Matrix(
    [
        sp.Poly(schur_numerator, x, y, z).coeff_monomial(
            monomial
        )
        for monomial in degree_14_monomials
    ]
)

# These are z^14, y*z^13, y^2*z^12, x*z^13, x*y*z^12,
# and x^2*z^12.  Their quotient matrix has determinant 4096*nu^12.
pivot_rows = (0, 1, 2, 15, 16, 29)
pivot_matrix = quotient_matrix[list(pivot_rows), :]
assert sp.factor(pivot_matrix.det()) == 4096 * nu**12

s = quartic_coefficients

# To avoid repeated rational-matrix simplification, store 2*nu^2 times
# the six exact quotient solutions.
quotient_numerators = (
    -mu * s[1] ** 2
    + 640 * nu**3 * s[0] ** 2
    - 96 * nu**2 * s[0] * s[9]
    + 2 * nu**2 * s[5] ** 2
    + 2 * nu * s[1] * s[10]
    + 6 * nu * s[12] * s[5]
    - 6 * nu * s[5] ** 2
    + nu * s[6] ** 2
    + 4 * nu * s[9] ** 2,
    -mu * s[5] ** 2
    + 640 * nu**3 * s[0] ** 2
    - 96 * nu**2 * s[0] * s[2]
    + 2 * nu**2 * s[1] ** 2
    - 6 * nu * s[1] ** 2
    + 6 * nu * s[1] * s[3]
    + 4 * nu * s[2] ** 2
    + 2 * nu * s[5] * s[7]
    + nu * s[6] ** 2,
    nu * (32 * nu * s[0] ** 2 + s[1] ** 2 + s[5] ** 2),
    -4
    * (
        mu * s[1] * s[5]
        + 24 * nu**2 * s[0] * s[6]
        - nu**2 * s[1] * s[5]
        - nu * s[1] * s[7]
        - nu * s[10] * s[5]
        - nu * s[2] * s[6]
        - nu * s[6] * s[9]
    ),
    -2
    * nu
    * (8 * nu * s[0] * s[5] - s[1] * s[6] - 2 * s[5] * s[9]),
    -2
    * nu
    * (8 * nu * s[0] * s[1] - 2 * s[1] * s[2] - s[5] * s[6]),
)

for pivot_position, row in enumerate(pivot_rows):
    pivot_remainder = sp.expand(
        2 * nu**2 * numerator_vector[row]
        - sum(
            quotient_matrix[row, column]
            * quotient_numerators[column]
            for column in range(6)
        )
    )
    assert pivot_remainder == 0, pivot_position

intrinsic_equations = []
for row in range(len(degree_14_monomials)):
    equation = sp.expand(
        2 * nu**2 * numerator_vector[row]
        - sum(
            quotient_matrix[row, column]
            * quotient_numerators[column]
            for column in range(6)
        )
    )
    if equation != 0:
        intrinsic_equations.append(equation)
assert len(intrinsic_equations) == 114


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact generic check")

program = f"""
ring rr=(0,mu,nu),({",".join(map(str, quartic_coefficients))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, intrinsic_equations))};
ideal G=slimgb(I);
ideal M={",".join(map(str, quartic_coefficients))};
ideal GM=std(M);
print(
  "GENERIC "
  +string(size(I))+" "
  +string(size(G))+" "
  +string(size(reduce(I,GM)))
);
int i;
poly cube;
for (i=1;i<=size(M);i++)
{{
  cube=M[i]^3;
  print("CUBE "+string(i)+" "+string(reduce(cube,G)==0));
}}
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
    timeout=180,
)
if completed.stderr.strip():
    raise RuntimeError(completed.stderr)

generic_marker = re.search(
    r"(?m)^GENERIC (\d+) (\d+) (\d+)$", completed.stdout
)
assert generic_marker is not None
assert tuple(map(int, generic_marker.groups())) == (114, 117, 0)
cube_markers = re.findall(
    r"(?m)^CUBE (\d+) ([01])$", completed.stdout
)
assert len(cube_markers) == 15
assert all(success == "1" for _, success in cube_markers)

print("PASS: six quotient pivots have determinant 4096*nu^12")
print("PASS: quotient elimination leaves 114 intrinsic equations")
print("PASS: the generic Groebner basis has 117 elements")
print("PASS: all fifteen quartic coefficient cubes reduce to zero")
print("SCOPE: generic rigidity only; exceptional nu!=0 fibers remain open")

#!/usr/bin/env python3
"""Classify even permutation-invariant Schur pairs on the HC4 surface.

Use

    R  = x^2+y^2+z^2,
    P2 = x^2*y^2+x^2*z^2+y^2*z^2,
    P3 = x^2*y^2*z^2,
    h6 = R^3/30 + A*R*P2 + B*P3.

The projective even symmetric quartics are a*R^2+b*P2.  Exact radical
calculations on a=1 and a=0 show that only the radial and Fermat points
carry such a nonzero quartic Schur norm.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, b, q = sp.symbols("A B b q")
variables = (x, y, z)

radius = x**2 + y**2 + z**2
pair_sum = x**2 * y**2 + x**2 * z**2 + y**2 * z**2
triple_product = x**2 * y**2 * z**2
sextic = (
    radius**3 / 30
    + A * radius * pair_sum
    + B * triple_product
)
hessian = sp.hessian(sextic, variables)


def coefficient_equations(quartic: sp.Expr) -> list[sp.Expr]:
    gradient = sp.Matrix(
        [sp.diff(quartic, variable) for variable in variables]
    )
    remainder = sp.expand(
        (gradient.T * hessian.adjugate() * gradient)[0]
        - hessian.det() * q * radius
    )
    return [
        sp.expand(sp.together(coefficient).as_numer_denom()[0])
        for coefficient in dict.fromkeys(
            sp.Poly(remainder, *variables).coeffs()
        )
    ]


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


affine_equations = coefficient_equations(radius**2 + b * pair_sum)
infinity_equations = coefficient_equations(pair_sum)
assert len(affine_equations) == 8
assert len(infinity_equations) == 8

singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact radical check")

program = f"""
LIB "primdec.lib";
ring affine_ring=0,(q,b,A,B),dp;
ideal I={",".join(map(singular_expression, affine_equations))};
ideal radial=q-16,b,A,B;
ideal fermat=q-16,b+2,10*A+1,10*B-1;
ideal expected=intersect(radial,fermat);
ideal actual=std(radical(I));
ideal left=reduce(actual,std(expected));
ideal right=reduce(std(expected),actual);
print(
  "AFFINE "
  +string(size(left))+" "
  +string(size(right))
);

ring infinity_ring=0,(q,A,B),dp;
ideal I={",".join(map(singular_expression, infinity_equations))};
ideal actual=std(radical(I));
print("INFINITY "+string(reduce(1,actual)==0));
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
    timeout=120,
)
if completed.stderr.strip():
    raise RuntimeError(completed.stderr)

affine_marker = re.search(
    r"(?m)^AFFINE (\d+) (\d+)$", completed.stdout
)
assert affine_marker is not None
assert tuple(map(int, affine_marker.groups())) == (0, 0)
infinity_marker = re.search(
    r"(?m)^INFINITY ([01])$", completed.stdout
)
assert infinity_marker is not None
assert infinity_marker.group(1) == "1"

# Translate the two deformation points back to (mu,nu):
#
#   A=nu-1/10,  B=mu-3*nu+1/10.
#
# The radial point A=B=0 is (mu,nu)=(1/5,1/10).
# The Fermat point A=-1/10, B=1/10 is (mu,nu)=(0,0).
print("PASS: the affine chart has exactly radial and Fermat support")
print("PASS: the pure-P2 projective chart is empty")
print("PASS: the even symmetric Schur locus has no parameter curve")
print("SCOPE: quartics in span(R^2,P2); nonsymmetric quartics remain open")

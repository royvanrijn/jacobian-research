#!/usr/bin/env python3
"""Verify the full symmetric-sextic quintic Schur obstruction.

For

    h6 = (x^6+y^6+z^6)/30 + mu*x^2*y^2*z^2

and a generic ternary quartic s, compare

    grad(s)^T adj(Hess(h6)) grad(s)

with det(Hess(h6)) times a generic quadratic quotient.  Over Q[mu], the
coefficient ideal saturated by mu has radical equal to the origin in all
quartic and quotient coefficients.  Thus every nonzero-mu member of the
pencil has no nonzero polynomial Schur norm.
"""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp


x, y, z, mu = sp.symbols("x y z mu")
quartic_coefficients = sp.symbols("s0:15")
quotient_coefficients = sp.symbols("q0:6")

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

h6 = (x**6 + y**6 + z**6) / 30 + mu * x**2 * y**2 * z**2
hessian = sp.hessian(h6, (x, y, z))
gradient = sp.Matrix(
    [sp.diff(quartic, variable) for variable in (x, y, z)]
)
schur_numerator = sp.expand(
    (gradient.T * hessian.adjugate() * gradient)[0]
)

q0, q1, q2, q3, q4, q5 = quotient_coefficients
quadratic_quotient = (
    q0 * x**2
    + q1 * y**2
    + q2 * z**2
    + q3 * x * y
    + q4 * x * z
    + q5 * y * z
)
remainder = sp.expand(
    schur_numerator - hessian.det() * quadratic_quotient
)
coefficient_equations = [
    coefficient
    for _, coefficient in sp.Poly(remainder, x, y, z).terms()
]
assert len(coefficient_equations) == 111


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact saturation")

all_unknowns = quartic_coefficients + quotient_coefficients
program = f"""
LIB "elim.lib";
ring rr=0,({",".join(map(str, (mu,) + all_unknowns))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, coefficient_equations))};
ideal M={",".join(map(str, all_unknowns))};
ideal MU=mu;
list saturation=sat(I,MU);
ideal S=saturation[1];
ideal G=slimgb(S);
ideal GM=std(M);
print(
  "SATURATION "
  +string(size(I))+" "
  +string(size(G))+" "
  +string(size(reduce(S,GM)))
);
int i;
poly fourthPower;
for (i=1;i<=size(M);i++)
{{
  fourthPower=M[i]^4;
  print(
    "FOURTH_POWER "+string(i)+" "
    +string(reduce(fourthPower,G)==0)
  );
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

saturation_marker = re.search(
    r"(?m)^SATURATION (\d+) (\d+) (\d+)$",
    completed.stdout,
)
assert saturation_marker is not None
assert tuple(map(int, saturation_marker.groups())) == (111, 261, 0)

power_markers = re.findall(
    r"(?m)^FOURTH_POWER (\d+) ([01])$",
    completed.stdout,
)
assert len(power_markers) == len(all_unknowns) == 21
assert all(success == "1" for _, success in power_markers)


# Keep the exact Singular transcript available only on explicit request;
# the deterministic marker data above are sufficient for the normal replay.
assert Path(singular).is_file()

print("PASS: the symmetric sextic Schur system has 111 exact equations")
print("PASS: saturation by mu has a 261-element rational Groebner basis")
print("PASS: its radical is the 21-variable coefficient origin")
print("PASS: fourth powers certify all quartic and quotient coefficients")
print("SCOPE: every nonzero-mu symmetric sextic has only s4=0")

#!/usr/bin/env python3
"""Verify the first two exact faces over the pure-seventh scalar top."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "artifacts" / "generated-results" / "hc4_pure_septic_opening.json"
)

x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def bordered_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in variables])
    hessian = sp.hessian(polynomial, variables)
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


def homogeneous_face(polynomial: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * x**monomial[0] * y**monomial[1] * z**monomial[2]
            for monomial, coefficient in sp.Poly(polynomial, *variables).terms()
            if sum(monomial) == degree
        )
    )


# The degree-twenty face is computed with every coefficient of the sextic
# correction independent.  It is precisely its passive Hessian determinant.
sextic_monomials = [
    x**i * y**j * z ** (6 - i - j)
    for i in range(7)
    for j in range(7 - i)
]
sextic_coefficients = sp.symbols(f"s0:{len(sextic_monomials)}")
generic_sextic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(sextic_coefficients, sextic_monomials)
)
first_face = homogeneous_face(bordered_invariant(x**7 + generic_sextic), 20)
assert sp.expand(
    first_face - 49 * x**12 * sp.hessian(generic_sextic, (y, z)).det()
) == 0


# The passive singular-Hessian classification gives
# c_6=H_6(x,y)+k*x^5*z.  Retain a completely generic quintic correction in
# the next face.  When H_yy is nonzero, this fixes its complete z^2 part.
a = sp.symbols("a0:7")
k = sp.symbols("k")
binary_sextic = sum(a[index] * x ** (6 - index) * y**index for index in range(7))
quintic_monomials = [
    x**i * y**j * z ** (5 - i - j)
    for i in range(6)
    for j in range(6 - i)
]
quintic_coefficients = sp.symbols(f"q0:{len(quintic_monomials)}")
generic_quintic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(quintic_coefficients, quintic_monomials)
)
second_face = homogeneous_face(
    bordered_invariant(x**7 + binary_sextic + k * x**5 * z + generic_quintic),
    19,
)
assert sp.expand(
    second_face
    - sp.Rational(49, 2)
    * x**12
    * sp.diff(binary_sextic, y, 2)
    * (2 * sp.diff(generic_quintic, z, 2) - sp.Rational(8, 7) * k**2 * x**3)
) == 0


payload = {
    "format": "hc4-pure-septic-opening-v1",
    "status": {
        "id": "HC4RSD34",
        "kind": "exact narrowing theorem",
        "scope": "the pure-seventh scalar leading chart left by HC4RSD33",
        "result": (
            "the sextic correction is passive singular-Hessian; in its binary "
            "normal form, when (H6)_yy is nonzero, the next face fixes the "
            "complete z^2 quintic tail"
        ),
    },
    "faces": {
        "degree_20": "49*x^12*det Hess_(y,z)(c6)",
        "degree_19": (
            "49*x^12*(H6)_yy*(2*(c5)_zz-(8/7)*k^2*x^3)/2"
        ),
    },
    "normal_form": {
        "c7": "x^7",
        "c6": "H6(x,y)+k*x^5*z",
        "c5_when_H6_yy_nonzero": "R5(x,y)+z*P4(x,y)+(2/7)*k^2*x^3*z^2",
    },
    "residual": (
        "the descendants of the moving k!=0, (H6)_yy!=0 chart and the "
        "passive-affine (H6)_yy=0 boundary"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: the pure-seventh first face is the passive sextic Hessian")
print("PASS: the next face fixes the complete quintic z^2 tail")
print("THEOREM: the pure-seventh chart has the stated two-face normal form")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")

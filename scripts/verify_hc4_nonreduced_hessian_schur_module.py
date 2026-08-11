#!/usr/bin/env python3
"""Exact calibrations for the nonreduced Hessian--Schur module note.

This checker verifies only polynomial identities and the numerical degree
table.  The duality, DVR, normalization, and line-bundle arguments are
written proofs in HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md.
"""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def hessian(form: sp.Expr) -> sp.Matrix:
    return sp.hessian(form, variables)


def adjugate(matrix: sp.Matrix) -> sp.Matrix:
    return matrix.adjugate()


def assert_zero(expression: sp.Expr) -> None:
    assert sp.expand(expression) == 0


# Fermat quintic: C=diag(x^3,y^3,z^3), with the full three-channel cubic.
a, b, c = sp.symbols("a b c")
h5_fermat = (x**5 + y**5 + z**5) / 20
s3_fermat = (a * x**3 + b * y**3 + c * z**3) / 3
C5 = hessian(h5_fermat)
d3 = sp.Matrix([sp.diff(s3_fermat, variable) for variable in variables])
Delta5 = sp.factor(C5.det())
norm5 = sp.expand((d3.T * adjugate(C5) * d3)[0])
quotient5 = a**2 * x + b**2 * y + c**2 * z
assert C5 == sp.diag(x**3, y**3, z**3)
assert Delta5 == x**3 * y**3 * z**3
assert_zero(norm5 - Delta5 * quotient5)

P5 = x * y * z
e5 = sp.Matrix([a * y * z, b * x * z, c * x * y])
assert all(sp.expand(value) == 0 for value in C5 * e5 - P5 * d3)
assert_zero((d3.T * e5)[0] - P5 * quotient5)

# Along x=0, the induced rank-two determinant is y^3*z^3, of degree six.
fermat_quintic_defect = sp.factor(C5[1:, 1:].det().subs(x, 0))
assert fermat_quintic_defect == y**3 * z**3
assert sp.Poly(fermat_quintic_defect, y, z).total_degree() == 6


# Fermat sextic: the same module calculation has defect degree eight.
h6_fermat = (x**6 + y**6 + z**6) / 30
s4_fermat = (a * x**4 + b * y**4 + c * z**4) / 4
C6 = hessian(h6_fermat)
d4 = sp.Matrix([sp.diff(s4_fermat, variable) for variable in variables])
Delta6 = sp.factor(C6.det())
norm6 = sp.expand((d4.T * adjugate(C6) * d4)[0])
quotient6 = a**2 * x**2 + b**2 * y**2 + c**2 * z**2
assert C6 == sp.diag(x**4, y**4, z**4)
assert Delta6 == x**4 * y**4 * z**4
assert_zero(norm6 - Delta6 * quotient6)

P6 = x * y * z
e6 = sp.Matrix([a * y * z, b * x * z, c * x * y])
assert all(sp.expand(value) == 0 for value in C6 * e6 - P6 * d4)
assert_zero((d4.T * e6)[0] - P6 * quotient6)

fermat_sextic_defect = sp.factor(C6[1:, 1:].det().subs(x, 0))
assert fermat_sextic_defect == y**4 * z**4
assert sp.Poly(fermat_sextic_defect, y, z).total_degree() == 8


# Radial sextic: C has a full factor R and is rank zero modulo R.
R = x**2 + y**2 + z**2
h6_radial = R**3 / 30
s4_radial = R**2
Cr = hessian(h6_radial)
dr = sp.Matrix([sp.diff(s4_radial, variable) for variable in variables])
Deltar = sp.factor(Cr.det())
normr = sp.factor((dr.T * adjugate(Cr) * dr)[0])
radial_matrix_formula = R / 5 * (
    R * sp.eye(3) + 4 * sp.Matrix(variables) * sp.Matrix(variables).T
)
assert all(
    sp.expand(value) == 0 for value in Cr - radial_matrix_formula
)
assert_zero(Deltar - R**6 / 25)
assert_zero(normr - Deltar * 16 * R)
assert all(sp.cancel(entry / R).is_polynomial(x, y, z) for entry in Cr)


# The quintic normalization table: use the largest possible local pole.
def defect_lower_bound(component_degree: int, rho: int, r: int = 3) -> int:
    return max(0, 2 * component_degree * (r + 1 - rho))


assert defect_lower_bound(1, 1) == 6  # line, multiplicity 2 or 3
assert defect_lower_bound(1, 2) == 4  # line, multiplicity 4 or 5
assert defect_lower_bound(1, 3) == 2  # line, multiplicity 6 or 7
assert defect_lower_bound(1, 4) == 0  # line, multiplicity 8 or 9
assert defect_lower_bound(2, 2) == 8  # conic, multiplicity 2 or 3
assert defect_lower_bound(2, 4) == 0  # conic, multiplicity 4
assert defect_lower_bound(3, 3) == 6  # cubic, multiplicity 2 or 3
assert defect_lower_bound(4, 4) == 0  # quartic, multiplicity 2


# Degree-nine determinant plus P^2 | Delta forces rho <= 4.  A clean
# component forces rho >= r+1=4, hence Delta=P^2*ell.
for rho in range(1, 5):
    assert 2 * rho <= 9
assert not 2 * 5 <= 9
assert 2 * 4 + 1 == 9


print("PASS: verified Fermat quintic and sextic Schur modules")
print("PASS: verified Fermat corank-two defect degrees 6 and 8")
print("PASS: verified the radial R^6 lower-Smith calibration")
print("PASS: verified the quintic normalization defect table")
print("THEOREM: a clean quintic Schur survivor has det(C)=P^2*ell, deg(P)=4")
print("SCOPE: normalization and duality steps are written proofs, not CAS checks")

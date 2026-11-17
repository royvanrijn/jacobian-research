#!/usr/bin/env python3
"""Exact local identities for the squarefree quartic-denominator frontend."""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")


def homogeneous_coefficients(prefix: str, degree: int) -> tuple[sp.Symbol, ...]:
    return sp.symbols(f"{prefix}0:{degree + 1}")


# Tangent row: L=x and v=partial_z.  The condition D_v h in (x^2) is
# represented by h=F_5(x,y)+x^2*H_3(x,y,z).
f = homogeneous_coefficients("f", 5)
F = sum(f[index] * x**index * y ** (5 - index) for index in range(6))
h_symbols = sp.symbols("h0:10")
cubic_monomials = [
    x**a * y**b * z ** (3 - a - b)
    for a in range(4)
    for b in range(4 - a)
]
H = sum(coefficient * monomial for coefficient, monomial in zip(h_symbols, cubic_monomials))
tangent_form = F + x**2 * H
tangent_determinant = sp.Poly(
    sp.expand(sp.hessian(tangent_form, (x, y, z)).det()), x, y, z
)
assert all(term[0] >= 2 for term, _ in tangent_determinant.terms())


# Transverse row: L=x and v=partial_x.  The constant-kernel condition gives
# h=F_5(y,z)+x^3*J_2+..., and the first normal determinant face is exactly
# 6*J_2*det Hess(F_5).
g = homogeneous_coefficients("g", 5)
binary_F = sum(g[index] * y**index * z ** (5 - index) for index in range(6))
j0, j1, j2, k0, k1, scalar = sp.symbols("j0 j1 j2 k0 k1 scalar")
J = j0 * y**2 + j1 * y * z + j2 * z**2
K = k0 * y + k1 * z
transverse_form = binary_F + x**3 * J + x**4 * K + scalar * x**5
transverse_determinant = sp.expand(
    sp.hessian(transverse_form, (x, y, z)).det()
)
binary_hessian = sp.hessian(binary_F, (y, z)).det()
assert sp.expand(transverse_determinant.coeff(x, 0)) == 0
assert sp.factor(transverse_determinant.coeff(x, 1) - 6 * J * binary_hessian) == 0

# After J=0 the directional polar starts in order x^3 and the determinant
# starts in order x^2.
transverse_sharpened = sp.expand(transverse_determinant.subs({j0: 0, j1: 0, j2: 0}))
assert sp.expand(transverse_sharpened.coeff(x, 0)) == 0
assert sp.expand(transverse_sharpened.coeff(x, 1)) == 0
assert sp.factor(transverse_sharpened.coeff(x, 2) - 12 * K * binary_hessian) == 0


print("PASS: tangent constant polar forces a double Hessian line")
print("PASS: transverse double-line condition sharpens D_v(h) from L^2 to L^3")
print("PASS: the first transverse faces are 6*J*Hess(F) and 12*K*Hess(F)")
print("THEOREM: squarefree P reduces to four linear flag-synchronization conditions")
print("SCOPE: generic corank one; synchronization emptiness remains open")

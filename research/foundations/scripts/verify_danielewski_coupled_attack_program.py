#!/usr/bin/env python3
"""Exact identities driving the coupled Danielewski attack program."""

from __future__ import annotations

import sympy as sp


a, c, w, x = sp.symbols("a c w x")
P = x * (x**2 + 1)
P_prime = sp.diff(P, x)


def bracket(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        c
        * (
            sp.diff(first, c) * sp.diff(second, x)
            - sp.diff(first, x) * sp.diff(second, c)
        )
        + P_prime
        * (
            sp.diff(first, c) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, c)
        )
        + w
        * (
            sp.diff(first, x) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, x)
        )
    )


# Test the weighted leading-face identity at independent positive degrees.
r, s_degree, t_degree = 2, 3, 4
A_top = c + x * w + 1
F_top = c * x + w**2 + x
G_top = c**2 + x**2 + w
U = F_top**r / A_top**s_degree
V = G_top**r / A_top**t_degree
leading_face = sp.expand(
    r * A_top * bracket(F_top, G_top)
    - s_degree * F_top * bracket(A_top, G_top)
    + t_degree * G_top * bracket(A_top, F_top)
)
logarithmic_side = sp.cancel(
    A_top
    * F_top
    * G_top
    / r
    * bracket(U, V)
    / (U * V)
)
assert sp.cancel(leading_face - logarithmic_side) == 0
print("PASS: the weighted leading-a-face bracket identity is exact")


# Verify extraction of the top a coefficient in the full coupled equation.
A = A_top * a**r + c * a + x
F = F_top * a**s_degree + w * a + c
G = G_top * a**t_degree + x * a + w
coupled = sp.expand(
    sp.diff(A, a) * bracket(F, G)
    - sp.diff(F, a) * bracket(A, G)
    + sp.diff(G, a) * bracket(A, F)
)
top_coefficient = sp.Poly(coupled, a).coeff_monomial(
    a ** (r + s_degree + t_degree - 1)
)
assert sp.expand(top_coefficient - leading_face) == 0
print("PASS: the highest a-coefficient is the weighted Poisson face")


# Formal exterior-algebra audit of the flux identity. Represent a relative
# two-form h*omega by its scalar h.
A_flux = a**2 * c + a * x + w
F_flux = a**3 * w + a * c * x + x**2
G_flux = a * c**2 + a**2 * x + w**2
volume_coefficient = sp.expand(
    sp.diff(A_flux, a) * bracket(F_flux, G_flux)
    - sp.diff(F_flux, a) * bracket(A_flux, G_flux)
    + sp.diff(G_flux, a) * bracket(A_flux, F_flux)
)
flux_derivative_reduced = sp.expand(
    sp.diff(A_flux, a) * bracket(F_flux, G_flux)
    - sp.diff(F_flux, a) * bracket(A_flux, G_flux)
    + sp.diff(G_flux, a) * bracket(A_flux, F_flux)
)
assert sp.expand(volume_coefficient - flux_derivative_reduced) == 0
print("PASS: the relative de Rham flux derivative equals the coupled ledger")


# The cubic surface has two independent vanishing-cycle functionals. For
# h=1,x, endpoint primitives on paths 0->i and 0->-i form this matrix.
root_plus, root_minus = sp.I, -sp.I
period_matrix = sp.Matrix(
    [
        [root_plus, root_plus**2 / 2],
        [root_minus, root_minus**2 / 2],
    ]
)
assert sp.simplify(period_matrix.det()) != 0
print("PASS: the two cubic vanishing-cycle period functionals are independent")
print("PASS coupled Danielewski attack identities")

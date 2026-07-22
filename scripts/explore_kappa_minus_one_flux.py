#!/usr/bin/env python3
"""Explore the missing kappa=-1 chart one shear monomial at a time."""

import argparse
import sympy as sp


parser = argparse.ArgumentParser()
parser.add_argument("--shear-degree", type=int, default=2)
args = parser.parse_args()

w, X, Q, Z, tau, coefficient = sp.symbols("w X Q Z tau coefficient")
rho = sp.symbols("rho")

H = sp.expand(
    w**2
    * (w - 1)
    * (
        tau * w**2
        + (sp.Rational(3, 2) - 2 * tau) * w
        + tau
        - sp.Rational(5, 2)
    )
)
p = sp.diff(H, w)
q = sp.expand(w * p - H)

psi = coefficient * Q ** args.shear_degree
W = Z + psi
source_u = 1 + X * W
source_gamma = 1 + X**2 * Q
marked = sp.expand(source_u * source_gamma)
S = sp.cancel(
    (source_u + q.subs(w, marked) / source_gamma**2) / (2 * X**2)
)
T = sp.cancel((1 + p.subs(w, marked) / source_gamma) / X)
R = 2 * X * source_gamma

SX, SQ, SZ = (sp.diff(S, variable) for variable in (X, Q, Z))
TX, TQ, TZ = (sp.diff(T, variable) for variable in (X, Q, Z))
w_family = sp.Matrix(
    [SZ * TQ - SQ * TZ, SX * TZ - SZ * TX, SQ * TX - SX * TQ]
).applyfunc(sp.cancel)
w_E = sp.Matrix(
    [(1 - 3 * X**2 * Q) / 2, 9 * X * Q**2 / 2, -6 * X * Q * Z]
)
difference = (w_family - w_E).applyfunc(sp.cancel)


def hamiltonian(f):
    return sp.Matrix(
        [
            -2 * X**3 * sp.diff(f, Z),
            (2 + 6 * X**2 * Q) * sp.diff(f, Z),
            2 * X**3 * sp.diff(f, X)
            - (2 + 6 * X**2 * Q) * sp.diff(f, Q),
        ]
    )


f_zero = sp.integrate(sp.cancel(difference[0] / (-2 * X**3)), Z)
residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

# At fixed rho=R, the remaining Hamiltonian vector field is
# 2*X^3*d/dX.  The transformed integrand is a finite Laurent polynomial, so
# integrate it termwise rather than invoking a general rational integrator.
Q_at_fixed_R = (rho - 2 * X) / (2 * X**3)
integrand = sp.cancel(
    residual[2].subs(Q, Q_at_fixed_R) / (2 * X**3)
)
numerator, denominator = sp.together(integrand).as_numer_denom()
denominator_degree = sp.Poly(denominator, X).degree()
denominator_unit = sp.factor(denominator / X**denominator_degree)
numerator_in_X = sp.Poly(numerator, X)
assert numerator_in_X.coeff_monomial(X ** (denominator_degree - 1)) == 0
h_fixed_R = sum(
    term_coefficient
    * X ** (monomial[0] - denominator_degree + 1)
    / (
        denominator_unit
        * (monomial[0] - denominator_degree + 1)
    )
    for monomial, term_coefficient in numerator_in_X.terms()
)
h = sp.cancel(h_fixed_R.subs(rho, R))
assert sp.cancel(residual[2] - hamiltonian(h)[2]) == 0
f = sp.cancel(f_zero + h)

f_numerator, f_denominator = sp.together(f).as_numer_denom()
negative_order = sp.Poly(f_denominator, X).degree()
f_numerator_in_X = sp.Poly(f_numerator, X)

print("shear degree =", args.shear_degree)
print("denominator =", sp.factor(f_denominator))
for exponent in range(negative_order):
    residue = sp.factor(f_numerator_in_X.coeff_monomial(X**exponent))
    if residue != 0:
        print(f"residue numerator [X^{exponent}] =", residue)

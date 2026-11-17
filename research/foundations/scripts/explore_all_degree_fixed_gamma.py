#!/usr/bin/env python3
"""Probe the fixed-kappa=-9 rank-two descent in a chosen degree.

This is an exploratory exact calculation, not an all-parameter certificate.
For degree N it uses the normalized seed

    H=w^2(w-1)P,  P=w^(N-3)-(N-1/2)w+(N-5/2),

so that P(1)=-1, P'(1)=-5/2 and hence H''(1)=-9.
"""

import argparse
import sympy as sp


parser = argparse.ArgumentParser()
parser.add_argument("degree", type=int)
args = parser.parse_args()
if args.degree < 5:
    parser.error("degree must be at least five")

w, X, Q, Z, s2 = sp.symbols("w X Q Z s2")
V, rho = sp.symbols("V rho")
d = args.degree - 3
P = w**d - sp.Rational(2 * d + 5, 2) * w + sp.Rational(2 * d + 1, 2)
H = sp.expand(w**2 * (w - 1) * P)
p = sp.diff(H, w)
q = sp.expand(w * p - H)
assert sp.degree(H, w) == args.degree
assert H.subs(w, 1) == 0
assert p.subs(w, 1) == -1
assert sp.diff(H, w, 2).subs(w, 1) == -9

a = -sp.Rational(8, 7)
W = Z + s2 * Q**2
Y = Q - X * W / 3
source_v = -3 * X * Y / (2 * a)
source_gamma = 1 - 3 * X * Q / 2
source_u = 1 + source_v
marked = sp.expand(source_u * source_gamma)
S = sp.cancel(
    -2
    * a
    * (source_u + q.subs(w, marked) / source_gamma**2)
    / (3 * X**2)
)
T = sp.cancel((1 + p.subs(w, marked) / source_gamma) / X)
R = sp.expand(2 * X * source_gamma)


def bracket(left, right):
    return sp.cancel(
        -3 * X**2 * (sp.diff(left, X) * sp.diff(right, Z) - sp.diff(left, Z) * sp.diff(right, X))
        + (6 * X * Q - 2) * (sp.diff(left, Q) * sp.diff(right, Z) - sp.diff(left, Z) * sp.diff(right, Q))
    )


assert sp.factor(sp.Matrix([S, T, R]).jacobian((X, Q, Z)).det()) == -1
assert bracket(S, T) == 1
assert bracket(R, S) == bracket(R, T) == 0

SX, SQ, SZ = (sp.diff(S, variable) for variable in (X, Q, Z))
TX, TQ, TZ = (sp.diff(T, variable) for variable in (X, Q, Z))
w_family = sp.Matrix(
    [SZ * TQ - SQ * TZ, SX * TZ - SZ * TX, SQ * TX - SX * TQ]
).applyfunc(sp.cancel)
w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
difference = (w_family - w_E).applyfunc(sp.cancel)


def hamiltonian(f):
    return sp.Matrix(
        [
            3 * X**2 * sp.diff(f, Z),
            (2 - 6 * X * Q) * sp.diff(f, Z),
            -3 * X**2 * sp.diff(f, X) + (6 * X * Q - 2) * sp.diff(f, Q),
        ]
    )


integrand_Z = sp.cancel(difference[0] / (3 * X**2))
num_Z, den_Z = sp.together(integrand_Z).as_numer_denom()
f_zero = sum(
    coefficient * Z ** (monomial[0] + 1) / (den_Z * (monomial[0] + 1))
    for monomial, coefficient in sp.Poly(num_Z, Z).terms()
)
residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

integrand_V = sp.cancel(
    residual[2].subs({X: 1 / V, Q: V * (2 - rho * V) / 3}) / 3
)
num_V, den_V = sp.together(integrand_V).as_numer_denom()
assert sp.Poly(den_V, V).degree() == 0
h_fixed_R = sum(
    coefficient * V ** (monomial[0] + 1) / (den_V * (monomial[0] + 1))
    for monomial, coefficient in sp.Poly(num_V, V).terms()
)
h = sp.cancel(h_fixed_R.subs({V: 1 / X, rho: R}))
assert sp.cancel(residual[2] - hamiltonian(h)[2]) == 0
f = sp.cancel(f_zero + h)

numerator, denominator = sp.together(f).as_numer_denom()
denominator_poly = sp.Poly(denominator, X)
pole_order = denominator_poly.degree()
assert sp.factor(denominator / X**pole_order).has(X) is False
poly_X = sp.Poly(numerator, X)
residues = [
    sp.factor(poly_X.coeff_monomial(X**exponent))
    for exponent in range(pole_order)
]
nonzero_residues = [residue for residue in residues if residue != 0]
candidate_solutions = sp.solve(nonzero_residues[0], s2)
candidate = candidate_solutions[0] if len(candidate_solutions) == 1 else None
completes = candidate is not None and all(
    sp.factor(residue.subs(s2, candidate)) == 0 for residue in residues
)

print(f"degree={args.degree} pole_order={pole_order}")
for exponent, residue in enumerate(residues):
    if residue != 0:
        print(f"residue numerator [X^{exponent}] = {residue}")
print("candidate s2 =", candidate)
print("single Q^2 shear completes sample =", completes)

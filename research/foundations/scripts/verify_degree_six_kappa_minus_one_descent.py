#!/usr/bin/env python3
"""Exact degree-six rank-two descent certificate on kappa=-1."""

import sympy as sp


w, X, Q, Z = sp.symbols("w X Q Z")
sigma, tau, c = sp.symbols("sigma tau c")
rho = sp.symbols("rho")

kappa = sp.Integer(-1)
H = sp.expand(
    w**2
    * (w - 1)
    * (
        sigma * w**3
        + tau * w**2
        + (kappa / 2 - 3 * sigma - 2 * tau + 2) * w
        - kappa / 2
        + 2 * sigma
        + tau
        - 3
    )
)
p = sp.diff(H, w)
q = sp.expand(w * p - H)
assert p.subs(w, 1) == -1
assert sp.diff(H, w, 2).subs(w, 1) == kappa
assert sp.Poly(H, w).LC() == sigma

W = Z + c * X * Q
source_u = 1 + X * W
source_gamma = 1 + X**2 * Q
marked = sp.expand(source_u * source_gamma)
S = sp.cancel(
    (source_u + q.subs(w, marked) / source_gamma**2) / (2 * X**2)
)
T = sp.cancel((1 + p.subs(w, marked) / source_gamma) / X)
R = sp.expand(2 * X * source_gamma)
assert R == 2 * X + 2 * X**3 * Q
assert all(
    not ({X, Q, Z} & sp.denom(component).free_symbols)
    for component in (S, T, R)
)

base_jacobian = sp.Matrix([S, T, R]).jacobian((X, Q, Z))
assert sp.factor(base_jacobian.det()) == -1


def quotient_bracket(left, right):
    return sp.cancel(
        2
        * X**3
        * (
            sp.diff(left, X) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, X)
        )
        - (2 + 6 * X**2 * Q)
        * (
            sp.diff(left, Q) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, Q)
        )
    )


assert quotient_bracket(S, T) == 1
assert quotient_bracket(R, S) == quotient_bracket(R, T) == 0

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


integrand_Z = sp.cancel(difference[0] / (-2 * X**3))
num_Z, den_Z = sp.together(integrand_Z).as_numer_denom()
poly_Z = sp.Poly(num_Z, Z)
f_zero = sum(
    coefficient * Z ** (monomial[0] + 1) / (den_Z * (monomial[0] + 1))
    for monomial, coefficient in poly_Z.terms()
)
residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

Q_at_fixed_R = (rho - 2 * X) / (2 * X**3)
integrand_X = sp.cancel(residual[2].subs(Q, Q_at_fixed_R) / (2 * X**3))
num_X, den_X = sp.together(integrand_X).as_numer_denom()
denominator_degree = sp.Poly(den_X, X).degree()
denominator_unit = sp.factor(den_X / X**denominator_degree)
poly_integrand = sp.Poly(num_X, X)
assert poly_integrand.coeff_monomial(X ** (denominator_degree - 1)) == 0
h_fixed_R = sum(
    coefficient
    * X ** (monomial[0] - denominator_degree + 1)
    / (denominator_unit * (monomial[0] - denominator_degree + 1))
    for monomial, coefficient in poly_integrand.terms()
)
h = sp.cancel(h_fixed_R.subs(rho, R))
assert sp.cancel(residual[2] - hamiltonian(h)[2]) == 0
f = sp.cancel(f_zero + h)

numerator, denominator = sp.together(f).as_numer_denom()
assert sp.factor(denominator - 18480 * X**4) == 0
obstruction = sp.expand(
    1155 * c
    - 251 * sigma**2
    - 209 * sigma * tau
    + 858 * sigma
    - 44 * tau**2
    + 330 * tau
    + 396
)
poly_X = sp.Poly(numerator, X)
assert sp.factor(poly_X.coeff_monomial(1) - obstruction) == 0
assert poly_X.coeff_monomial(X) == 0
assert sp.factor(poly_X.coeff_monomial(X**2) + 4 * Q * obstruction) == 0
assert poly_X.coeff_monomial(X**3) == 0

candidate = sp.factor(
    (
        251 * sigma**2
        + 209 * sigma * tau
        - 858 * sigma
        + 44 * tau**2
        - 330 * tau
        - 396
    )
    / 1155
)
assert sp.factor(obstruction.subs(c, candidate)) == 0

W_back = sp.symbols("W_back")
assert sp.expand(
    (Z + candidate * X * Q).subs(Z, W_back - candidate * X * Q) - W_back
) == 0

print("PASS: full normalized degree-six kappa=-1 seed divisor")
print("PASS: complete principal part has one scalar obstruction")
print("PASS: unique X*Q shear gives a polynomial rank-two Hamiltonian")

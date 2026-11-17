#!/usr/bin/env python3
"""Exact rank-two descent certificate for the generic degree-six seed space.

The calculation is performed over Q(a,sigma,tau,s2), where

    kappa = -(1+2*a)/(1+a).

Thus this chart is exactly kappa != -1,-2.  The exceptional admissible
divisor kappa=-1 is handled by verify_degree_six_kappa_minus_one_descent.py.
"""

import sympy as sp


w, X, Q, Z = sp.symbols("w X Q Z")
a, sigma, tau, s2 = sp.symbols("a sigma tau s2")
V, rho = sp.symbols("V rho")

kappa = -(1 + 2 * a) / (1 + a)
H = sp.factor(
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
assert sp.factor(H.subs(w, 1)) == 0
assert sp.factor(p.subs(w, 1) + 1) == 0
assert sp.factor(sp.diff(H, w, 2).subs(w, 1) - kappa) == 0
assert sp.factor(sp.Poly(H, w).LC() - sigma) == 0

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
assert R == 2 * X - 3 * X**2 * Q
assert all(
    not ({X, Q, Z} & sp.denom(component).free_symbols)
    for component in (S, T, R)
)

base_jacobian = sp.Matrix([S, T, R]).jacobian((X, Q, Z))
assert sp.factor(base_jacobian.det()) == -1


def quotient_bracket(left, right):
    return sp.cancel(
        -3
        * X**2
        * (
            sp.diff(left, X) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, X)
        )
        + (6 * X * Q - 2)
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
w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
difference = (w_family - w_E).applyfunc(sp.cancel)


def hamiltonian(f):
    return sp.Matrix(
        [
            3 * X**2 * sp.diff(f, Z),
            (2 - 6 * X * Q) * sp.diff(f, Z),
            -3 * X**2 * sp.diff(f, X)
            + (6 * X * Q - 2) * sp.diff(f, Q),
        ]
    )


# Integrate the Z-component termwise.  This avoids heuristic rational
# integration and keeps every identity in the exact rational function field.
integrand_Z = sp.cancel(difference[0] / (3 * X**2))
num_Z, den_Z = sp.together(integrand_Z).as_numer_denom()
poly_Z = sp.Poly(num_Z, Z)
f_zero = sum(
    coefficient * Z ** (monomial[0] + 1) / (den_Z * (monomial[0] + 1))
    for monomial, coefficient in poly_Z.terms()
)
residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

# At fixed rho=R, V=X^-1 turns the remaining Hamiltonian field into 3*d/dV.
integrand_V = sp.cancel(
    residual[2].subs({X: 1 / V, Q: V * (2 - rho * V) / 3}) / 3
)
num_V, den_V = sp.together(integrand_V).as_numer_denom()
assert sp.Poly(den_V, V).degree() == 0
poly_V = sp.Poly(num_V, V)
h_fixed_R = sum(
    coefficient * V ** (monomial[0] + 1) / (den_V * (monomial[0] + 1))
    for monomial, coefficient in poly_V.terms()
)
h = sp.cancel(h_fixed_R.subs({V: 1 / X, rho: R}))
assert sp.cancel(residual[2] - hamiltonian(h)[2]) == 0
f = sp.cancel(f_zero + h)

numerator, denominator = sp.together(f).as_numer_denom()
assert sp.factor(
    denominator - 681246720 * X**4 * a**10 * (a + 1) ** 2
) == 0

obstruction = sp.expand(
    308 * (a + 1) ** 2 * s2
    + a * (a + 1) ** 2 * (753 * sigma**2 + 627 * sigma * tau + 132 * tau**2)
    - a
    * (a + 1)
    * ((2079 * a + 2574) * sigma + (792 * a + 990) * tau)
    - 99 * (24 * a**3 + 72 * a**2 + 82 * a + 35)
)
poly_X = sp.Poly(numerator, X)
expected = (
    16384 * a**10 * obstruction,
    98304 * Q * a**10 * obstruction,
    368640 * Q**2 * a**10 * obstruction,
    1105920 * Q**3 * a**10 * obstruction,
)
for exponent, coefficient in enumerate(expected):
    assert sp.factor(poly_X.coeff_monomial(X**exponent) - coefficient) == 0

candidate = sp.factor(
    (
        -a
        * (a + 1) ** 2
        * (753 * sigma**2 + 627 * sigma * tau + 132 * tau**2)
        + a
        * (a + 1)
        * ((2079 * a + 2574) * sigma + (792 * a + 990) * tau)
        + 99 * (24 * a**3 + 72 * a**2 + 82 * a + 35)
    )
    / (308 * (a + 1) ** 2)
)
assert sp.factor(obstruction.subs(s2, candidate)) == 0
assert all(
    sp.factor(poly_X.coeff_monomial(X**exponent).subs(s2, candidate)) == 0
    for exponent in range(4)
)

fixed_gamma_candidate = sp.factor(
    (
        6024 * sigma**2
        + 5016 * sigma * tau
        + 11088 * sigma
        + 1056 * tau**2
        + 4752 * tau
        - 16929
    )
    / 2156
)
assert sp.factor(candidate.subs(a, -sp.Rational(8, 7)) - fixed_gamma_candidate) == 0

# The shear is triangular and polynomially invertible on the parameter chart.
W_back = sp.symbols("W_back")
assert sp.expand(
    (Z + candidate * Q**2).subs(Z, W_back - candidate * Q**2) - W_back
) == 0

print("PASS: full generic normalized degree-six seed space")
print("PASS: complete principal part has one scalar obstruction")
print("PASS: unique Q^2 shear gives a polynomial rank-two Hamiltonian")
print("PASS: generic chart is kappa != -1,-2")

#!/usr/bin/env python3
"""Exact rank-two descent certificate for the degree-six kappa=-9 slice."""

import sympy as sp

w, X, Q, Z = sp.symbols("w X Q Z")
sigma, tau, s2 = sp.symbols("sigma tau s2")
V, rho = sp.symbols("V rho")

kappa = sp.Integer(-9)
a = -sp.Rational(8, 7)
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

# Exact witness that the fixed-kappa surface meets the full-cover faithful,
# exact-double-zero boundary-clean open used by the stable-moduli theorem.
witness = sp.factor(H.subs({sigma: 1, tau: 0}))
assert witness == w**2 * (w - 1) * (2 * w**3 - 11 * w + 7) / 2
assert sp.gcd(witness, sp.diff(witness, w)) == w
assert sp.diff(witness, w, 2).subs(w, 0) == -7
assert sp.discriminant(sp.cancel(witness / w**2), w) == sp.Rational(1339, 4)
assert sp.discriminant(sp.diff(witness, w, 2), w) == 19184247360

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
R = 2 * X * source_gamma

base_jacobian = sp.Matrix([S, T, R]).jacobian((X, Q, Z))
assert sp.factor(base_jacobian.det()) == -1


def quotient_bracket(left, right):
    return sp.cancel(
        -3 * X**2
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


integrand_Z = sp.cancel(difference[0] / (3 * X**2))
num_Z, den_Z = sp.together(integrand_Z).as_numer_denom()
poly_Z = sp.Poly(num_Z, Z)
f_zero = sum(
    coefficient * Z ** (monomial[0] + 1) / (den_Z * (monomial[0] + 1))
    for monomial, coefficient in poly_Z.terms()
)
residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

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
assert denominator == 2560190835043860480 * X**4
obstruction = (
    2156 * s2
    - 6024 * sigma**2
    - 5016 * sigma * tau
    - 11088 * sigma
    - 1056 * tau**2
    - 4752 * tau
    + 16929
)
poly_X = sp.Poly(numerator, X)
expected = (
    8796093022208 * obstruction,
    52776558133248 * Q * obstruction,
    197912092999680 * Q**2 * obstruction,
    593736278999040 * Q**3 * obstruction,
)
for exponent, coefficient in enumerate(expected):
    assert sp.factor(poly_X.coeff_monomial(X**exponent) - coefficient) == 0

candidate = sp.factor(
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
assert sp.factor(obstruction.subs(s2, candidate)) == 0

print("PASS: full normalized degree-six kappa=-9 seed slice")
print("PASS: complete principal part has one scalar obstruction")
print("PASS: unique Q^2 shear gives a polynomial rank-two Hamiltonian")
print("PASS: fixed-gamma surface meets the full-cover faithful stable open")

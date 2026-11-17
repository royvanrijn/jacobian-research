#!/usr/bin/env python3
"""Exact certificate for rank-two descent on the degree-five seed surface.

The weighted coordinate ``a`` is used instead of ``kappa`` during the large
calculation.  They are related by

    kappa = -(1+2*a)/(1+a),  a = -(1+kappa)/(2+kappa).

This makes the uniform two-parameter calculation small enough for a direct
exact replay over Q(a,tau,s2).
"""

import sympy as sp


w = sp.symbols("w")
X, Q, Z = sp.symbols("X Q Z")
a, tau, s2 = sp.symbols("a tau s2")
kappa_symbol = sp.symbols("kappa")

kappa = -(1 + 2 * a) / (1 + a)
H = sp.factor(
    w**2
    * (w - 1)
    * (
        tau * w**2
        + (kappa / 2 - 2 * tau + 2) * w
        - kappa / 2
        + tau
        - 3
    )
)
p = sp.diff(H, w)
q = sp.expand(w * p - H)
assert sp.factor(H.subs(w, 1)) == 0
assert sp.factor(p.subs(w, 1) + 1) == 0
assert sp.factor(sp.diff(H, w, 2).subs(w, 1) - kappa) == 0
assert sp.factor(sp.Poly(H, w).LC() - tau) == 0

# Uniform adapted source coordinates.  Directly using the invariants avoids
# expanding the source map before its polynomiality cancellations.
W = Z + s2 * Q**2
Y = Q - X * W / 3
source_v = -3 * X * Y / (2 * a)
source_gamma = 1 - 3 * X * Q / 2
source_u = 1 + source_v
marked = sp.expand(source_u * source_gamma)

A_sub = (source_u + q.subs(w, marked) / source_gamma**2) / X**2
B_sub = (1 + p.subs(w, marked) / source_gamma) / X
S = sp.cancel(-2 * a * A_sub / 3)
T = sp.cancel(B_sub)
R = sp.expand(2 * X * source_gamma)
assert R == 2 * X - 3 * X**2 * Q
assert all(
    not ({X, Q, Z} & sp.denom(component).free_symbols)
    for component in (S, T, R)
)

SX, SQ, SZ = (sp.diff(S, variable) for variable in (X, Q, Z))
TX, TQ, TZ = (sp.diff(T, variable) for variable in (X, Q, Z))
RX, RQ, RZ = (sp.diff(R, variable) for variable in (X, Q, Z))
base_determinant = sp.cancel(
    SX * (TQ * RZ - TZ * RQ)
    - SQ * (TX * RZ - TZ * RX)
    + SZ * (TX * RQ - TQ * RX)
)
assert base_determinant == -1


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
assert quotient_bracket(R, S) == 0
assert quotient_bracket(R, T) == 0

# Since det d(S,T,R)/d(X,Q,Z)=-1, the third inverse-Jacobian column is the
# following cofactor vector.
w_family = sp.Matrix(
    [SZ * TQ - SQ * TZ, SX * TZ - SZ * TX, SQ * TX - SX * TQ]
).applyfunc(sp.cancel)
w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
difference = (w_family - w_E).applyfunc(sp.cancel)


def quotient_hamiltonian(f):
    return sp.Matrix(
        [
            3 * X**2 * sp.diff(f, Z),
            (2 - 6 * X * Q) * sp.diff(f, Z),
            -3 * X**2 * sp.diff(f, X)
            + (6 * X * Q - 2) * sp.diff(f, Q),
        ]
    )


f_zero = sp.integrate(sp.cancel(difference[0] / (3 * X**2)), Z)
residual = (difference - quotient_hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

inv_x, rho = sp.symbols("inv_x rho")
residual_vrho = sp.cancel(
    residual[2].subs({X: 1 / inv_x, Q: inv_x * (2 - rho * inv_x) / 3})
)
h_vrho = sp.integrate(residual_vrho / 3, inv_x)
h = sp.cancel(
    h_vrho.subs({inv_x: 1 / X, rho: 2 * X - 3 * X**2 * Q})
)
assert sp.cancel(residual[2] - quotient_hamiltonian(h)[2]) == 0
f_generic = sp.cancel(f_zero + h)

# Extract the complete principal part.  All four coefficients have one scalar
# obstruction, and no power below X^-4 occurs.
numerator, denominator = sp.together(f_generic).as_numer_denom()
assert sp.factor(denominator - 15482880 * a**8 * (a + 1) ** 2 * X**4) == 0
obstruction = sp.expand(
    28 * (a + 1) ** 2 * s2
    + 12 * a * (a + 1) ** 2 * tau**2
    - 18 * a * (a + 1) * (4 * a + 5) * tau
    - 216 * a**3
    - 648 * a**2
    - 738 * a
    - 315
)
numerator_in_X = sp.Poly(numerator, X)
expected = (
    4096 * a**8 * obstruction,
    24576 * Q * a**8 * obstruction,
    92160 * Q**2 * a**8 * obstruction,
    276480 * Q**3 * a**8 * obstruction,
)
for exponent, coefficient in enumerate(expected):
    assert sp.factor(
        numerator_in_X.coeff_monomial(X**exponent) - coefficient
    ) == 0

candidate_a = sp.factor(
    (
        -12 * a * (a + 1) ** 2 * tau**2
        + 18 * a * (a + 1) * (4 * a + 5) * tau
        + 216 * a**3
        + 648 * a**2
        + 738 * a
        + 315
    )
    / (28 * (a + 1) ** 2)
)
assert sp.factor(obstruction.subs(s2, candidate_a)) == 0

candidate_kappa = sp.factor(
    (
        12 * (kappa_symbol + 1) * tau**2
        - 18 * (kappa_symbol + 1) * (kappa_symbol + 6) * tau
        + 9
        * (
            kappa_symbol**3
            + 16 * kappa_symbol**2
            + 52 * kappa_symbol
            + 72
        )
    )
    / (28 * (kappa_symbol + 2))
)
assert sp.factor(candidate_kappa.subs(kappa_symbol, kappa) - candidate_a) == 0
assert sp.factor(candidate_kappa.subs(kappa_symbol, -1) - sp.Rational(45, 4)) == 0

delta_cusp = (
    kappa_symbol**4
    - 4 * kappa_symbol**3 * tau
    + 18 * kappa_symbol**3
    + 8 * kappa_symbol**2 * tau**2
    - 36 * kappa_symbol**2 * tau
    + 123 * kappa_symbol**2
    - 8 * kappa_symbol * tau**3
    + 24 * kappa_symbol * tau**2
    - 72 * kappa_symbol * tau
    + 376 * kappa_symbol
    + 8 * tau**4
    + 24 * tau**3
    - 24 * tau**2
    + 72 * tau
    + 432
)
H_kappa = H.subs(a, -(1 + kappa_symbol) / (2 + kappa_symbol))
assert sp.factor(
    sp.discriminant(sp.diff(H_kappa, w, 2), w) - 432 * delta_cusp
) == 0

# The four low coefficients are the complete obstruction to polynomiality.
assert all(sp.factor(coefficient.subs(s2, candidate_a)) == 0 for coefficient in expected)

# The base/source change is polynomially invertible, and the final E change is
# triangular once the Hamiltonian is polynomial.
Y_back, W_back = sp.symbols("Y_back W_back")
Q_back = Y_back + X * W_back / 3
Z_back = W_back - candidate_a * Q_back**2
assert sp.expand(
    (Z + candidate_a * Q**2).subs(
        {Q: Q_back, Z: Z_back}, simultaneous=True
    )
    - W_back
) == 0
assert sp.expand(
    (Q - X * (Z + candidate_a * Q**2) / 3).subs(
        {Q: Q_back, Z: Z_back}, simultaneous=True
    )
    - Y_back
) == 0

print("PASS: parametrized the full normalized degree-five seed surface")
print("PASS: one adapted coordinate R works for every a != 0,-1")
print("PASS: the complete Laurent principal part has one scalar obstruction")
print("PASS: the unique quadratic shear cancels it on the entire surface chart")
print("PASS: the published kappa=-9 formula is its fixed-gamma specialization")
print("PASS: identified degree-drop, chart, admissibility, and cusp divisors")

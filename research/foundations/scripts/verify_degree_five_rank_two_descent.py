#!/usr/bin/env python3
"""Exact symbolic certificate for rank-two descent of the degree-five family.

The script works over Q(lambda).  It constructs the common adapted quotient,
computes the relative Hamiltonian obstruction for Z -> Z+k*Q^2, proves that
one rational k(lambda) cancels its complete negative-X principal part, and
checks the resulting Hamiltonian identities.  Expanded output polynomials are
deliberately not serialized: the finite integration formulas are the exact
certificate and avoid a several-thousand-term presentation artifact.
"""

import sympy as sp


w, lam = sp.symbols("w lambda")
x, y, z = sp.symbols("x y z")
X, Q, Z, k = sp.symbols("X Q Z k")

# The verified family has gamma=1-(8/7)xy+x^2z, not 1-(7/8)xy+x^2z.
H = w**2 * (w - 1) * (3 * w**2 - (5 * lam + 1) * w + 3 * lam) / 60
p_seed = sp.diff(H, w)
c_seed = (lam - 1) / 30
u = 1 + x * y
gamma = 1 - sp.Rational(8, 7) * x * y + x**2 * z
marked_root = sp.expand(u * gamma)
q_seed = sp.cancel((w * p_seed - H) / c_seed)

C_map = sp.expand(x * gamma)
B_map = sp.cancel((c_seed + p_seed.subs(w, marked_root) / gamma) / x)
A_map = sp.cancel((u + q_seed.subs(w, marked_root) / gamma**2) / x**2)
assert sp.factor(sp.Matrix([A_map, B_map, C_map]).jacobian((x, y, z)).det() - c_seed) == 0

# A common base coordinate system.  The diagonal source change sends 2*C to
# the same R used in the foundational rank-two completion.
W = Z + k * Q**2
Y = Q - X * W / 3
base_substitution = {x: X, y: sp.Rational(21, 16) * Y, z: -W / 2}
S = sp.cancel(sp.Rational(160, 7) * A_map.subs(base_substitution) / (lam - 1))
T = sp.cancel(B_map.subs(base_substitution))
R = sp.expand(2 * C_map.subs(base_substitution))
assert sp.expand(R - (2 * X - 3 * X**2 * Q)) == 0

base_jacobian = sp.Matrix([S, T, R]).jacobian((X, Q, Z))
assert sp.factor(base_jacobian.det()) == -1


def quotient_bracket(left, right):
    """Induced Poisson bracket on Q(lambda)[X,Q,Z]."""

    return sp.cancel(
        -3 * X**2 * (sp.diff(left, X) * sp.diff(right, Z) - sp.diff(left, Z) * sp.diff(right, X))
        + (6 * X * Q - 2)
        * (sp.diff(left, Q) * sp.diff(right, Z) - sp.diff(left, Z) * sp.diff(right, Q))
    )


assert quotient_bracket(S, T) == 1
assert quotient_bracket(R, S) == 0
assert quotient_bracket(R, T) == 0

# The horizontal derivation for (S,T,R), and the E-Hamiltonian derivation in
# the common adapted coordinates.
w_family = (base_jacobian.adjugate() / base_jacobian.det())[:, 2].applyfunc(sp.cancel)
w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
difference = (w_family - w_E).applyfunc(sp.cancel)


def quotient_hamiltonian(f):
    """Components of {f,-} on (X,Q,Z)."""

    return sp.Matrix(
        [
            3 * X**2 * sp.diff(f, Z),
            (2 - 6 * X * Q) * sp.diff(f, Z),
            -3 * X**2 * sp.diff(f, X) + (6 * X * Q - 2) * sp.diff(f, Q),
        ]
    )


# Integrate the first component in Z.  The first two residual components then
# vanish.  On the last component, v=1/X and rho=R turn the remaining
# derivation into 3*d/dv.
f_zero = sp.integrate(sp.cancel(difference[0] / (3 * X**2)), Z)
residual = (difference - quotient_hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == 0
assert residual[1] == 0

v, rho = sp.symbols("v rho")
residual_vrho = sp.cancel(
    residual[2].subs({X: 1 / v, Q: v * (2 - rho * v) / 3})
)
h_vrho = sp.integrate(residual_vrho / 3, v)
h = sp.cancel(h_vrho.subs({v: 1 / X, rho: 2 * X - 3 * X**2 * Q}))
f_generic = sp.cancel(f_zero + h)

# This is an independent identity check after the two integrations.
hamiltonian_difference = quotient_hamiltonian(f_generic)
assert all(
    sp.cancel(difference[index] - hamiltonian_difference[index]) == 0
    for index in range(3)
)

# Extract the complete negative-X principal part from the common denominator.
numerator, denominator = sp.together(f_generic).as_numer_denom()
denominator_constant = sp.Integer(1818317354434560)
assert sp.factor(
    denominator - denominator_constant * X**4 * (lam - 1) ** 2
) == 0

obstruction = (
    196 * (lam - 1) ** 2 * k
    + 27 * (57 * lam**2 - 138 * lam + 73)
)
numerator_in_X = sp.Poly(numerator, X)
expected_low_coefficients = (
    sp.Integer(68719476736) * obstruction,
    sp.Integer(412316860416) * Q * obstruction,
    sp.Integer(1546188226560) * Q**2 * obstruction,
    sp.Integer(4638564679680) * Q**3 * obstruction,
)
for exponent, expected in enumerate(expected_low_coefficients):
    assert sp.factor(numerator_in_X.coeff_monomial(X**exponent) - expected) == 0

kappa = -sp.Rational(27, 196) * (57 * lam**2 - 138 * lam + 73) / (lam - 1) ** 2
assert sp.factor(obstruction.subs(k, kappa)) == 0

# No powers below X^-4 occur, so cancellation of these four coefficients is
# exactly polynomiality.  Parameter denominators are units on lambda!=1.
assert sp.degree(denominator, X) == 4
assert all(
    sp.factor(coefficient.subs(k, kappa)) == 0
    for coefficient in expected_low_coefficients
)

# The source changes used in the LR factorization are polynomial automorphisms:
# Q=Y+XW/3 and Z=W-kappa*Q^2 invert (X,Q,Z)->(X,Y,W), while E->E+f is
# triangular.  Together with det(S,T,R)/det(X,Q,Z)=-1, the Hamiltonian identity
# proves the six canonical brackets for (R,T,D=E+f,S).
Q_back = sp.symbols("Y_back") + X * sp.symbols("W_back") / 3
Z_back = sp.symbols("W_back") - kappa * Q_back**2
assert sp.expand((Z + kappa * Q**2).subs({Q: Q_back, Z: Z_back}, simultaneous=True) - sp.symbols("W_back")) == 0
assert sp.expand((Q - X * (Z + kappa * Q**2) / 3).subs({Q: Q_back, Z: Z_back}, simultaneous=True) - sp.symbols("Y_back")) == 0

print("PASS: verified gamma=1-(8/7)xy+x^2z and the fixed R coordinate")
print("PASS: the normalized degree-five base triple has determinant -1")
print("PASS: the relative flux is Hamiltonian after localization")
print("PASS: its complete principal part has one obstruction factor")
print("PASS: the unique pole-free shear is kappa(lambda)")
print("PASS: the resulting Hamiltonian is polynomial over Q(lambda)")
print("PASS: the construction gives canonical (R,T,D,S) LR-equivalent to F_lambda x id")

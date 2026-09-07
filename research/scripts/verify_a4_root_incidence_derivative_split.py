#!/usr/bin/env python3
"""Exact quotient-algebra checks for the A4 derivative-unit split."""

import sympy as sp


# ---------------------------------------------------------------------------
# 1. A compact three-coefficient form of the scaled JLY quartic
# ---------------------------------------------------------------------------

A, B, C, T, Q, Omega = sp.symbols("A B C T Q Omega")

P = (
    T**4
    - 6 * A * B * T**2
    - 8 * B**3 * T
    + B**2 * (9 * A**2 - 12 * C * B)
)
dP = sp.diff(P, T)
R = (
    4 * A**3 * B
    - 3 * A**2 * C**2
    - 6 * A * B**2 * C
    + B**4
    + 4 * B * C**3
)
disc = -110592 * B**8 * R

assert sp.factor(sp.discriminant(P, T) - disc) == 0

# The inverse of P' in the rank-four algebra Q(A,B,C)[T]/(P) has this
# compact representative.  Keeping the unreduced discriminant denominator
# is what exposes the two equal orientation factors.
inverse_numerator_reduced = (
    (-A * C + B**2) * T**3
    + 2 * B * (A**2 - B * C) * T**2
    + B * (3 * A**2 * C - 7 * A * B**2 + 4 * B * C**2) * T
    - 6 * B**2 * (A**3 - 2 * A * B * C + B**3)
)
inverse_denominator_reduced = 48 * B**4 * R

assert sp.rem(
    sp.expand(dP * inverse_numerator_reduced - inverse_denominator_reduced),
    P,
    T,
) == 0

adjugate = sp.expand(-2304 * B**4 * inverse_numerator_reduced)
assert sp.rem(sp.expand(dP * adjugate - disc), P, T) == 0


# ---------------------------------------------------------------------------
# 2. Split the square orientation denominator between two coordinates
# ---------------------------------------------------------------------------

primitive = (
    576 * B**4 * (A * C - B**2) * T**4
    - 1536 * B**5 * (A**2 - B * C) * T**3
    - 1152
    * B**5
    * (3 * A**2 * C - 7 * A * B**2 + 4 * B * C**2)
    * T**2
    + 13824 * B**6 * (A**3 - 2 * A * B * C + B**3) * T
)
assert sp.expand(sp.diff(primitive, T) - adjugate) == 0

# On Omega^2=Disc(P), the two coordinates primitive/Omega and Q/Omega
# have relative Jacobian 1/P' in the quotient algebra.
split_pair = sp.Matrix([primitive / Omega, Q / Omega])
split_jacobian = sp.factor(
    split_pair.jacobian((T, Q)).det().subs(Omega**2, disc)
)
assert sp.cancel(
    split_jacobian
    - inverse_numerator_reduced / inverse_denominator_reduced
) == 0

# The first split coordinate still generates the generic root algebra.
# It is enough to show that the power-basis change determinant is not the
# zero polynomial.  At A=C=0, B=1, primitive mod P is 9216*T.
primitive_mod_P = sp.rem(sp.expand(primitive), P, T)
power_columns = []
for exponent in range(4):
    reduced_power = sp.Poly(
        sp.rem(sp.expand(primitive_mod_P**exponent), P, T),
        T,
    )
    power_columns.append([reduced_power.nth(index) for index in range(4)])
power_matrix = sp.Matrix(
    4,
    4,
    lambda row, column: power_columns[column][row],
)
power_determinant = sp.expand(power_matrix.det())
assert power_determinant.subs({A: 0, B: 1, C: 0}) == 9216**6


# ---------------------------------------------------------------------------
# 3. Specialize to the two-parameter generic A4 presentation
# ---------------------------------------------------------------------------

a, b = sp.symbols("a b")

A0 = a**3 - b**3 - 9 * b**2 - 27 * b - 54
B0 = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
C0 = a**3 - b**3 + 27
rho = b**2 + 3 * b + 9
sigma = (
    2 * a**3 * b
    + 3 * a**3
    - 3 * a**2 * b**2
    - 9 * a**2 * b
    - 27 * a**2
    + b**4
    + 6 * b**3
    + 27 * b**2
    + 54 * b
    + 81
)
orientation = 1728 * rho * B0**4 * sigma
jly_substitution = {A: A0, B: B0, C: C0}

# The discriminant becomes the literal square orientation^2.
assert sp.factor(R.subs(jly_substitution) + 27 * rho**2 * sigma**2) == 0
assert sp.factor(disc.subs(jly_substitution) - orientation**2) == 0

# The compact inverse numerator acquires one factor rho.  This explains the
# smaller reduced denominator seen after the A4 specialization.
special_inverse_numerator = sp.Poly(
    sp.expand(inverse_numerator_reduced.subs(jly_substitution)),
    T,
)
for coefficient in special_inverse_numerator.all_coeffs():
    quotient, remainder = sp.div(coefficient, -3 * rho, a, b)
    assert remainder == 0
    assert sp.Poly(quotient, a, b).is_zero is False

special_inverse_denominator = sp.factor(
    inverse_denominator_reduced.subs(jly_substitution) / (-3 * rho)
)
assert sp.factor(
    special_inverse_denominator - 432 * rho * B0**4 * sigma**2
) == 0

# Target-only polynomiality already fails for the primitive split.  The
# T^3 coefficient of primitive mod P is
#
#   -4608 * B0^5 * rho * residual,
#
# so division by orientation leaves residual/sigma.  The gcd check proves
# that this is a genuine pole, not a presentation artifact.
residual = (
    a**4
    - a**3 * b
    - 6 * a**3
    - a * b**3
    + 27 * a
    + b**4
    + 6 * b**3
    + 27 * b**2
    + 54 * b
    + 81
)
assert sp.factor(A0**2 - B0 * C0 - 3 * rho * residual) == 0
assert sp.gcd(
    sp.Poly(sigma, a, b),
    sp.Poly(B0 * residual, a, b),
).as_expr() == 1

special_primitive_mod_P = sp.Poly(
    sp.expand(primitive_mod_P.subs(jly_substitution)),
    T,
)
assert sp.factor(
    special_primitive_mod_P.nth(3)
    + 4608 * B0**5 * rho * residual
) == 0


# ---------------------------------------------------------------------------
# 4. Compare the derivative unit with the ordinary U,V incidence chart
# ---------------------------------------------------------------------------

U, V = sp.symbols("U V")

H = (
    8 * U**3
    - 6 * U * V**2
    - 18 * U * V
    - 54 * U
    - 2 * V**3
    - 9 * V**2
    - 27 * V
    - 27
)
K = 4 * U**2 + 4 * U * V + 6 * U + V**2 + 3 * V + 9
M = U**2 + 2 * V**2 + 6 * V + 18
L = (
    U**3
    - 3 * U * V**2
    - 9 * U * V
    - 27 * U
    + 2 * V**3
    + 9 * V**2
    + 27 * V
    + 27
)
N1 = sp.expand(M * K)
N2 = (
    8 * U**3 * V
    + 12 * U**2 * V**2
    + 36 * U**2 * V
    + 108 * U**2
    + 6 * U * V**3
    + 36 * U * V**2
    + 108 * U * V
    + 162 * U
    + V**4
    + 9 * V**3
    + 27 * V**2
    + 54 * V
)
source_A = U**3 - V**3 - 9 * V**2 - 27 * V - 54
source_rho = V**2 + 3 * V + 9
source_sigma = (
    2 * U**3 * V
    + 3 * U**3
    - 3 * U**2 * V**2
    - 9 * U**2 * V
    - 27 * U**2
    + V**4
    + 6 * V**3
    + 27 * V**2
    + 54 * V
    + 81
)

# Write A0(alpha,beta)=K^3*A6/H^3 and similarly for C0.  The target B0
# pullback is the already known K^3*L^2/H^3.
A6 = sp.cancel(
    (
        N1**3
        - N2**3
        - 9 * N2**2 * H
        - 27 * N2 * H**2
        - 54 * H**3
    )
    / K**3
)
C6 = sp.cancel((N1**3 - N2**3 + 27 * H**3) / K**3)
assert sp.denom(A6) == 1
assert sp.denom(C6) == 1

B0_pullback_numerator = (
    N1**3
    - 3 * N1 * N2**2
    + 2 * N2**3
    - 9 * N1 * N2 * H
    + 9 * N2**2 * H
    - 27 * N1 * H**2
    + 27 * N2 * H**2
    + 27 * H**3
)
assert sp.expand(B0_pullback_numerator - K**3 * L**2) == 0

# The selected scaled root has a compact rational expression.
root_numerator = 3 * source_A * K**3 * L
root_identity = (
    81 * source_A**4
    - 54 * A6 * source_A**2
    - 24 * source_A * L**3
    + 9 * A6**2
    - 12 * C6 * L**2
)
assert sp.expand(root_identity) == 0

# Theta is polynomial and is the new divisor appearing in P'(T) on this
# root branch.
Theta = sp.cancel(
    (27 * source_A**3 - 9 * A6 * source_A - 2 * L**3) / (2 * L)
)
assert sp.denom(Theta) == 1

root_derivative_numerator = (
    4 * root_numerator**3
    - 12 * (K**3 * A6) * (K**3 * L**2) * root_numerator
    - 8 * (K**3 * L**2) ** 3
)
assert sp.expand(
    root_derivative_numerator - 8 * K**9 * L**4 * Theta
) == 0
root_derivative = sp.factor(8 * K**9 * L**4 * Theta / H**9)
target_map = sp.Matrix([N1 / H, N2 / H])
target_map_jacobian = sp.factor(
    target_map.jacobian((U, V)).det()
)
assert sp.factor(target_map_jacobian - 4 * K**3 * L / H**3) == 0

# Pull back the three factors of the square-discriminant orientation.
rho_pullback_numerator = N2**2 + 3 * N2 * H + 9 * H**2
assert sp.expand(
    rho_pullback_numerator - source_rho * K**3
) == 0

sigma_pullback_numerator = (
    2 * N1**3 * N2
    + 3 * N1**3 * H
    - 3 * N1**2 * N2**2
    - 9 * N1**2 * N2 * H
    - 27 * N1**2 * H**2
    + N2**4
    + 6 * N2**3 * H
    + 27 * N2**2 * H**2
    + 54 * N2 * H**3
    + 81 * H**4
)
assert sp.expand(
    sigma_pullback_numerator
    - K**3 * source_sigma * Theta
) == 0

orientation_pullback = sp.factor(
    1728
    * source_rho
    * source_sigma
    * Theta
    * K**18
    * L**8
    / H**18
)
assert sp.factor(
    orientation_pullback / root_derivative
    - 216 * source_rho * source_sigma * K**9 * L**4 / H**9
) == 0

# The ordinary threefold suspension (alpha,beta,Q/P') has this residual
# Jacobian.  Hence 1/P' does not cancel the U,V chart Jacobian.
ordinary_suspension_jacobian = sp.factor(
    target_map_jacobian / root_derivative
)
assert sp.factor(
    ordinary_suspension_jacobian
    - H**6 / (2 * Theta * K**6 * L**3)
) == 0

print("PASS: the scaled generic quartic has the compact rank-four inverse basis")
print("PASS: a quartic primitive splits 1/P' across two orientation coordinates")
print("PASS: the primitive retains the generic four-sheet root field")
print("PASS: the JLY specialization has square discriminant Omega^2")
print("PASS: the first split coordinate retains a genuine sigma pole")
print("PASS: the selected root and orientation have exact U,V pullback ledgers")
print("PASS: 1/P' leaves H^6/(2*Theta*K^6*L^3) on the ordinary chart")
print("NOTE: this is a localized incidence identity, not a polynomial Keller map")

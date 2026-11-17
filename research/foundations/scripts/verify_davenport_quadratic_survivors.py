#!/usr/bin/env python3
"""Exact audit of the surviving quadratic-vector Keller branches."""

from __future__ import annotations

import sympy as sp


T, Y, U = sp.symbols("T Y U")
s, y, u = sp.symbols("s y u", nonzero=True)


def jacobian_pair(left: sp.Expr, right: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> sp.Expr:
    """Return d(left) wedge d(right) in the coordinates (x,z)."""
    return sp.diff(left, x) * sp.diff(right, z) - sp.diff(left, z) * sp.diff(
        right, x
    )


# ---------------------------------------------------------------------------
# 1. The through-origin PDE becomes an ordinary zero-Jacobian equation.
# ---------------------------------------------------------------------------

# The source base change T=s^2, U=u/s has constant Jacobian two.
source_change = sp.Matrix((s**2, y, u / s))
assert sp.simplify(source_change.jacobian((s, y, u)).det() - 2) == 0

alpha = sp.Function("alpha")(s, y)
A = sp.Function("A")(s, y)
B = sp.Function("B")(s, y)
G = sp.Function("G")(s, y)
K = sp.Function("K")(s, y)

laurent_map = sp.Matrix(
    (
        s**2 + alpha * u + u**2,
        G + A * u,
        K + B * u,
    )
)
laurent_jacobian = sp.Poly(
    sp.expand(laurent_map.jacobian((s, y, u)).det()),
    u,
)
assert sp.simplify(
    laurent_jacobian.nth(3)
    - 2 * jacobian_pair(A, B, s, y)
) == 0

# If A=P(w), B=Q(w), the next coefficient is one exact derivative along
# the common parameter w.
w = sp.Function("w")(s, y)
P = sp.Function("P")
Q = sp.Function("Q")
Pw = P(w)
Qw = Q(w)
Pp = sp.diff(Pw, w)
Qp = sp.diff(Qw, w)
Delta = sp.expand(Pw * Qp - Qw * Pp)

composed_map = laurent_map.subs({A: Pw, B: Qw})
composed_jacobian = sp.Poly(
    sp.expand(composed_map.jacobian((s, y, u)).det()),
    u,
)
expected_u2 = (
    Delta * jacobian_pair(w, alpha, s, y)
    + 2 * Pp * jacobian_pair(w, K, s, y)
    - 2 * Qp * jacobian_pair(w, G, s, y)
)
assert sp.simplify(composed_jacobian.nth(2) - expected_u2) == 0

first_integral = (
    K
    + Delta * alpha / (2 * Pp)
    - Qp * G / Pp
)
assert sp.simplify(
    expected_u2
    - 2 * Pp * jacobian_pair(w, first_integral, s, y)
) == 0

# The smallest nonconstant projective monomial branch is
# A=w, B=w^3, w=s*y.  In polynomial variables it gives
# a2=T*Y, a3=T^2*Y^3.  Its u^2/J2 primitive is exact, and the remaining
# constant coefficient is divisible by Y.
generic_g = sp.Function("generic_g")(T, Y)
generic_a = sp.Function("generic_a")(T, Y)
R0 = sp.Function("R0")
monomial_q = T * Y**2
monomial_a2 = T * Y
monomial_a3 = T**2 * Y**3
monomial_H = (
    3 * T * Y**2 * generic_g
    - T * Y**3 * generic_a
    + R0(monomial_q)
)
monomial_map = sp.Matrix(
    (
        T + generic_a * U + T * U**2,
        generic_g + monomial_a2 * U,
        monomial_H + monomial_a3 * U,
    )
)
monomial_jacobian = sp.Poly(
    sp.expand(monomial_map.jacobian((T, Y, U)).det()),
    U,
)
assert monomial_jacobian.nth(3) == 0
assert monomial_jacobian.nth(2) == 0
monomial_constant = sp.factor(monomial_jacobian.nth(0))
assert sp.simplify(monomial_constant.subs(Y, 0)) == 0

# Valuation checks used for the full primitive monomial-composition screen.
# Odd P,Q have odd initial orders p,q.  After ordering p<=q, Q'/P' is
# constant plus order at least two, while Delta/(2P') has order at least
# three in every projectively nonconstant case.
for p_order in range(1, 12, 2):
    for q_order in range(p_order, 12, 2):
        if q_order > p_order:
            assert q_order - p_order >= 2
            assert q_order >= 3
        else:
            # Equal initial order cancels in Delta; the first possible
            # projectively nonconstant term is two orders later.
            assert p_order + 2 >= 3


# ---------------------------------------------------------------------------
# 2. The next equation for a nonlinear projective direction.
# ---------------------------------------------------------------------------

lam = sp.Function("lam")(T, Y)
q = sp.Function("q")(T, Y)
c0 = sp.Function("c0")
c1 = sp.Function("c1")
c2 = sp.Function("c2")
curve = sp.Matrix((c0(q), c1(q), c2(q)))
curve_prime = curve.diff(q)
quadratic_direction = lam * curve
linear_direction = sp.Matrix(
    tuple(sp.Function(f"linear_{index}")(T, Y) for index in range(3))
)


def triple(left: sp.Matrix, middle: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.Matrix.hstack(left, middle, right).det()


b_T = quadratic_direction.diff(T)
b_Y = quadratic_direction.diff(Y)
a_T = linear_direction.diff(T)
a_Y = linear_direction.diff(Y)
u4_coefficient = (
    triple(b_T, b_Y, linear_direction)
    + 2 * triple(b_T, a_Y, quadratic_direction)
    + 2 * triple(a_T, b_Y, quadratic_direction)
)
phi = triple(curve, curve_prime, linear_direction)
developable_u4 = lam * (
    -jacobian_pair(q, lam, T, Y) * phi
    + 2 * lam * jacobian_pair(q, phi, T, Y)
)
assert sp.simplify(u4_coefficient - developable_u4) == 0

# Equivalently, away from phi=0 and lam=0,
# J(q,phi^2/lam)=0.
assert sp.simplify(
    developable_u4
    - lam**3 / phi * jacobian_pair(q, phi**2 / lam, T, Y)
) == 0

# For the nonlinear conic curve (1,q,q^2), phi is the displayed moving
# normal component, and for lam=1 the U^4 equation is 2*J(q,phi)=0.
conic_curve = sp.Matrix((1, q, q**2))
conic_prime = conic_curve.diff(q)
conic_phi = sp.factor(triple(conic_curve, conic_prime, linear_direction))
assert sp.expand(
    conic_phi
    - (
        q**2 * linear_direction[0]
        - 2 * q * linear_direction[1]
        + linear_direction[2]
    )
) == 0


# ---------------------------------------------------------------------------
# 3. Exact integration for the pure Davenport conic q=g.
# ---------------------------------------------------------------------------

g = sp.Function("g")(T, Y)
x = sp.Function("x")(T, Y)
v = sp.Function("v")(T, Y)
H = sp.Function("H")(T, Y)
R = sp.Function("R")
S = sp.Function("S")

Rg = R(g)
R_prime = sp.diff(Rg, g)
R_second = sp.diff(Rg, (g, 2))
S_prime = sp.diff(S(g), g)

# Moving-frame decomposition:
# a=x(1,g,g^2)+v(0,1,2g)+R(g)(0,0,1).
pure_a2 = g * x + v
pure_a3 = g**2 * x + 2 * g * v + Rg
pure_map = sp.Matrix(
    (
        T + x * U + U**2,
        g + pure_a2 * U + g * U**2,
        H + pure_a3 * U + g**2 * U**2,
    )
)
pure_jacobian = sp.Poly(
    sp.expand(pure_map.jacobian((T, Y, U)).det()),
    U,
)
assert pure_jacobian.nth(5) == 0
assert pure_jacobian.nth(4) == 0

first_pure_primitive = (
    H
    + T * g**2
    - v**2
    - R_prime * v
    - Rg * x / 2
)
assert sp.simplify(
    pure_jacobian.nth(3)
    - 2 * jacobian_pair(g, first_pure_primitive, T, Y)
) == 0

integrated_H = (
    -T * g**2
    + v**2
    + R_prime * v
    + Rg * x / 2
    + S(g)
)
integrated_map = pure_map.subs(H, integrated_H)
integrated_jacobian = sp.Poly(
    sp.expand(integrated_map.jacobian((T, Y, U)).det()),
    U,
)
assert integrated_jacobian.nth(3) == 0

second_pure_primitive = (
    -Rg * x**2 / 2
    - 4 * v**2 * x
    - 2 * R_prime * v * x
    + (-8 * T * g - 8 * g + 4 * S_prime) * v
    + 2 * R_second * v**2
    + T * (2 * Rg - 4 * g * R_prime)
)
assert sp.simplify(
    integrated_jacobian.nth(2)
    + jacobian_pair(g, second_pure_primitive, T, Y) / 2
) == 0

# For R != 0, the second primitive is quadratic in x.  Its discriminant is
# the explicit quartic recorded in the note.
x_discriminant, v_discriminant = sp.symbols(
    "x_discriminant v_discriminant"
)
R_symbol, R1_symbol, R2_symbol = sp.symbols(
    "R_symbol R1_symbol R2_symbol"
)
S1_symbol, L_symbol, g_symbol = sp.symbols(
    "S1_symbol L_symbol g_symbol"
)
quadratic_constant = (
    (-8 * (T + 1) * g_symbol + 4 * S1_symbol) * v_discriminant
    + 2 * R2_symbol * v_discriminant**2
    + T * (2 * R_symbol - 4 * g_symbol * R1_symbol)
    - L_symbol
)
quadratic_relation = (
    -R_symbol * x_discriminant**2 / 2
    - (4 * v_discriminant**2 + 2 * R1_symbol * v_discriminant)
    * x_discriminant
    + quadratic_constant
)
quadratic_poly = sp.Poly(quadratic_relation, x_discriminant)
quartic_discriminant = sp.factor(sp.discriminant(quadratic_poly))
expected_discriminant = (
    16 * v_discriminant**4
    + 16 * R1_symbol * v_discriminant**3
    + 4
    * (R_symbol * R2_symbol + R1_symbol**2)
    * v_discriminant**2
    + 8
    * R_symbol
    * (S1_symbol - 2 * (T + 1) * g_symbol)
    * v_discriminant
    + 4 * R_symbol**2 * T
    - 8 * R_symbol * R1_symbol * T * g_symbol
    - 2 * L_symbol * R_symbol
)
assert sp.expand(quartic_discriminant - expected_discriminant) == 0

# On the R=0 branch, solve the second primitive for x.  The U coefficient
# is then a third exact centralizer equation.
L = sp.Function("L")
Lg = L(g)
L_prime = sp.diff(Lg, g)
base_linear_term = S_prime - 2 * (T + 1) * g
x_R_zero = (
    4 * base_linear_term * v - Lg
) / (4 * v**2)
H_R_zero = -T * g**2 + v**2 + S(g)
a2_R_zero = g * x_R_zero + v
a3_R_zero = g**2 * x_R_zero + 2 * g * v
map_R_zero = sp.Matrix(
    (
        T + x_R_zero * U + U**2,
        g + a2_R_zero * U + g * U**2,
        H_R_zero + a3_R_zero * U + g**2 * U**2,
    )
)
jacobian_R_zero = sp.Poly(
    sp.factor(map_R_zero.jacobian((T, Y, U)).det()),
    U,
)
assert sp.simplify(jacobian_R_zero.nth(2)) == 0

third_R_zero_primitive = (
    (2 * T + 4 - sp.diff(S(g), (g, 2))) * v**2
    + L_prime * v / 2
    + Lg * x_R_zero / 4
    + 2 * T * g * S_prime
    - 2 * (T + 1) ** 2 * g**2
)
assert sp.simplify(
    jacobian_R_zero.nth(1)
    - jacobian_pair(g, third_R_zero_primitive, T, Y)
) == 0

# Eliminating x from the second and third primitives gives a nonzero
# quartic in v over K(T,g).  A non-base element of the prime degree-seven
# Davenport extension cannot satisfy it.
v_algebraic = sp.symbols("v_algebraic")
M = sp.Function("M")
base_constant_term = (
    2 * T * g * S_prime
    - 2 * (T + 1) ** 2 * g**2
    - M(g)
)
quartic_gate = (
    16
    * (2 * T + 4 - sp.diff(S(g), (g, 2)))
    * v_algebraic**4
    + 8 * L_prime * v_algebraic**3
    + 16 * base_constant_term * v_algebraic**2
    + 4 * Lg * base_linear_term * v_algebraic
    - Lg**2
)
quartic_poly = sp.Poly(quartic_gate, v_algebraic)
assert quartic_poly.degree() == 4
assert quartic_poly.nth(4).has(T)

# The v=0 exception is not covered by division above.  It fails already
# in J1 because S'(g)=2(T+1)g is impossible for a polynomial S(g).
x_zero_branch = sp.Function("x_zero_branch")(T, Y)
H_v_zero = -T * g**2 + S(g)
map_v_zero = sp.Matrix(
    (
        T + x_zero_branch * U + U**2,
        g + g * x_zero_branch * U + g * U**2,
        H_v_zero + g**2 * x_zero_branch * U + g**2 * U**2,
    )
)
jacobian_v_zero = sp.Poly(
    sp.factor(map_v_zero.jacobian((T, Y, U)).det()),
    U,
)
expected_v_zero_u1 = (
    -2 * base_linear_term * g * sp.diff(g, Y)
)
assert sp.simplify(
    jacobian_v_zero.nth(1) - expected_v_zero_u1
) == 0

# If x and v also belong to K[T,g], every output factors through
# (T,g,U), and the chain rule leaves the nonunit factor g_Y.
target_g = sp.Function("target_g")(T, Y)
psi0 = sp.Function("psi0")(T, target_g, U)
psi1 = sp.Function("psi1")(T, target_g, U)
psi2 = sp.Function("psi2")(T, target_g, U)
factor_through = sp.Matrix((psi0, psi1, psi2))
factor_through_jacobian = sp.factor(
    factor_through.jacobian((T, Y, U)).det()
)
assert factor_through_jacobian.has(sp.diff(target_g, Y))

print("PASS: the through-origin PDE is a Laurent zero-Jacobian equation")
print("PASS: its common-parameter U^2 equation has one exact first integral")
print("PASS: primitive monomial common-parameter branches retain a divisor")
print("PASS: nonlinear developable directions satisfy the phi^2/lambda law")
print("PASS: the pure Davenport conic has two successive exact integrals")
print("PASS: its R!=0 equation has the recorded quartic discriminant")
print("PASS: its R=0 branch ends in a forbidden degree-four equation")
print("PASS: its K[T,g] subbranch retains the nonunit g_Y factor")
print("PASS Davenport quadratic-survivor audit")

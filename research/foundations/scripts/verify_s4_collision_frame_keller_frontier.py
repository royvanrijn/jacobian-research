#!/usr/bin/env python3
"""Exact checks for the S4 edge-action collision-frame Keller frontier.

The script verifies the displayed polynomial identities, the decomposable
proper-nonabelian Keller example, the two-mask rational/logarithmic lifts,
and the small affine-source obstruction models.  Group-theoretic and class
group arguments which are uniform in characteristic zero are written out in
the companion note; the finite-field counts here are independent motivic
regressions, not substitutes for those proofs.
"""

from __future__ import annotations

from itertools import combinations, permutations, product

import sympy as sp


def assert_zero(expression: sp.Expr) -> None:
    assert sp.factor(expression) == 0


# ---------------------------------------------------------------------------
# 1. A literal decomposable proper-nonabelian Keller map.
# ---------------------------------------------------------------------------

x, y, z = sp.symbols("x y z")
u = 1 + x * y
foundational = sp.Matrix(
    [
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    ]
)
assert sp.factor(foundational.jacobian((x, y, z)).det()) == -2

composite = sp.Matrix(
    [
        sp.expand(entry.xreplace(dict(zip((x, y, z), foundational))))
        for entry in foundational
    ]
)
assert sp.factor(composite.jacobian((x, y, z)).det()) == 4
assert [sp.Poly(entry, x, y, z).total_degree() for entry in composite] == [43, 37, 25]

# The generic inverse cubic of the foundational map has group S3.  The
# specialization T^3-2T^2-1 is irreducible and has nonsquare discriminant
# -59, while the generic discriminant is visibly nonsquare.
T, A, B, C = sp.symbols("T A B C")
inverse_cubic = C * T**3 - 2 * T**2 + B * T - 2 * A
inverse_cubic_discriminant = sp.factor(sp.discriminant(inverse_cubic, T))
assert_zero(inverse_cubic_discriminant + 4 * (
    27 * A**2 * C**2 - 18 * A * B * C + 16 * A + B**3 * C - B**2
))
cubic_witness = sp.Poly(T**3 - 2 * T**2 - 1, T, domain=sp.QQ)
assert sp.factor_list(cubic_witness.as_expr(), T)[1] == [(cubic_witness.as_expr(), 1)]
assert sp.discriminant(cubic_witness.as_expr(), T) == -59

# The composition tower has three blocks of size three.  Its block quotient
# is S3 and its monodromy lies in the proper wreath product S3 wr S3.
assert 6**4 < sp.factorial(9)


# ---------------------------------------------------------------------------
# 2. Quartic factorization and the six-edge action.
# ---------------------------------------------------------------------------

p, a, m = sp.symbols("p a m")
quadratic_a = T**2 + a * T + (p + a**2 - m) / 2
quadratic_b = T**2 - a * T + (p + a**2 + m) / 2
quartic_product = sp.Poly(sp.expand(quadratic_a * quadratic_b), T)
q = a * m
r = ((p + a**2) ** 2 - m**2) / 4
assert quartic_product == sp.Poly(T**4 + p * T**2 + q * T + r, T)

J = m**2 + 2 * a**2 * (p + a**2)
phi = sp.Matrix([p, q, r])
phi_jacobian = phi.jacobian((p, a, m))
assert_zero(phi_jacobian.det() + J / 2)
assert_zero(sp.resultant(quadratic_a, quadratic_b, T) - J)

primitive = T**6 + 2 * p * T**4 + (p**2 - 4 * r) * T**2 - q**2
assert_zero(primitive.subs(T, a))
assert_zero(sp.diff(primitive, T).subs(T, a) - 2 * a * J)

quartic = T**4 + p * T**2 + q * T + r
quartic_discriminant = sp.factor(sp.discriminant(quartic, T))
d_minus = a**2 + 2 * p - 2 * m
d_plus = a**2 + 2 * p + 2 * m
assert_zero(quartic_discriminant - d_minus * d_plus * J**2)
assert_zero(sp.discriminant(quadratic_a, T) + d_minus)
assert_zero(sp.discriminant(quadratic_b, T) + d_plus)

# One good specialization forces the generic depressed quartic group to S4:
# X^4-X-1 is irreducible modulo 2 (a 4-cycle) and has type (1,3) modulo 7.
quartic_witness = T**4 - T - 1


def factor_degrees_mod(poly: sp.Expr, prime: int) -> list[int]:
    _, factors = sp.factor_list(poly, T, modulus=prime)
    return sorted(
        degree
        for factor, exponent in factors
        for degree in [sp.degree(factor, T)]
        for _ in range(exponent)
    )


assert factor_degrees_mod(quartic_witness, 2) == [4]
assert factor_degrees_mod(quartic_witness, 7) == [1, 3]
assert sp.discriminant(quartic_witness, T) == -283

# Enumerate the faithful S4 action on the six two-subsets.  The stabilizer of
# one edge has order four, and complementary edges form three blocks.
vertices = tuple(range(4))
edges = tuple(combinations(vertices, 2))
edge_index = {edge: index for index, edge in enumerate(edges)}


def edge_action(permutation: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        edge_index[tuple(sorted((permutation[i], permutation[j])))]
        for i, j in edges
    )


edge_actions = {edge_action(permutation) for permutation in permutations(vertices)}
assert len(edge_actions) == 24
assert {action[0] for action in edge_actions} == set(range(6))
assert sum(action[0] == 0 for action in edge_actions) == 4

complement = {
    edge_index[edge]: edge_index[tuple(sorted(set(vertices) - set(edge)))]
    for edge in edges
}
edge_blocks = {
    frozenset((edge_number, complement[edge_number])) for edge_number in range(6)
}
assert len(edge_blocks) == 3
for action in edge_actions:
    assert {
        frozenset((action[next(iter(block))], action[complement[next(iter(block))]]))
        for block in edge_blocks
    } == edge_blocks

# The complementary-edge involution sends (a,m) to (-a,-m).  Its basic
# invariants u=a^2, q=am, v=m^2 satisfy the singular quadric relation q^2=uv.
edge_u, edge_q, edge_v = sp.symbols("edge_u edge_q edge_v")
edge_quotient_relation = edge_q**2 - edge_u * edge_v
assert_zero(q**2 - a**2 * m**2)
assert all(
    sp.diff(edge_quotient_relation, variable).subs(
        {edge_u: 0, edge_q: 0, edge_v: 0}
    ) == 0
    for variable in (edge_u, edge_q, edge_v)
)


# ---------------------------------------------------------------------------
# 3. Tiny normal form and conductor ledger.
# ---------------------------------------------------------------------------

zz, xx, yy = sp.symbols("z x y")
p_tiny = (zz - 2 * xx**2 - yy**2) / 2
q_tiny = xx**2 * yy
c_tiny = xx**2 * (zz - xx**2 - 2 * yy**2)
theta = sp.Matrix([p_tiny, q_tiny, c_tiny])
assert_zero(theta.jacobian((zz, xx, yy)).det() + xx**3 * zz)
assert_zero(J.subs({p: p_tiny, a: xx, m: xx * yy}) - xx**2 * zz)

source_rechart = sp.Matrix([p_tiny, xx, xx * yy])
assert_zero(source_rechart.jacobian((zz, xx, yy)).det() - xx / 2)
assert_zero((T**6 + 2 * p_tiny * T**4 - c_tiny * T**2 - q_tiny**2).subs(T, xx))

P, Q, R, S = sp.symbols("P Q R S")
delta_target = sp.factor(sp.discriminant(T**4 + P * T**2 + Q * T + R, T))
delta_c_target = sp.factor(delta_target.subs(R, (S + P**2) / 4))
sextic_target = T**6 + 2 * P * T**4 - S * T**2 - Q**2
assert_zero(sp.discriminant(sextic_target, T) - 64 * Q**2 * delta_c_target**2)

# The all-degree one-normal zero-section rank obstruction has this most
# general derivative pattern at a=m=w=0.
c1, c2, c3, d1, d2, d3, d4 = sp.symbols("c1 c2 c3 d1 d2 d3 d4")
one_normal_derivative = sp.Matrix(
    [
        [1, 0, 0, c1],
        [0, 0, 0, c2],
        [p / 2, 0, 0, c3],
        [d1, d2, d3, d4],
    ]
)
assert one_normal_derivative.det() == 0

# If all three coefficient outputs are retained in a two-mask extension, the
# full determinant retains det(D Phi), independently of the lower-left block.
l11, l12, l13, l21, l22, l23 = sp.symbols("l11 l12 l13 l21 l22 l23")
n11, n12, n21, n22 = sp.symbols("n11 n12 n21 n22")
block_derivative = phi_jacobian.row_join(sp.zeros(3, 2)).col_join(
    sp.Matrix([[l11, l12, l13, n11, n12], [l21, l22, l23, n21, n22]])
)
assert_zero(block_derivative.det() - phi_jacobian.det() * (n11 * n22 - n12 * n21))

tau, sigma, normal_n = sp.symbols("tau sigma normal_n")
normal_substitution = {
    m: -2 * a * tau + sigma,
    p: -a**2 - 2 * tau**2 + normal_n,
}
assert_zero(
    J.subs(normal_substitution)
    - (sigma**2 - 4 * a * tau * sigma + 2 * a**2 * normal_n)
)


# ---------------------------------------------------------------------------
# 4. Rational cotangent and polynomial logarithmic lifts.
# ---------------------------------------------------------------------------

z1, z2 = sp.symbols("z1 z2")
A0 = sp.Matrix([[m, a], [a * (p + a**2), -m / 2]])
assert A0 == sp.Matrix([q, r]).jacobian((a, m))
mask_rational = sp.Matrix(
    [
        (m * z1 + 2 * a * (p + a**2) * z2) / J,
        (2 * a * z1 - 2 * m * z2) / J,
    ]
)
# With the coefficient block fixed, neither adjugate numerator is divisible
# by J for independent masks.  Any polynomialization must therefore change
# the coefficient block or impose a non-affine mask relation.
for numerator in (
    m * z1 + 2 * a * (p + a**2) * z2,
    2 * a * z1 - 2 * m * z2,
):
    assert sp.Poly(numerator, m).rem(sp.Poly(J, m)).as_expr() != 0
assert sp.simplify(mask_rational - A0.inv().T * sp.Matrix([z1, z2])) == sp.zeros(2, 1)
rational_lift = sp.Matrix([p, q, r, *mask_rational])
assert sp.factor(rational_lift.jacobian((p, a, m, z1, z2)).det()) == 1

target_matrix = sp.Matrix(
    [
        [2 * P * (P**2 - 4 * R) + 9 * Q**2, Q * (P**2 + 12 * R)],
        [
            Q * (P**2 + 12 * R),
            16 * R**2 - 4 * P**2 * R + sp.Rational(3, 2) * P * Q**2,
        ],
    ]
)
assert_zero(target_matrix.det() + P * delta_target / 2)

# Its columns are relative logarithmic derivations in (Q,R).  It is a Saito
# basis only after P is inverted, since det(T) has the extra factor P.
delta_q = sp.diff(delta_target, Q)
delta_r = sp.diff(delta_target, R)
assert_zero(target_matrix[0, 0] * delta_q + target_matrix[1, 0] * delta_r - 36 * Q * delta_target)
assert_zero(
    target_matrix[0, 1] * delta_q
    + target_matrix[1, 1] * delta_r
    + 4 * (P**2 - 12 * R) * delta_target
)

source_matrix = sp.Matrix(
    [
        [a**2 + 2 * p, 2 * a * m],
        [2 * a * m, -2 * (a**2 * p - 2 * m**2 + 2 * p**2)],
    ]
)
target_pullback = target_matrix.subs({P: p, Q: q, R: r})
assert sp.simplify(target_pullback - A0 * source_matrix * A0.T) == sp.zeros(2, 2)

# Multiplying the rational masks by the target matrix gives the polynomial
# mask A0*H*z.  Its determinant is the pullback of -P*Delta/2.
polynomial_mask_matrix = A0 * source_matrix
assert all(sp.denom(sp.cancel(entry)) == 1 for entry in polynomial_mask_matrix)
assert_zero(
    phi_jacobian.det() * polynomial_mask_matrix.det()
    - target_matrix.det().subs({P: p, Q: q, R: r})
)


# ---------------------------------------------------------------------------
# 5. Affine-source recognition and failed equal-quadratic slices.
# ---------------------------------------------------------------------------

# In difference coordinates, two general monic quadratics are A and A+hT+e.
# Their resultant is rho=e^2-aeh+bh^2.
aa, bb, ee, hh = sp.symbols("aa bb ee hh")
base_quadratic = T**2 + aa * T + bb
shifted_quadratic = base_quadratic + hh * T + ee
rho = ee**2 - aa * ee * hh + bb * hh**2
assert_zero(sp.resultant(base_quadratic, shifted_quadratic, T) - rho)

for level_equation in (rho - 1, rho - hh):
    smoothness_ideal = [level_equation] + [
        sp.diff(level_equation, variable) for variable in (aa, bb, ee, hh)
    ]
    assert sp.groebner(smoothness_ideal, aa, bb, ee, hh).polys == [sp.Poly(1, aa, bb, ee, hh)]

# The asymmetric level has e^2=h(1+ae-bh); along E=(e,h), h has divisor
# 2E.  Localizing h eliminates b and gives a Laurent polynomial ring.
assert_zero((rho - hh) - (ee**2 - hh * (1 + aa * ee - bb * hh)))
localized_b = (hh + aa * ee * hh - ee**2) / hh**2
assert_zero((rho - hh).subs(bb, localized_b))

# The symmetric unimodular completion is exactly SL2 after b=2a.  Four root
# locally nilpotent derivations preserve the determinant relation; the first
# two kernels already contain all four coordinate generators in pairs.
mm, bc, cc, vv = sp.symbols("mm bc cc vv")
sl2_relation = mm * vv - bc * cc - 1


def derivation(poly: sp.Expr, values: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    return sp.expand(sum(sp.diff(poly, variable) * value for variable, value in values.items()))


left_upper = {mm: cc, bc: vv, cc: 0, vv: 0}
left_lower = {mm: 0, bc: 0, cc: mm, vv: bc}
right_upper = {mm: 0, bc: mm, cc: 0, vv: cc}
right_lower = {mm: bc, bc: 0, cc: vv, vv: 0}
for lnd in (left_upper, left_lower, right_upper, right_lower):
    assert derivation(sl2_relation, lnd) == 0
assert all(derivation(variable, left_upper) == 0 for variable in (cc, vv))
assert all(derivation(variable, left_lower) == 0 for variable in (mm, bc))

# Point counts realize the three Grothendieck/Hodge ledgers over every tested
# odd finite field: [SL2]=L^3-L, [rho=1]=L^3+L^2, [rho=h]=L^3.
for prime in (3, 5, 7):
    values = range(prime)
    sl2_count = sum(
        (m0 * v0 - b0 * c0 - 1) % prime == 0
        for m0, b0, c0, v0 in product(values, repeat=4)
    )
    rho_one_count = sum(
        (e0 * e0 - a0 * e0 * h0 + b0 * h0 * h0 - 1) % prime == 0
        for a0, b0, e0, h0 in product(values, repeat=4)
    )
    rho_h_count = sum(
        (e0 * e0 - a0 * e0 * h0 + b0 * h0 * h0 - h0) % prime == 0
        for a0, b0, e0, h0 in product(values, repeat=4)
    )
    assert sl2_count == prime**3 - prime
    assert rho_one_count == prime**3 + prime**2
    assert rho_h_count == prime**3

# For UA+VB=T^3, put W=U+V=T-tau and use values/slopes of the two
# linears V,D=B-A at T=tau.  Polynomial division leaves xy-tau^3, giving
# the exact A2 surface singularity xy=tau^3 times an affine plane.
alpha, beta, x_value, y_value = sp.symbols("alpha beta x_value y_value")
W_linear = T - tau
V_linear = alpha * (T - tau) + x_value
D_linear = beta * (T - tau) + y_value
quotient_a, remainder_a = sp.div(T**3 - V_linear * D_linear, W_linear, T)
assert_zero(remainder_a - (tau**3 - x_value * y_value))
assert sp.LC(sp.Poly(quotient_a, T)) == 1
U_linear = W_linear - V_linear
B_quadratic = quotient_a + D_linear
bezout_identity = sp.expand(U_linear * quotient_a + V_linear * B_quadratic - T**3)
assert_zero(bezout_identity - (x_value * y_value - tau**3))

# For UA+VB=T^2 the coefficient model has a positive-dimensional singular
# locus vh=1, u=e=b=0.  The two gradients have rank one there.
v_slope, u_value, h_slope, e_value, coeff_a, coeff_b = sp.symbols(
    "v_slope u_value h_slope e_value coeff_a coeff_b"
)
w_constant = 1 - v_slope * h_slope
f1 = w_constant * coeff_a + v_slope * e_value + u_value * h_slope
f2 = w_constant * coeff_b + u_value * e_value
variables = (coeff_a, coeff_b, v_slope, u_value, h_slope, e_value)
singular_jacobian = sp.Matrix([f1, f2]).jacobian(variables).subs(
    {h_slope: 1 / v_slope, u_value: 0, e_value: 0, coeff_b: 0}
)
assert singular_jacobian.rank() == 1
assert_zero(f1.subs({h_slope: 1 / v_slope, u_value: 0, e_value: 0, coeff_b: 0}))
assert_zero(f2.subs({h_slope: 1 / v_slope, u_value: 0, e_value: 0, coeff_b: 0}))

print("PASS: F composed with itself is a degree-nine decomposable Keller map with determinant four")
print("PASS: the quartic collision frame has determinant -J/2 and exact S4 six-edge monodromy")
print("PASS: discriminant, primitive-conductor, tiny-normal-form, and two-normal ledgers are exact")
print("PASS: the rational cotangent lift has determinant one and the relative Saito factorization is polynomial")
print("PASS: SL2, resultant-level, asymmetric-divisor, and repeated-root Bezout models agree")

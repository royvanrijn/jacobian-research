#!/usr/bin/env python3
"""Exact degree-six Ritt atlas on the normalized admissible seed slice.

The computation is deliberately made only after passing to the Hessian
divisor: H is obtained by integrating a generic quartic K=H'' twice and the
two endpoint conditions are imposed before either decomposition is
eliminated.  All equations are over Q and are understood on the exact-degree
open h6 != 0.
"""

from __future__ import annotations

import sympy as sp


w, z = sp.symbols("w z")
h2, h3, h4, h5, h6 = sp.symbols("h2 h3 h4 h5 h6")
k0, k1, k2, k3, k4 = sp.symbols("k0 k1 k2 k3 k4")

# Hessian-first normalized sextic slice.
K = k0 + k1 * w + k2 * w**2 + k3 * w**3 + k4 * w**4
H_from_K = sp.integrate(sp.integrate(K, w), w)
assert H_from_K.subs(w, 0) == sp.diff(H_from_K, w).subs(w, 0) == 0

H = h2 * w**2 + h3 * w**3 + h4 * w**4 + h5 * w**5 + h6 * w**6
hessian_substitution = {
    k0: 2 * h2,
    k1: 6 * h3,
    k2: 12 * h4,
    k3: 20 * h5,
    k4: 30 * h6,
}
assert sp.expand(H_from_K.subs(hessian_substitution) - H) == 0

endpoint_zero = sp.expand(H.subs(w, 1))
endpoint_derivative = sp.expand(sp.diff(H, w).subs(w, 1) + 1)
assert endpoint_zero == h2 + h3 + h4 + h5 + h6
assert endpoint_derivative == 2 * h2 + 3 * h3 + 4 * h4 + 5 * h5 + 6 * h6 + 1

# Eliminate h2 only after both Hessian-integral endpoint constraints are in
# place.  The remaining affine threefold chart has one linear equation L=0.
h2_normal = -(h3 + h4 + h5 + h6)
normalization = h3 + 2 * h4 + 3 * h5 + 4 * h6 + 1
assert sp.expand(endpoint_derivative.subs(h2, h2_normal) - normalization) == 0

# -------------------------------------------------------------------------
# The 2 o 3 atlas.  Normalize the inner cubic to B(0)=0 and leading term w^3.
# The constant and linear terms are invisible to H'', so coefficient matching
# below is exactly matching the already constrained Hessian divisor.
a, b, p, q = sp.symbols("a b p q")
B23 = w**3 + p * w**2 + q * w
f23 = sp.expand(a * B23**2 + b * B23)
H23 = sp.expand(f23 - sp.diff(f23, w).subs(w, 0) * w)
coeff23 = {
    h6: a,
    h5: 2 * a * p,
    h4: a * (p**2 + 2 * q),
    h3: 2 * a * p * q + b,
    h2: a * q**2 + b * p,
}
assert sp.expand(sp.diff(H23, w, 2) - sp.diff(H, w, 2).subs(coeff23)) == 0

phi23 = sp.expand(
    32 * h3 * h5 * h6**2
    + 64 * h3 * h6**3
    + 16 * h4**2 * h6**2
    - 24 * h4 * h5**2 * h6
    + 64 * h4 * h6**3
    + 5 * h5**4
    + 64 * h5 * h6**3
    + 64 * h6**4
)
assert sp.factor(
    phi23.subs(coeff23) - 64 * a**3 * H23.subs(w, 1)
) == 0

# Rational reconstruction proves sufficiency on D(h6), not just necessity.
p_back = h5 / (2 * h6)
q_back = sp.cancel((h4 / h6 - p_back**2) / 2)
b_back = sp.cancel(h3 - 2 * h6 * p_back * q_back)
reconstructed_h2_23 = sp.cancel(h6 * q_back**2 + b_back * p_back)
num23 = sp.factor(sp.together(h2 - reconstructed_h2_23).as_numer_denom()[0])
assert sp.expand(num23.subs(h2, h2_normal) + phi23) == 0

# A dense normalized parameter chart.
D23 = 2 * p**2 + 5 * p - q + 3
a23 = sp.cancel(-(p + 1) / ((p + q + 1) * D23))
b23 = sp.cancel((p + q + 1) / D23)
H23_normal = sp.factor(H23.subs({a: a23, b: b23}))
assert sp.factor(H23_normal.subs(w, 1)) == 0
assert sp.factor(sp.diff(H23_normal, w).subs(w, 1) + 1) == 0

# -------------------------------------------------------------------------
# The 3 o 2 atlas, with monic B(0)=0.
c = sp.symbols("c")
B32 = w**2 + p * w
f32 = sp.expand(a * B32**3 + b * B32**2 + c * B32)
H32 = sp.expand(f32 - sp.diff(f32, w).subs(w, 0) * w)
coeff32 = {
    h6: a,
    h5: 3 * a * p,
    h4: 3 * a * p**2 + b,
    h3: a * p**3 + 2 * b * p,
    h2: b * p**2 + c,
}
assert sp.expand(sp.diff(H32, w, 2) - sp.diff(H, w, 2).subs(coeff32)) == 0

phi32 = 27 * h3 * h6**2 - 18 * h4 * h5 * h6 + 5 * h5**3
assert sp.expand(phi32.subs(coeff32)) == 0

p_back_32 = h5 / (3 * h6)
b_back_32 = sp.cancel(h4 - 3 * h6 * p_back_32**2)
c_back_32 = sp.cancel(h2 - b_back_32 * p_back_32**2)
num32 = sp.factor(
    sp.together(
        h3 - (h6 * p_back_32**3 + 2 * b_back_32 * p_back_32)
    ).as_numer_denom()[0]
)
assert sp.expand(num32 - phi32) == 0

b32 = sp.cancel(
    -(a * p**3 + 6 * a * p**2 + 9 * a * p + 4 * a + 1)
    / (2 * (p + 1))
)
c32 = sp.cancel(
    a * p**4 / 2
    + 5 * a * p**3 / 2
    + 9 * a * p**2 / 2
    + 7 * a * p / 2
    + a
    + p / 2
    + sp.Rational(1, 2)
)
H32_normal = sp.factor(H32.subs({b: b32, c: c32}))
assert sp.factor(H32_normal.subs(w, 1)) == 0
assert sp.factor(sp.diff(H32_normal, w).subs(w, 1) + 1) == 0

# -------------------------------------------------------------------------
# Ritt intersection.  The coefficient ideal is (normalization,phi23,phi32).
# On the dense 2 o 3 chart it is a single explicit plane curve.
ritt_curve = sp.expand(
    4 * p**4
    + 4 * p**3
    - 18 * p**2 * q
    - 27 * p**2
    - 72 * p * q
    - 54 * p
    - 27 * q**2
    - 54 * q
    - 27
)
phi32_on_23 = sp.factor(phi32.subs(coeff23).subs({a: a23, b: b23}))
expected_phi32_on_23 = sp.factor(
    -(p + 1) ** 2
    * ritt_curve
    / ((p + q + 1) ** 3 * D23**3)
)
assert sp.factor(phi32_on_23 - expected_phi32_on_23) == 0

# A rational point on the Ritt curve which is Hessian-clean and boundary-clean.
intersection_witness = sp.factor(H23_normal.subs({p: -3, q: -2}))
assert intersection_witness == -w**2 * (w - 1) * (w**3 - 5 * w**2 + 20) / 16
assert sp.factor(sp.discriminant(sp.diff(intersection_witness, w, 2), w)) != 0
assert sp.gcd(
    sp.Poly(intersection_witness / w**2, w),
    sp.Poly(sp.diff(intersection_witness, w) / w, w),
).degree() == 0

# -------------------------------------------------------------------------
# Omitted-value components.
# Every quadratic outer polynomial becomes a square after adding a constant:
# hence closure(O_222)=D_23.  Exact type (2,2,2) is the cubic-squarefree open.
r23 = sp.cancel(b / (2 * a))
s23 = -b * q
t23 = sp.cancel(b**2 / (4 * a))
E222 = sp.expand(H23 - s23 * w + t23)
assert sp.factor(E222 - a * (B23 + r23) ** 2) == 0

# A cubic outer polynomial becomes a cube exactly when 3ac=b^2.  In seed
# coordinates this is psi33=0 inside D_32.
psi33 = 9 * h2 * h6**2 - 3 * h4**2 * h6 + h4 * h5**2
assert sp.factor(psi33.subs(coeff32) - 3 * a * (3 * a * c - b**2)) == 0
r32 = sp.cancel(b / (3 * a))
s32 = -c * p
t32 = sp.cancel(b**3 / (27 * a**2))
E33 = sp.expand(H32 - s32 * w + t32)
assert sp.factor(
    E33.subs(c, b**2 / (3 * a)) - a * (B32 + r32) ** 3
) == 0

# Direct normalization of O_33 supplies a useful rational exact-(3,3) witness.
u, v, lam = sp.symbols("u v lam")
Q33 = w**2 + u * w + v
H33 = sp.expand(lam * (Q33**3 - v**3 - 3 * u * v**2 * w))
tangent33 = sp.expand((1 + u + v) ** 3 - v**3 - 3 * u * v**2)
derivative33 = sp.expand(3 * (1 + u + v) ** 2 * (2 + u) - 3 * u * v**2)
assert sp.expand(H33.subs(w, 1) - lam * tangent33) == 0
assert sp.expand(sp.diff(H33, w).subs(w, 1) - lam * derivative33) == 0
o33_witness = sp.factor(
    H33.subs({u: sp.Rational(1, 2), v: -sp.Rational(3, 2), lam: sp.Rational(8, 27)})
)
assert o33_witness == w**2 * (w - 1) * (8 * w**3 + 20 * w**2 - 10 * w - 45) / 27
assert sp.discriminant(Q33.subs({u: sp.Rational(1, 2), v: -sp.Rational(3, 2)}), w) != 0

# O_33 intersect D_23.  Uniqueness of the omitted value forces a sixth
# power.  The tangent condition gives four reduced algebraic seeds.
rho = sp.symbols("rho")
type6_tangent = 15 * rho**4 - 20 * rho**3 + 15 * rho**2 - 6 * rho + 1
type6_derivative = 6 * (5 * rho**4 - 10 * rho**3 + 10 * rho**2 - 5 * rho + 1)
H6 = sp.expand(
    -((w - rho) ** 6 - rho**6 + 6 * rho**5 * w) / type6_derivative
)
assert sp.factor(H6.subs(w, 1) + type6_tangent / type6_derivative) == 0
assert sp.factor(sp.diff(H6, w).subs(w, 1) + 1) == 0
assert sp.gcd(type6_tangent, type6_derivative) == 1
assert sp.discriminant(type6_tangent, rho) == 10800
type6_coefficients = {
    hj: sp.Poly(H6, w).coeff_monomial(w**j)
    for hj, j in ((h2, 2), (h3, 3), (h4, 4), (h5, 5), (h6, 6))
}
for equation in (normalization, phi23, phi32, psi33):
    numerator = sp.together(equation.subs(type6_coefficients)).as_numer_denom()[0]
    assert sp.rem(sp.Poly(numerator, rho), sp.Poly(type6_tangent, rho)) == 0
type6_admissibility = sp.factor(sp.diff(H6, w, 2).subs(w, 1) + 2)
assert sp.gcd(
    type6_tangent,
    sp.together(type6_admissibility).as_numer_denom()[0],
) == 1
assert sp.discriminant(sp.diff(H6, w, 2), w) == 0

# The same four points appear from the dense Ritt chart plus psi33.
psi33_on_23 = sp.factor(psi33.subs(coeff23).subs({a: a23, b: b23}))
psi33_numerator = sp.factor(sp.together(psi33_on_23).as_numer_denom()[0])
intersection_resultant = sp.factor(sp.resultant(ritt_curve, psi33_numerator, q))
type6_in_p = 5 * p**4 + 20 * p**3 + 45 * p**2 + 54 * p + 27
assert intersection_resultant == 9 * (p + 1) ** 6 * type6_in_p**2
assert sp.factor(type6_in_p - 27 * type6_tangent.subs(rho, -p / 3)) == 0

# -------------------------------------------------------------------------
# Degeneration of the affine root-sheet distinction.  On exact degree this
# is exactly failure of the exact-double, boundary-clean condition:
# P=H/w^2 and 2P+wP' acquire a common root.  Root zero accounts for h2=0;
# nonzero common roots are additional multiple primitive roots.
P = sp.cancel(H / w**2)
affine_sheet_boundary = sp.factor(
    sp.resultant(P, 2 * P + w * sp.diff(P, w), w).subs(h2, h2_normal)
)
basic_boundary_factors = (
    h6
    * (h3 + h4 + h5 + h6)
    * (h3 + 2 * h4 + 3 * h5 + 4 * h6) ** 2
)
boundary_quotient = sp.cancel(affine_sheet_boundary / basic_boundary_factors)
assert boundary_quotient.as_numer_denom()[1] == 1
assert sp.Poly(boundary_quotient, h3, h4, h5, h6).total_degree() == 4

# Compact factorizations after pulling this boundary back to the two dense
# Ritt charts.  Denominators and exact-degree units are intentionally removed.
delta23_zero_cluster = sp.expand(
    p**3 + 2 * p**2 * q + 2 * p**2 + 2 * p * q + p - q**2
)
delta23_extra_root = sp.expand(
    4 * p**2 * q + p**2 - 10 * p * q - 2 * p - 27 * q**2 - 14 * q - 3
)
boundary23 = sp.factor(
    affine_sheet_boundary.subs(coeff23).subs({a: a23, b: b23})
)
boundary23_expected = sp.factor(
    -(p + 1) ** 3
    * delta23_zero_cluster
    * delta23_extra_root
    / ((p + q + 1) ** 4 * D23**6)
)
assert sp.factor(boundary23 - boundary23_expected) == 0

delta32_extra_1 = a * p**3 + 4 * a * p**2 + 5 * a * p + 2 * a + 1
delta32_zero_cluster = 5 * a * p**3 + 12 * a * p**2 + 9 * a * p + 2 * a + 2 * p + 1
delta32_extra_2 = sp.expand(
    25 * a**2 * p**6
    + 155 * a**2 * p**5
    + 379 * a**2 * p**4
    + 457 * a**2 * p**3
    + 272 * a**2 * p**2
    + 64 * a**2 * p
    + 23 * a * p**3
    + 47 * a * p**2
    + 24 * a * p
    - 2
)
boundary32 = sp.factor(
    affine_sheet_boundary.subs(coeff32).subs({b: b32, c: c32})
)
boundary32_expected = sp.factor(
    -a**2
    * delta32_extra_1
    * delta32_zero_cluster
    * delta32_extra_2
    / (8 * (p + 1) ** 4)
)
assert sp.factor(boundary32 - boundary32_expected) == 0

# The common Ritt curve is not trapped in the boundary: the rational witness
# above is the promised concrete clean point.
assert sp.factor(ritt_curve.subs({p: -3, q: -2})) == 0
assert sp.factor(delta23_zero_cluster.subs({p: -3, q: -2})) != 0
assert sp.factor(delta23_extra_root.subs({p: -3, q: -2})) != 0

# Boundary cuts on the common Ritt curve.  These are the exact projected
# equations on the dense p,q chart; the coefficient ideal retains p=-1.
common_zero_projection = sp.factor(
    sp.resultant(ritt_curve, delta23_zero_cluster, q)
)
common_extra_projection = sp.factor(
    sp.resultant(ritt_curve, delta23_extra_root, q)
)
assert common_zero_projection == -(
    (p + 1) ** 2
    * (
        560 * p**6
        + 2520 * p**5
        + 3699 * p**4
        + 540 * p**3
        - 3402 * p**2
        - 2916 * p
        - 729
    )
)
assert sp.expand(
    common_extra_projection
    - 108
    * (p - 3)
    * (p + 1) ** 3
    * (20 * p**4 + 100 * p**3 - 49 * p**2 - 624 * p - 576)
) == 0

# Four small exact witnesses separate the open pieces of the atlas.
d23_only = sp.factor(H23_normal.subs({p: 1, q: 1}))
assert sp.expand(
    d23_only - w**2 * (-2 * w**4 - 4 * w**3 - 6 * w**2 + 5 * w + 7) / 27
) == 0
d23_only_coefficients = {
    h2: sp.Poly(d23_only, w).coeff_monomial(w**2),
    h3: sp.Poly(d23_only, w).coeff_monomial(w**3),
    h4: sp.Poly(d23_only, w).coeff_monomial(w**4),
    h5: sp.Poly(d23_only, w).coeff_monomial(w**5),
    h6: sp.Poly(d23_only, w).coeff_monomial(w**6),
}
assert sp.factor(phi32.subs(d23_only_coefficients)) != 0

d32_only = sp.factor(H32_normal.subs({a: 1, p: 0}))
assert sp.expand(
    d32_only - w**2 * (w - 1) * (w + 1) * (2 * w**2 - 3) / 2
) == 0
d32_only_coefficients = {
    h2: sp.Poly(d32_only, w).coeff_monomial(w**2),
    h3: sp.Poly(d32_only, w).coeff_monomial(w**3),
    h4: sp.Poly(d32_only, w).coeff_monomial(w**4),
    h5: sp.Poly(d32_only, w).coeff_monomial(w**5),
    h6: sp.Poly(d32_only, w).coeff_monomial(w**6),
}
assert sp.factor(phi23.subs(d32_only_coefficients)) != 0

for witness in (d23_only, d32_only, intersection_witness, o33_witness):
    assert witness.subs(w, 1) == 0
    assert sp.diff(witness, w).subs(w, 1) == -1
    assert sp.diff(witness, w, 2).subs(w, 1) != -2
    assert sp.discriminant(sp.diff(witness, w, 2), w) != 0
    assert sp.gcd(
        sp.Poly(witness / w**2, w),
        sp.Poly(sp.diff(witness, w) / w, w),
    ).degree() == 0

print("PASS degree-six Ritt atlas: Hessian-first normalized seed chart")
print("PASS degree-six Ritt atlas: exact saturated equations for 2o3 and 3o2")
print("PASS degree-six Ritt atlas: explicit common Ritt curve and clean rational witness")
print("PASS degree-six Ritt atlas: O_222=D_23 and O_33 is a curve in D_32")
print("PASS degree-six Ritt atlas: omitted intersection has four type-(6) support points, each doubled")
print("PASS degree-six Ritt atlas: affine-sheet degeneration cuts factor on both charts")

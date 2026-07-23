#!/usr/bin/env python3
"""Exact regression for the Hessian-incidence solution of degree-six OP-RITT."""

import sympy as sp


W, x = sp.symbols("W x")
k0, k1, k2, k3, k4 = sp.symbols("k0 k1 k2 k3 k4")
p, q, r = sp.symbols("p q r")


# The 3 o 2 condition: after translating a quartic to remove its cubic term,
# its linear term must also vanish.
center = -k3 / (4 * k4)
K = k4 * W**4 + k3 * W**3 + k2 * W**2 + k1 * W + k0
phi32 = k3**3 - 4 * k4 * k3 * k2 + 8 * k4**2 * k1
assert sp.factor(8 * k4**2 * sp.diff(K, W).subs(W, center)) == phi32


# The 2 o 3 condition.  Normalize the leading coefficient of the cubic and
# of its Hessian, then solve successively for p,q,r.  The last coefficient is
# the quartic implicit equation displayed in the note.
C = W**3 + p * W**2 + q * W + r
KC = sp.Poly(sp.diff(C**2, W, 2) / 30, W)
assert KC.nth(4) == 1
assert KC.nth(3) == 4 * p / 3
assert KC.nth(2) == 2 * p**2 / 5 + 4 * q / 5
assert KC.nth(1) == 2 * p * q / 5 + 2 * r / 5
assert KC.nth(0) == q**2 / 15 + 2 * p * r / 15

p_sol = 3 * k3 / 4
q_sol = sp.expand((5 * k2 - 2 * p_sol**2) / 4)
r_sol = sp.expand((5 * k1 - 2 * p_sol * q_sol) / 2)
last_relation = sp.factor(
    15 * k0 - q_sol**2 - 2 * p_sol * r_sol
)
phi23_affine = (
    3072 * k0
    - 768 * k1 * k3
    - 320 * k2**2
    + 432 * k2 * k3**2
    - 81 * k3**4
)
assert sp.factor(last_relation - 5 * phi23_affine / 1024) == 0

phi23 = (
    3072 * k0 * k4**3
    - 768 * k1 * k3 * k4**2
    - 320 * k2**2 * k4**2
    + 432 * k2 * k3**2 * k4
    - 81 * k3**4
)
assert sp.Poly(phi23, k0, k1, k2, k3, k4).total_degree() == 4


# On the 3 o 2 locus, translate to the center.  The 2 o 3 equation becomes
# 48*k0*k4-5*k2^2, and its solutions are exactly Hessians of
# (x^3+q*x)^2.
centered_phi23 = sp.factor(phi23.subs({k3: 0, k1: 0}))
assert centered_phi23 == 64 * k4**2 * (48 * k0 * k4 - 5 * k2**2)

common = (x**3 + q * x) ** 2
common_hessian = sp.expand(sp.diff(common, x, 2))
assert common_hessian == 30 * x**4 + 24 * q * x**2 + 2 * q**2
assert sp.expand(common - x**2 * (x**2 + q) ** 2) == 0


# Pull the two projective equations back to the normalized sextic seed chart.
h3, h4, h5, h6 = sp.symbols("h3 h4 h5 h6")
seed_substitution = {
    k0: -2 * (h3 + h4 + h5 + h6),
    k1: 6 * h3,
    k2: 12 * h4,
    k3: 20 * h5,
    k4: 30 * h6,
}
seed_phi32 = 27 * h3 * h6**2 - 18 * h4 * h5 * h6 + 5 * h5**3
seed_phi23 = (
    32 * h3 * h5 * h6**2
    + 64 * h3 * h6**3
    + 16 * h4**2 * h6**2
    - 24 * h4 * h5**2 * h6
    + 64 * h4 * h6**3
    + 5 * h5**4
    + 64 * h5 * h6**3
    + 64 * h6**4
)
assert sp.factor(phi32.subs(seed_substitution) - 1600 * seed_phi32) == 0
assert sp.factor(phi23.subs(seed_substitution) + 2592000 * seed_phi23) == 0

# The two hypersurfaces meet in the expected codimension two at the clean
# rational witness below (h3,h4,h5,h6)=(4,-25/2,15,-25/4).
normalized_h3 = -1 - 2 * h4 - 3 * h5 - 4 * h6
normalized_equations = [
    seed_phi32.subs(h3, normalized_h3),
    seed_phi23.subs(h3, normalized_h3),
]
witness_coefficients = {
    h4: -sp.Rational(25, 2),
    h5: 15,
    h6: -sp.Rational(25, 4),
}
intersection_jacobian = sp.Matrix(
    [
        [sp.diff(equation, variable).subs(witness_coefficients)
         for variable in (h4, h5, h6)]
        for equation in normalized_equations
    ]
)
assert intersection_jacobian.rank() == 2


# The marked normalized intersection curve.
c = sp.symbols("c")
g = ((W - c) ** 3 + q * (W - c)) ** 2
g0 = g.subs(W, 0)
gp0 = sp.diff(g, W).subs(W, 0)
endpoint = sp.factor(g.subs(W, 1) - g0 - gp0)
expected_endpoint = (
    15 * c**4
    - 20 * c**3
    + 12 * c**2 * q
    + 15 * c**2
    - 8 * c * q
    - 6 * c
    + q**2
    + 2 * q
    + 1
)
assert sp.expand(endpoint - expected_endpoint) == 0
assert sp.factor(sp.discriminant(endpoint, q)) == 4 * c * (3 * c - 1) * (
    7 * c**2 - 7 * c + 2
)

D = sp.factor(sp.diff(g, W).subs(W, 1) - gp0)
Q = sp.cancel((g - g0 - gp0 * W) / W**2)
assert sp.factor(Q.subs(W, 0)) == 15 * c**4 + 12 * c**2 * q + q**2
assert sp.factor(sp.discriminant(Q, W)) == (
    16
    * c**2
    * (c**2 + q) ** 2
    * (3 * c**2 + q)
    * (225 * c**4 + 285 * c**2 * q + 64 * q**2)
)
assert sp.factor(sp.discriminant(sp.diff(g, W, 2), W)) == 108380160 * q**6


# A rational point on the clean curve, including both decompositions and all
# normalization conditions.
witness = {c: sp.Rational(2, 5), q: -sp.Rational(1, 5)}
assert endpoint.subs(witness) == 0
assert D.subs(witness) == sp.Rational(4, 25)

lambda_witness = -1 / D.subs(witness)
H = sp.factor(
    lambda_witness * (g - g0 - gp0 * W).subs(witness)
)
s0 = sp.factor((gp0 / D).subs(witness))
expected_H = (
    -sp.Rational(25, 4) * W**6
    + 15 * W**5
    - sp.Rational(25, 2) * W**4
    + 4 * W**3
    - sp.Rational(1, 4) * W**2
)
assert sp.expand(H - expected_H) == 0
assert s0 == sp.Rational(7, 125)
assert H.subs(W, 0) == 0
assert sp.diff(H, W).subs(W, 0) == 0
assert H.subs(W, 1) == 0
assert sp.diff(H, W).subs(W, 1) == -1
assert sp.diff(H, W, 2).subs(W, 1) == -14
assert sp.diff(H, W, 2).subs(W, 0) != 0
assert sp.discriminant(sp.diff(H, W, 2), W) != 0
assert sp.discriminant(Q.subs(witness), W) != 0

C_witness = (W - sp.Rational(2, 5)) ** 3 - (W - sp.Rational(2, 5)) / 5
f = sp.expand(H - s0 * W)
outer_quadratic = -sp.Rational(25, 4) * (
    C_witness**2 - sp.Rational(4, 15625)
)
Z = sp.symbols("Z")
outer_cubic = -sp.Rational(25, 4) * (
    Z * (Z - sp.Rational(1, 5)) ** 2 - sp.Rational(4, 15625)
)
inner_quadratic = (W - sp.Rational(2, 5)) ** 2
assert sp.expand(f - outer_quadratic) == 0
assert sp.expand(f - outer_cubic.subs(Z, inner_quadratic)) == 0

print("PASS: vertical decomposition is detected entirely by the projective Hessian")
print("PASS: the degree-six 3o2 and 2o3 Hessian hypersurfaces are exact")
print("PASS: their intersection is precisely the monomial/Chebyshev Ritt collision")
print("PASS: the normalized marked intersection curve has a rational clean witness")

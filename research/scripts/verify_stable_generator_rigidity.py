#!/usr/bin/env python3
"""Exact regressions for low poles, affine pencil transport, and collisions."""

from __future__ import annotations

import sympy as sp


w, r = sp.symbols("w r")
s_h, t_h, s_g, t_g = sp.symbols("s_H t_H s_G t_G")
lam, mu, p = sp.symbols("lambda mu p", nonzero=True)
a = sp.cancel(lam / mu)
d = sp.expand(lam * p)


# A symbolic affine identity.  The geometric target map goes H -> G, so its
# contravariant pullback is
# beta^*(s_G,t_G)=(lambda*s_H,mu*t_H+lambda*p*s_H).
h0, h1, h2, h3, h4 = sp.symbols("h0:5")
H_symbolic = h0 + h1 * w + h2 * w**2 + h3 * w**3 + h4 * w**4
G_symbolic = sp.expand(mu * H_symbolic.subs(w, a * (w - p)))
E_g = G_symbolic - s_g * w + t_g
beta_pullback_E_g = sp.expand(
    E_g.subs({s_g: lam * s_h, t_g: mu * t_h + d * s_h})
)
scaled_E_h = sp.expand(
    mu * (H_symbolic.subs(w, a * (w - p)) - s_h * a * (w - p) + t_h)
)
assert sp.cancel(beta_pullback_E_g - scaled_E_h) == 0

# The inverse substitution is the convention used when old target coordinates
# are expressed in the new presentation (as in G=beta F alpha).
inverse_target_expression = sp.expand(
    mu
    * (
        H_symbolic.subs(w, a * (w - p))
        - (s_g / lam) * a * (w - p)
        + (t_g - d * s_g / lam) / mu
    )
)
assert sp.cancel(E_g - inverse_target_expression) == 0

# The same orientation on the discriminant normalization is
# r_G=(mu/lambda)r_H+p, not its inverse.
r_g = sp.cancel(mu * r / lam + p)
G_prime = sp.diff(G_symbolic, w)
assert sp.cancel(G_prime.subs(w, r_g) - lam * sp.diff(H_symbolic, w).subs(w, r)) == 0
G_t = sp.expand(r_g * G_prime.subs(w, r_g) - G_symbolic.subs(w, r_g))
H_t = sp.expand(r * sp.diff(H_symbolic, w).subs(w, r) - H_symbolic.subs(w, r))
assert sp.cancel(G_t - (mu * H_t + d * sp.diff(H_symbolic, w).subs(w, r))) == 0


# Exact normalized rerooting at the extra root 3.
H = sp.expand(sp.Rational(1, 2) * w**2 * (w - 1) * (w - 3))
root = sp.Integer(3)
kappa = sp.cancel(-1 / (root * sp.diff(H, w).subs(w, root)))
G = sp.expand(kappa * H.subs(w, root * w))
assert H.subs(w, 0) == sp.diff(H, w).subs(w, 0) == H.subs(w, 1) == 0
assert sp.diff(H, w).subs(w, 1) == -1
assert G.subs(w, 0) == sp.diff(G, w).subs(w, 0) == G.subs(w, 1) == 0
assert sp.diff(G, w).subs(w, 1) == -1

lambda_root = sp.expand(kappa * root)
mu_root = kappa
E_G = G - s_g * w + t_g
reroot_pullback = sp.expand(
    E_G.subs({s_g: lambda_root * s_h, t_g: mu_root * t_h})
)
E_H_rerooted = sp.expand(
    mu_root * (H.subs(w, root * w) - s_h * root * w + t_h)
)
assert sp.expand(reroot_pullback - E_H_rerooted) == 0
assert sp.cancel(mu_root / lambda_root) == sp.Rational(1, 3)

# The root-one G sheet maps to the extra H root 3.  It is simple, hence a
# polar extra-root component rather than H's distinguished affine root-one.
assert H.subs(w, root) == 0
assert sp.diff(H, w).subs(w, root) != 0
assert root != 1


# The universal discriminant is the point-in-LL-divisor incidence for the
# line P_s=H-sW: its critical value at r is -t.
P_s = sp.expand(H - s_h * w)
assert sp.expand(sp.diff(P_s, w).subs({w: r, s_h: sp.diff(H, w).subs(w, r)})) == 0
critical_value = sp.expand(P_s.subs({w: r, s_h: sp.diff(H, w).subs(w, r)}))
assert sp.expand(critical_value + (r * sp.diff(H, w).subs(w, r) - H.subs(w, r))) == 0


# Higher zero clusters.  For m=2 both roots lie in the W=kC chart.  For
# m>=3 that chart contains one affine root and the other m-1 roots have the
# dicritical valuation C=delta^(m-1), W=k*delta.
B, A, C, c, h, k, delta = sp.symbols("B A C c h k delta", nonzero=True)
for multiplicity in range(2, 7):
    local_zero = h * w**multiplicity - B * C * w + c * A * C**2
    finite_chart = sp.expand(local_zero.subs(w, k * C))
    coefficient_c2 = finite_chart.coeff(C, 2)
    if multiplicity == 2:
        assert sp.expand(coefficient_c2 - (h * k**2 - B * k + c * A)) == 0
    else:
        assert sp.expand(coefficient_c2 - (-B * k + c * A)) == 0
        boundary_chart = sp.expand(
            local_zero.subs({C: delta ** (multiplicity - 1), w: k * delta})
        )
        leading = boundary_chart.coeff(delta, multiplicity)
        assert sp.factor(leading - k * (h * k ** (multiplicity - 1) - B)) == 0


# Every nonzero multiple-root cluster is polar.  Its roots have
# C=delta^mu, W=rho+k*delta, so y=(W-gamma)/C starts with rho/delta^mu.
rho = sp.symbols("rho", nonzero=True)
for multiplicity in range(2, 7):
    local_multiple = h * (w - rho) ** multiplicity - B * C * w + c * A * C**2
    collision_chart = sp.expand(
        local_multiple.subs({C: delta**multiplicity, w: rho + k * delta})
    )
    assert sp.factor(
        collision_chart.coeff(delta, multiplicity) - (h * k**multiplicity - B * rho)
    ) == 0
    local_gamma = sp.expand(-sp.diff(local_multiple, w) / c)
    gamma_chart = sp.expand(
        local_gamma.subs({C: delta**multiplicity, w: rho + k * delta})
    )
    assert sp.expand((rho + k * delta) - gamma_chart).subs(delta, 0) == rho


# A simple extra root is still polar when its y-pole happens to cancel: then
# gamma and W both tend to rho, x tends to zero, and the z numerator tends to
# rho-1.  Only the distinguished root rho=1 can be regular.
hprime, a0 = sp.symbols("hprime a0", nonzero=True)
gamma_limit = -hprime / c
z_numerator_limit = sp.expand(
    gamma_limit - 1 - a0 * (rho / gamma_limit - 1)
)
cancelled_y_numerator = sp.factor((rho - gamma_limit).subs(hprime, -c * rho))
assert cancelled_y_numerator == 0
assert sp.factor(z_numerator_limit.subs(hprime, -c * rho) - (rho - 1)) == 0


# At a transverse collision u^mu=epsilon, normalization epsilon=delta^mu
# makes the derivative vanish to order mu-1 and cycles the mu reroot choices.
epsilon, u = sp.symbols("epsilon u")
for multiplicity in range(2, 7):
    collision_polynomial = u**multiplicity - epsilon
    normalized = collision_polynomial.subs({u: delta, epsilon: delta**multiplicity})
    assert sp.expand(normalized) == 0
    derivative = sp.diff(collision_polynomial, u).subs(u, delta)
    assert derivative == multiplicity * delta ** (multiplicity - 1)


# Low-pole filtration and implicit-degree checks in several degrees.  The
# arithmetic part is universal; resultants additionally protect the claimed
# monic-normal-form degree in t for exact sample pencils.
s, t = sp.symbols("s t")
for degree in range(4, 9):
    n = degree - 1
    sample = sp.expand(sum(sp.Integer(j + 1) * r**j for j in range(2, degree + 1)))
    sample_s = sp.diff(sample, r)
    sample_t = sp.expand(r * sample_s - sample)
    implicit = sp.resultant(s - sample_s, t - sample_t, r)
    implicit_poly = sp.Poly(implicit, t)
    assert implicit_poly.degree() == n

    weights = {}
    for j in range(n):
        for i in range(degree + 1):
            weight = i * n + j * degree
            assert weight not in weights
            weights[weight] = (i, j)

    below = {pair for weight, pair in weights.items() if weight <= degree - 2}
    through_s = {pair for weight, pair in weights.items() if weight <= degree - 1}
    through_t = {pair for weight, pair in weights.items() if weight <= degree}
    assert below == {(0, 0)}
    assert through_s == {(0, 0), (1, 0)}
    assert through_t == {(0, 0), (1, 0), (0, 1)}


print("PASS low-pole filtration: bases 1; 1,s; 1,s,t in degrees 4..8")
print("PASS LL incidence: D_H is the critical-value incidence of H-sW")
print("PASS affine target orientation: triangular pullback and inverse agree")
print("PASS affine pencil identity: beta^*E_G=mu*E_H(a*(W-p))")
print("PASS rerooting identity: normalized root 3 maps W_G=1 to W_H=3")
print("PASS collision models: higher zero clusters and nonzero multiple roots")
print("PASS simple-root reconstruction: only root one survives both y and z")

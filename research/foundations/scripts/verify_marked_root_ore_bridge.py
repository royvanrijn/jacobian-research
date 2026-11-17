#!/usr/bin/env python3
"""Exact bridge from marked-root coordinates to the Ore/Darboux chart.

The generic rank-two source chart is

    W = Z + psi(Q),  Y = Q - X*W/3,
    (x_seed,y_seed,z_seed) = (X,-3*Y/(2*a),-W/2).

On the affine marked-root chart put

    u = y_seed + 1/x_seed,  v = 1/x_seed.

The script proves that the Ore coordinate v is literally the marked-root
derivative coordinate, derives u in (X,Q,Z), and identifies the weighted
seed argument with R*u/2.  It also records the root-at-infinity valuation
relations.  No quantum boundary descent is asserted.
"""

import sympy as sp


a, X, Q, Z, psi = sp.symbols("a X Q Z psi", nonzero=True)
rho, u, v, P = sp.symbols("rho u v P")
s, d, marked = sp.symbols("s d marked", nonzero=True)

W = Z + psi
Y = Q - X * W / 3
seed_x = X
seed_y = -3 * Y / (2 * a)
seed_z = -W / 2

root_u = sp.factor(seed_y + 1 / seed_x)
root_v = 1 / seed_x
source_unit = sp.factor(1 + seed_x * seed_y)
gamma = sp.factor(1 - 3 * X * Q / 2)
R = sp.factor(2 * X - 3 * X**2 * Q)
seed_marked_argument = sp.factor(source_unit * gamma)

expected_root_u = v - 3 * Q / (2 * a) + (Z + psi) / (2 * a * v)
assert sp.factor(root_u.subs(X, 1 / v) - expected_root_u) == 0
assert sp.factor(source_unit - X * root_u) == 0
assert sp.factor(gamma - R / (2 * X)) == 0
assert sp.factor(seed_marked_argument - R * root_u / 2) == 0

# In exact localized Darboux coordinates P=-Z/3.  The displayed equation is
# the inverse bridge: it recovers P from the marked root and the central
# parameter after Q=v(2-rho*v)/3 is substituted.
q_vrho = v * (2 - rho * v) / 3
p_from_root = sp.factor(
    -(
        2 * a * v * u
        - 2 * a * v**2
        + 3 * v * q_vrho
        - psi
    )
    / 3
)
bridge_residual = expected_root_u.subs(
    {
        Q: q_vrho,
        Z: -3 * P,
    }
)
assert sp.factor(bridge_residual.subs(P, p_from_root) - u) == 0

# Root-at-infinity chart: s=1/u and d=v/u.  A finite weighted marked
# argument forces rho=2*s*marked.  Although v has a pole, gamma=rho*v/2 and
# the seed argument remain regular.
infinity_substitution = {
    u: 1 / s,
    v: d / s,
    rho: 2 * s * marked,
}
assert sp.factor((1 / v).subs(infinity_substitution) - s / d) == 0
assert sp.factor((u / v).subs(infinity_substitution) - 1 / d) == 0
assert sp.factor((rho * v / 2).subs(infinity_substitution) - d * marked) == 0
assert sp.factor((rho * u / 2).subs(infinity_substitution) - marked) == 0

# The original marked-root reconstruction has x=s/d and y=(1-d)/s.
assert sp.factor(
    (u - v).subs(infinity_substitution) - (1 - d) / s
) == 0

print("PASS: Ore v=X^-1 is exactly the affine marked-root derivative coordinate")
print("PASS: u=v-3Q/(2a)+(Z+psi(Q))/(2av) is the marked-root/Ore bridge")
print("PASS: the weighted seed argument is R*u/2")
print("PASS: at root infinity, R=2*s*marked and gamma=d*marked are regular")
print("SCOPE: exact classical coordinate and valuation identities only")

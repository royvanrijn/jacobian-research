#!/usr/bin/env python3
"""Exact certificate for the monodromy of the quartic inverse cover."""

import sympy as sp


W, h, epsilon = sp.symbols("W h epsilon")
s, t = sp.symbols("s t")

E = W**2 - W**4 - 2 * s * W + t

# The generic cover is connected.  Indeed, E is linear and monic in t, hence
# irreducible in Q[s, W, t]; this independent fraction-field check certifies
# the corresponding statement over Q(s,t)[W].
K = sp.QQ.frac_field(s, t)
assert sp.Poly(E, W, domain=K).degree() == 4
assert sp.Poly(E, W, domain=K).is_irreducible

# A transverse slice through the smooth discriminant point associated to r=1.
# Its unique double root splits with a square-root parameter, hence gives a
# transposition.  The simple zero of the slice discriminant checks transversality.
smooth = sp.factor(E.subs({s: -1, t: -2}))
assert smooth == -(W - 1) ** 2 * (W**2 + 2 * W + 2)
smooth_shift = sp.expand(E.subs({W: 1 + h, s: -1, t: -2 + epsilon}))
assert smooth_shift.subs({h: 0, epsilon: 0}) == 0
assert sp.diff(smooth_shift, h).subs({h: 0, epsilon: 0}) == 0
assert sp.diff(smooth_shift, h, 2).subs({h: 0, epsilon: 0}) != 0
assert sp.diff(smooth_shift, epsilon).subs({h: 0, epsilon: 0}) == 1
smooth_disc = sp.factor(sp.discriminant(E.subs({s: -1, t: -2 + epsilon}), W))
assert smooth_disc == -16 * epsilon * (16 * epsilon**2 - 88 * epsilon + 125)

# At either cusp the repeated root is triple, while the fourth root remains
# simple.  A t-slice has local form epsilon + c*h^3 + higher terms, so a loop
# around epsilon=0 cyclically permutes the three nearby roots: a 3-cycle.
for sign in (-1, 1):
    r = sign / sp.sqrt(6)
    s0 = sp.simplify(r - 2 * r**3)
    t0 = sp.simplify(r**2 - 3 * r**4)
    cusp = sp.factor(E.subs({s: s0, t: t0}), extension=sp.sqrt(6))
    assert sp.simplify(cusp + (W - r) ** 3 * (W + 3 * r)) == 0
    cusp_shift = sp.expand(E.subs({W: r + h, s: s0, t: t0 + epsilon}))
    origin = {h: 0, epsilon: 0}
    assert cusp_shift.subs(origin) == 0
    assert sp.diff(cusp_shift, h).subs(origin) == 0
    assert sp.diff(cusp_shift, h, 2).subs(origin) == 0
    assert sp.diff(cusp_shift, h, 3).subs(origin) != 0
    assert sp.diff(cusp_shift, epsilon).subs(origin) == 1

# At the node two independent double roots split simultaneously.  A loop in
# this slice exchanges both pairs, giving cycle type (2,2).
node_shift = sp.expand(E.subs({s: 0, t: -sp.Rational(1, 4) + epsilon}))
assert sp.expand(node_shift - (-(W**2 - sp.Rational(1, 2)) ** 2 + epsilon)) == 0
node_disc = sp.factor(sp.discriminant(node_shift, W))
assert node_disc == -64 * epsilon**2 * (4 * epsilon - 1)

# Irreducibility makes the geometric monodromy G <= S4 transitive, so 4 | |G|.
# The cusp cycle gives 3 | |G|, hence 12 | |G|.  Thus |G| is 12 or 24.  The
# smooth-point transposition is odd, excluding the unique index-two subgroup
# A4.  Therefore G=S4.  Since the model is over Q, arithmetic monodromy contains
# geometric monodromy and is consequently S4 as well.

print("quartic inverse irreducible over Q(s,t): connected four-sheet cover")
print("local branch cycles certified: transposition, 3-cycle, double transposition")
print("geometric and arithmetic monodromy: S4")

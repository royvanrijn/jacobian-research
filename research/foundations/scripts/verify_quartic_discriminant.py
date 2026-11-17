#!/usr/bin/env python3
"""Exact discriminant normalization and cusp/node audit for the quartic map."""

import sympy as sp


W, r, h, s, t = sp.symbols("W r h s t")
E = W**2 - W**4 - 2*s*W + t
Delta = 27*s**4 - 36*s**2*t - s**2 + 16*t**3 + 8*t**2 + t

assert sp.factor(sp.discriminant(E, W) + 16*Delta) == 0

# Repeated roots parameterize the discriminant.
s_norm = r - 2*r**3
t_norm = r**2 - 3*r**4
assert sp.factor(E.subs({W: r, s: s_norm, t: t_norm})) == 0
assert sp.factor(sp.diff(E, W).subs({W: r, s: s_norm})) == 0
assert sp.factor(Delta.subs({s: s_norm, t: t_norm})) == 0

elimination = sp.groebner(
    [s-s_norm, t-t_norm], r, s, t, order="lex"
)
target_only = [
    sp.factor(polynomial.as_expr())
    for polynomial in elimination.polys
    if not polynomial.as_expr().has(r)
]
assert target_only == [Delta]

# The parameter r is integral and generically recovered rationally, proving
# that A^1_r is the normalization of the irreducible discriminant curve.
assert sp.factor((r**2 - 3*s*r + 2*t).subs({s: s_norm, t: t_norm})) == 0
inverse_denominator = 18*s**2 - 4*t - 1
inverse_numerator = s*(12*t - 1)
assert sp.factor(
    (inverse_denominator*r - inverse_numerator).subs({s: s_norm, t: t_norm})
) == 0
assert sp.factor(Delta).is_Pow is False
assert sp.Poly(Delta, s, t).is_irreducible


# The complete singular scheme has one node and two cusps.
singular_gb = sp.groebner(
    [Delta, sp.diff(Delta, s), sp.diff(Delta, t)], s, t, order="lex"
)
assert [sp.factor(polynomial.as_expr()) for polynomial in singular_gb.polys] == [
    36*s**2 - 48*t**2 - 16*t - 1,
    s*(12*t - 1)**2,
    (4*t + 1)*(12*t - 1)**2,
]

# Node: r=+/-1/sqrt(2) map to (0,-1/4), with independent tangent vectors.
r_node_plus = 1/sp.sqrt(2)
r_node_minus = -r_node_plus
node = {s: 0, t: -sp.Rational(1, 4)}
for root in (r_node_plus, r_node_minus):
    assert sp.simplify(s_norm.subs(r, root)) == node[s]
    assert sp.simplify(t_norm.subs(r, root)) == node[t]
tangent = sp.Matrix([sp.diff(s_norm, r), sp.diff(t_norm, r)])
node_tangents = sp.Matrix.hstack(
    tangent.subs(r, r_node_plus), tangent.subs(r, r_node_minus)
)
assert sp.simplify(node_tangents.det()) != 0
assert sp.factor(E.subs(node) + (W**2-sp.Rational(1, 2))**2) == 0

# Cusps: r=+/-1/sqrt(6). The second and third parameter derivatives are
# independent, giving local form (h^2,h^3); the quartic has a triple root and
# one remaining simple root.
for root in (1/sp.sqrt(6), -1/sp.sqrt(6)):
    cusp_s = sp.simplify(s_norm.subs(r, root))
    cusp_t = sp.simplify(t_norm.subs(r, root))
    assert cusp_t == sp.Rational(1, 12)
    assert sp.simplify(cusp_s**2-sp.Rational(2, 27)) == 0
    second = sp.Matrix([
        sp.diff(s_norm, r, 2).subs(r, root),
        sp.diff(t_norm, r, 2).subs(r, root),
    ])
    third = sp.Matrix([
        sp.diff(s_norm, r, 3).subs(r, root),
        sp.diff(t_norm, r, 3).subs(r, root),
    ])
    assert sp.simplify(sp.Matrix.hstack(second, third).det()) != 0
    expected = -(W-root)**3*(W+3*root)
    assert sp.simplify(sp.expand(E.subs({s: cusp_s, t: cusp_t})-expected)) == 0

print("PASS: quartic discriminant is Delta=0 and A^1 normalizes it")
print("PASS: the normalization is integral and generically birational")
print("PASS: singularities are one double-double node and two triple-root cusps")
print("PASS: node tangents separate and cusp jets have orders two and three")

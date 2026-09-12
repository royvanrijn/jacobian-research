#!/usr/bin/env python3
"""Exact regressions for the simultaneous-multicluster LL comparison."""

from __future__ import annotations

import sympy as sp


w, r, s, t, u = sp.symbols("w r s t u")


def root_valuation(polynomial: sp.Expr, variable: sp.Symbol, root: sp.Expr) -> int:
    """Return the exact multiplicity of ``variable-root`` in a polynomial."""

    quotient = sp.Poly(sp.factor(polynomial), variable)
    divisor = sp.Poly(variable - root, variable)
    valuation = 0
    while quotient.eval(root) == 0:
        quotient, remainder = sp.div(quotient, divisor)
        assert remainder.is_zero
        valuation += 1
    return valuation


# Two distinct multiple primitive roots give branches with tangent lines
# t=rho*s.  The model branch with ramification e=m-1 is
# (s,t)=(u^e,rho*u^e+u^(e+1)).  Evaluating the other branch equation shows
# intersection multiplicity e_i*e_j whenever rho_i != rho_j.
rho_i, rho_j = sp.symbols("rho_i rho_j", nonzero=True)
for e_i in range(1, 5):
    for e_j in range(1, 5):
        branch_j = (t - rho_j * s) ** e_j - s ** (e_j + 1)
        on_branch_i = sp.expand(
            branch_j.subs({s: u**e_i, t: rho_i * u**e_i + u ** (e_i + 1)})
        )
        assert root_valuation(on_branch_i, u, 0) == e_i * e_j


# For normalized seeds with two and three simultaneous clusters, compute the
# global implicit discriminant and its adjunction conductor.  On branch i the
# predicted conductor exponent is e_i(E-1), where E=sum_j e_j.
cluster_profiles = (
    ((0, 2), (2, 2)),
    ((0, 2), (2, 3)),
    ((0, 3), (2, 3)),
    ((0, 2), (2, 4)),
    ((0, 2), (2, 2), (3, 2)),
    ((0, 3), (2, 2), (3, 2)),
)
for profile in cluster_profiles:
    raw_seed = w - 1
    for center, multiplicity in profile:
        raw_seed *= (w - center) ** multiplicity
    scale = sp.cancel(-1 / sp.diff(raw_seed, w).subs(w, 1))
    seed = sp.expand(scale * raw_seed)

    assert seed.subs(w, 0) == sp.diff(seed, w).subs(w, 0) == 0
    assert seed.subs(w, 1) == 0
    assert sp.diff(seed, w).subs(w, 1) == -1

    parameter_s = sp.diff(seed, w).subs(w, r)
    parameter_t = sp.expand(r * parameter_s - seed.subs(w, r))
    implicit = sp.factor(sp.resultant(s - parameter_s, t - parameter_t, r))
    assert sp.degree(implicit, t) == sp.degree(seed, w) - 1

    conductor = sp.cancel(
        sp.diff(implicit, t).subs({s: parameter_s, t: parameter_t})
        / sp.diff(seed, w, 2).subs(w, r)
    )
    conductor_numerator, conductor_denominator = sp.together(conductor).as_numer_denom()
    assert sp.Poly(conductor_denominator, r).degree() == 0

    total_e = sum(multiplicity - 1 for _, multiplicity in profile)
    for center, multiplicity in profile:
        ramification = multiplicity - 1
        assert root_valuation(conductor_numerator, r, center) == (
            ramification * (total_e - 1)
        )

    # The local parameter orders are e_i and e_i+1, and the tangent slopes
    # are the root centers 0 and 2.
    for center, multiplicity in profile:
        local_s = sp.expand(parameter_s.subs(r, center + u))
        local_q = sp.expand((parameter_t - center * parameter_s).subs(r, center + u))
        assert root_valuation(local_s, u, 0) == multiplicity - 1
        assert root_valuation(local_q, u, 0) == multiplicity


# The full marked-root incidence E=H-sW+t is smooth independently of the
# collision partition, because its t derivative is one.  Completing at a
# root W=rho and eliminating t therefore gives k[[s,u]].
generic_seed = sp.Function("H")
incidence = generic_seed(w) - s * w + t
assert sp.diff(incidence, t) == 1


print("PASS multicluster tangents: I(C_i,C_j)=e_i*e_j")
print("PASS multicluster conductor: c_i=e_i*(sum(e_j)-1)")
print("PASS marked-root incidence: every collision chart remains regular")

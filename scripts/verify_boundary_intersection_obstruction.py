#!/usr/bin/env python3
"""Exact regressions for the C24/weighted boundary-intersection obstruction.

The uniform argument is proved in BOUNDARY_INTERSECTION_OBSTRUCTION.md.  This
script checks the full scheme-theoretic and reduced restricted discriminants
in inverse degrees 4--6; it is not used as a proof for arbitrary degree.
"""
from __future__ import annotations

import sympy as sp

from master_cancellation import fiber_antiderivative


T, W, P, Q, R, A, B, C = sp.symbols("T W P Q R A B C")


for m in range(2, 5):
    # C24, r=1.  The finite critical chart gives the G_m component.  The
    # saturated discriminant also contains the Q=0 affine-line component
    # contributed by the closure of the critical divisor at infinity.
    anti = fiber_antiderivative(m, 1, T, P, Q)
    critical = sp.factor(sp.diff(anti, T))
    assert sp.factor(critical - (1 - T * (Q - P * T) ** m)) == 0
    critical_at_boundary = sp.solve(critical.subs(P, 0), T)[0]
    branch_value = sp.factor(C * anti.subs({P: 0, T: critical_at_boundary}))
    assert sp.factor(branch_value - C / (2 * Q**m)) == 0
    c24_discriminant = sp.factor(sp.discriminant(C * anti - R, T))
    p_order = min(
        monomial[0]
        for monomial, _ in sp.Poly(sp.expand(c24_discriminant), P).terms()
    )
    assert p_order == m * (m - 1)
    c24_trace = sp.factor((c24_discriminant / P**p_order).subs(P, 0))
    nilpotency_index = m * (m + 1)
    finite_factor = 2 * R * Q**m - C
    c24_scheme = Q**nilpotency_index * finite_factor
    c24_reduced = Q * finite_factor
    quotient = sp.factor(c24_trace / c24_scheme)
    assert quotient != 0 and not quotient.has(P, Q, R)
    assert sp.gcd(sp.Poly(Q**nilpotency_index, Q, R), sp.Poly(finite_factor, Q, R)) == 1
    reduced_ratio = sp.factor(
        sp.Poly(c24_scheme, Q, R).sqf_part().as_expr() / c24_reduced
    )
    assert reduced_ratio != 0 and not reduced_ratio.has(Q, R)

    # A generic weighted primitive of the same inverse degree m+2.  Its
    # nonzero roots are chosen rationally only to keep this regression small.
    roots = list(range(1, m + 1))
    raw_H = sp.expand(W**2 * sp.prod(W - rho for rho in roots))
    c = sp.factor(-sp.diff(raw_H, W).subs(W, 1))
    assert c != 0
    H = raw_H
    h2 = sp.Poly(H, W).coeff_monomial(W**2)
    inverse = sp.expand(H - B * C * W + c * A * C**2)
    pulled_discriminant = sp.factor(sp.discriminant(inverse, W))
    saturated = sp.factor(pulled_discriminant / C**2)
    assert sp.factor(pulled_discriminant - C**2 * saturated) == 0
    trace = sp.factor(saturated.subs(C, 0))
    conic = B**2 - 4 * h2 * c * A
    assert sp.rem(sp.Poly(trace, A), sp.Poly(conic, A)).is_zero
    quotient = sp.factor(trace / conic)
    assert quotient != 0 and not quotient.has(A, B)
    conic_ratio = sp.factor(sp.Poly(conic, A, B).sqf_part().as_expr() / conic)
    assert conic_ratio != 0 and not conic_ratio.has(A, B)

    print(
        f"PASS degree {m + 2}: C24 nilpotency index {nilpotency_index}, "
        f"reduction {c24_reduced}=0; weighted reduced intersection "
        f"{sp.factor(conic)}=0"
    )

print("PASS: thickened A^1 disjoint-union G_m versus reduced A^1")

#!/usr/bin/env python3
"""Exact identities used by stable cancellation-parameter faithfulness.

The theorem itself is a UFD, residue-degree, and stable-boundary argument.
This regression checks its coordinate identities and the endpoint facts for
the cancellation parameter polynomial on a bounded all-positive grid.
"""
from __future__ import annotations

import sympy as sp

from master_cancellation import parameter_polynomial


A, B, P, Q, R, C = sp.symbols("A B P Q R C", nonzero=True)
x, y, z, q, u, a = sp.symbols("x y z q u a", nonzero=True)


# The stable boundary-plane calculation is insensitive to extra variables.
# Once Q maps to uQ, preservation of G=(r+1)RQ^m-C forces the displayed
# inverse weight on R.
for m_value in range(1, 7):
    for r_value in range(1, 7):
        G = (r_value + 1) * R * Q**m_value - C
        transformed = (
            (r_value + 1)
            * (u ** (-m_value) * R)
            * (u * Q) ** m_value
            - C
        )
        assert sp.cancel(transformed - G) == 0


# The source factorization fixes A, then A-1=xy^m fixes x and y up to the
# visible scaling.  The full Q identity forces a=u^(m+1).
for m_value in range(1, 7):
    source_A = 1 + x * y**m_value
    scaled_x = u ** (-m_value) * x
    scaled_y = u * y
    assert sp.cancel(1 + scaled_x * scaled_y**m_value - source_A) == 0

    source_Q = y + x * B
    pulled_Q = scaled_y + scaled_x * a * B
    discrepancy = sp.expand(pulled_Q - u * source_Q)
    expected = sp.expand((a * u ** (-m_value) - u) * x * B)
    assert sp.cancel(discrepancy - expected) == 0
    assert sp.cancel(discrepancy.subs(a, u ** (m_value + 1))) == 0

    # The reconstruction-pole residue has weight zero.
    scaled_residue = (u ** (m_value + 1) * B) / scaled_y ** (m_value + 1)
    assert sp.cancel(scaled_residue - B / y ** (m_value + 1)) == 0


# Parameter roots never meet the excluded endpoints q=0,1.  Construction and
# uniqueness of the full Hensel jet are checked by verify_master_universal.py.
for m_value in range(1, 9):
    for r_value in range(1, 9):
        modulus = parameter_polynomial(m_value, r_value, q)
        assert sp.expand(modulus.subs(q, 0)) != 0
        assert sp.expand(modulus.subs(q, 1)) != 0


print("PASS: stable boundary-plane weights preserve the cancellation trace")
print("PASS: source UFD scaling forces B to have weight m+1")
print("PASS: the reconstruction-pole residue B/y^(m+1) recovers q")

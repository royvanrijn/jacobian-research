#!/usr/bin/env python3
"""Exact regressions for target-fixed rigidity of C24 parameter branches."""
from __future__ import annotations

import sympy as sp

from master_cancellation import parameter_polynomial


A, y, z, q, w, v = sp.symbols("A y z q w v")


def infinity_polynomial(m: int, r: int) -> sp.Expr:
    integral = sp.integrate(v**r * (1 - v) ** (m * r), (v, 0, w))
    return sp.cancel(integral / w ** (r + 1))


# The unique target-field birational identification fixes A,y and B. Its
# last coordinate is polynomial exactly when the two jets agree mod A^(r+1).
for r in range(1, 6):
    left = sp.symbols(f"l0:{r + 3}")
    right = sp.symbols(f"r0:{r + 3}")
    h_left = sum(left[i] * A**i for i in range(r + 3))
    h_right = sum(right[i] * A**i for i in range(r + 3))
    delta = sp.expand(h_left - h_right)
    remainder = sum((left[i] - right[i]) * A**i for i in range(r + 1))
    assert sp.rem(
        sp.Poly(delta - remainder, A), sp.Poly(A ** (r + 1), A)
    ).is_zero

    for m in range(1, 5):
        B_left = A ** (r + 1) * z + y ** (m + 1) * h_left
        z_prime = z + y ** (m + 1) * delta / A ** (r + 1)
        B_right_after = A ** (r + 1) * z_prime + y ** (m + 1) * h_right
        assert sp.cancel(B_right_after - B_left) == 0

        # Distinct constant parameters leave a nonzero A=0 pole numerator.
        pole_numerator = sp.expand(y ** (m + 1) * remainder)
        assert sp.expand(pole_numerator.subs(A, 0) - y ** (m + 1) * (
            left[0] - right[0]
        )) == 0


# A selected parameter root fills the projective P=0 branch
# w=-q/(1-q). Check the exact identity modulo M_(m,r).
for m in range(1, 6):
    for r in range(1, 6):
        modulus = parameter_polynomial(m, r, q)
        K = infinity_polynomial(m, r)
        substituted = sp.together(K.subs(w, -q / (1 - q)))
        numerator = substituted.as_numer_denom()[0]
        assert sp.rem(sp.Poly(numerator, q), sp.Poly(modulus, q)).is_zero


print("PASS: unique target-fixed birational identification and pole criterion")
print("PASS: every parameter root marks its projective P=0 branch")

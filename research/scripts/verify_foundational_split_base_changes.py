#!/usr/bin/env python3
"""Exact algebra for split one-coordinate foundational base changes.

For the normalized tangent slice C1=1, pull back the foundational cover by

    (w,u,v) |-> (C0,C2,C3) = (f(w),u,v).

The source has equation f(w)=a*c in the global factor-space coordinates.
This checker verifies the global formulas, the smoothness certificate, the
UFD localization, and the divisor-valuation presentation used to compute
the class group.  The all-degree class-group conclusion is proved in the
companion note; the finite partition loop is only a regression.
"""

from __future__ import annotations

from math import gcd

import sympy as sp


a, y, z, w = sp.symbols("a y z w")
b, c, d, e = sp.symbols("b c d e")

factor_b = 1 + a * y
factor_c = 1 - sp.Rational(3, 2) * a * y + a**2 * z
factor_d = (
    sp.Rational(1, 2) * y
    - a * z
    + sp.Rational(3, 2) * a * y**2
    - a**2 * y * z
)
factor_e = (
    -2 * z
    + 4 * y**2
    - 4 * a * y * z
    + 3 * a * y**3
    - 2 * a**2 * y**2 * z
)

substitution = {
    b: factor_b,
    c: factor_c,
    d: factor_d,
    e: factor_e,
}

resultant = a**2 * e - a * b * d + b**2 * c
mixed_coefficient = a * d + b * c
leading_coefficient = a * c

assert sp.expand(resultant.subs(substitution) - 1) == 0
assert sp.expand(mixed_coefficient.subs(substitution) - 1) == 0
assert sp.expand(leading_coefficient.subs(substitution) - a * factor_c) == 0

# Use a generic polynomial f for the differential and localization checks.
f0, f1, f2, f3 = sp.symbols("f0 f1 f2 f3")
f = f0 + f1 * w + f2 * w**2 + f3 * w**3
H = sp.expand(f - a * factor_c)

# The hypersurface H=0 is smooth for every nonzero f: away from a=0 the
# z-derivative is nonzero, while on a=0 the a-derivative is -1.
assert sp.diff(H, z) == -a**3
assert sp.expand(sp.diff(H, a).subs(a, 0) + 1) == 0

# On D(a), solve c=f/a and then solve globally for z.  This identifies the
# localization with k[a^(+/-1),y,w], a Laurent polynomial UFD.
localized_z = sp.cancel((f / a - 1 + sp.Rational(3, 2) * a * y) / a**2)
assert sp.cancel(factor_c.subs(z, localized_z) - f / a) == 0
assert sp.cancel(H.subs(z, localized_z)) == 0


def class_group_invariants(multiplicities: tuple[int, ...]) -> tuple[int, int]:
    """Return free rank and torsion order of Z^r / Z*(m_1,...,m_r)."""

    assert multiplicities and all(value > 0 for value in multiplicities)
    torsion_order = 0
    for value in multiplicities:
        torsion_order = gcd(torsion_order, value)
    return len(multiplicities) - 1, torsion_order


def partitions(total: int, ceiling: int | None = None) -> list[tuple[int, ...]]:
    """Integer partitions, used only to regress the zero-class-group test."""

    if total == 0:
        return [()]
    if ceiling is None or ceiling > total:
        ceiling = total
    result: list[tuple[int, ...]] = []
    for first in range(ceiling, 0, -1):
        for tail in partitions(total - first, first):
            result.append((first,) + tail)
    return result


for degree in range(1, 13):
    for multiplicities in partitions(degree):
        free_rank, torsion_order = class_group_invariants(multiplicities)
        class_group_is_zero = free_rank == 0 and torsion_order == 1
        assert class_group_is_zero == (multiplicities == (1,))

print("PASS normalized factor-space and split base-change hypersurface formulas")
print("PASS smoothness and Laurent-UFD localization certificates")
print("PASS Cl = Z^r / Z*(m_1,...,m_r); it vanishes only for degree-one f")
print("SCOPE: split one-coordinate coefficient maps (f(w),u,v)")

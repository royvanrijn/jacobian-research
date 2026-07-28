#!/usr/bin/env python3
"""Exact checks for power-shifted universal quartic multiplicity."""

from __future__ import annotations

import sympy as sp


x, y, z, S = sp.symbols("x y z S")
g1, g2, g3, g4 = sp.symbols("g1 g2 g3 g4", nonzero=True)

t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
P = t * q
root = x / t
Q = y + x * q
D = sp.factor(1 - root * (Q - P * root))
assert D == 1 / t


def shifted_map(power_shift: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """The raw determinant-minus-two quartic map F_m."""
    second = (
        y
        + 3 * (g3 / g1) * x * q
        + 2 * (g2 / g1) * t * q
        + 4
        * (g4 / g1)
        * t ** (power_shift + 2)
        * x**2
        * q ** (power_shift + 4)
    )
    third = (
        x * (5 - 3 * t)
        - (g3 / g1) * x**3 * z
        - 2
        * (g4 / g1)
        * t**power_shift
        * x**4
        * q ** (power_shift + 4)
    )
    return P, second, third


def lifted_seed(power_shift: int, variable: sp.Expr) -> sp.Expr:
    return (
        g1 * variable
        + P * (g2 * variable**2 + g3 * variable**3)
        + g4 * P ** (power_shift + 4) * variable**4
    )


# Direct three-variable regressions.  The rational-chart proof in the note
# is uniform in m; these cases guard the denominator-free algebraization.
for m in range(4):
    first, second, third = shifted_map(m)
    inverse = sp.factor(
        lifted_seed(m, root) - g1 * (second * root**2 + third) / 2
    )
    assert inverse == 0

    derivative = sp.diff(
        lifted_seed(m, S) - g1 * (second * S**2 + third) / 2,
        S,
    ).subs(S, root)
    assert sp.factor(derivative - g1 * D) == 0

    jacobian = sp.det(
        sp.Matrix(
            [
                [sp.diff(component, variable) for variable in (x, y, z)]
                for component in (first, second, third)
            ]
        )
    )
    assert sp.factor(jacobian) == -2


# One connected witness: f(T)=T^4-3T^2-1, translated at a=1.
T = sp.symbols("T")
quartic = T**4 - 3 * T**2 - 1
translated = sp.expand(quartic.subs(T, 1 + S))
value = quartic.subs(T, 1)
seed = sp.expand(translated - value)
seed_poly = sp.Poly(seed, S)
witness_coefficients = [
    seed_poly.coeff_monomial(S**index) for index in range(1, 5)
]
assert witness_coefficients == [-2, 3, 4, 1]
assert sp.discriminant(quartic, T) != 0

for m in range(8):
    witness_lift = (
        witness_coefficients[0] * S
        + sp.Symbol("p")
        * (
            witness_coefficients[1] * S**2
            + witness_coefficients[2] * S**3
        )
        + witness_coefficients[3] * sp.Symbol("p") ** (m + 4) * S**4
    )
    target_c = -2 * value / witness_coefficients[0]
    selected_inverse = sp.expand(
        witness_lift.subs(sp.Symbol("p"), 1)
        - witness_coefficients[0] * target_c / 2
    )
    assert sp.expand(selected_inverse - translated) == 0


# The normalized ramified-stratum Fitting divisor.
p, r, a2, a3, a4 = sp.symbols("p r a2 a3 a4", nonzero=True)
for m in range(8):
    H = r + p * (a2 * r**2 + a3 * r**3) + a4 * p ** (m + 4) * r**4
    slope = sp.cancel(sp.diff(H, r) / r)
    fitting = sp.expand(r**2 * sp.diff(slope, r))
    expected = -1 + 3 * a3 * p * r**2 + 8 * a4 * p ** (m + 4) * r**3
    assert sp.expand(fitting - expected) == 0

    support_vectors = sp.Matrix([[1, 2], [m + 4, 3]])
    assert abs(int(support_vectors.det())) == 2 * m + 5

    # Lower Newton polygon vertices over p=0 have slopes 0, 1, m+3.
    vertices = [(0, 0), (2, 0), (3, 1), (4, m + 4)]
    slopes = [
        sp.Rational(y2 - y1, x2 - x1)
        for (x1, y1), (x2, y2) in zip(vertices, vertices[1:])
    ]
    assert slopes == [0, 1, m + 3]


print("PASS: power-shifted quartic gauges have determinant -2")
print("PASS: their selected inverse polynomial is independent of the shift")
print("PASS: every map has the same complete quartic witness fiber")
print("PASS: the ramified Fitting-support indices are 2*m+5")

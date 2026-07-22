#!/usr/bin/env python3
"""Verify the exact linear normalization used in Long's BCW discussion."""

import sympy as sp


x, y, z = sp.symbols("x y z")
u = 1 + x * y

# Canonical determinant -2 presentation used throughout this repository.
F = sp.Matrix(
    [
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    ]
)

# Determinant-one presentation transcribed from arXiv:2607.18186v1.
v = 1 + 2 * x * y
L = sp.Matrix(
    [
        v**3 * z + 4 * y**2 * v * (2 + 3 * x * y),
        y + 3 * x * v**2 * z + 12 * x * y**2 * (2 + 3 * x * y),
        -x + 3 * x**2 * y + x**3 * z,
    ]
)

source_substitution = {x: x, y: 2 * y, z: 2 * z}
target_scaling = sp.diag(sp.Rational(1, 2), sp.Rational(1, 2), -sp.Rational(1, 2))
assert sp.simplify(L - target_scaling * F.subs(source_substitution)) == sp.zeros(3, 1)
assert sp.factor(L.jacobian((x, y, z)).det()) == 1

# The repository's three-point collision transports, including the third point
# not displayed in Long's introductory two-point list.
points = [
    (0, 0, -sp.Rational(1, 8)),
    (1, -sp.Rational(3, 4), sp.Rational(13, 4)),
    (-1, sp.Rational(3, 4), sp.Rational(13, 4)),
]
target = sp.Matrix([-sp.Rational(1, 8), 0, 0])
for point in points:
    assert sp.simplify(L.subs(dict(zip((x, y, z), point)))) == target

print("PASS Long normalization: L = diag(1/2,1/2,-1/2) o F o diag(1,2,2)")
print("PASS Long normalization: determinant 1 and transported three-point collision")

#!/usr/bin/env python3
"""Exact checks for the foundational map's p-adic inverse-branch note.

The symbolic part verifies the two inverse charts, discriminant, omitted
triple-root curve, and the integral collision obstructing a global Tate
inverse.  The finite-field part checks the rational-root/fiber formula at
three good primes; the geometric 3/1/0 tube proof is in the note.
"""

from __future__ import annotations

from collections import Counter

import sympy as sp


x, y, z, a, b, c, t, s, T = sp.symbols("x y z a b c t s T")

u = 1 + x * y
F = (
    u**3 * z + y**2 * u * (4 + 3 * x * y),
    y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)

jacobian = sp.factor(sp.det(sp.Matrix(F).jacobian((x, y, z))))
assert jacobian == -2

P = c * T**3 - 2 * T**2 + b * T - 2 * a
Q = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
assert sp.factor(sp.discriminant(P, T) + 4 * Q) == 0

# Finite-root reconstruction, with a eliminated by P(t)=0.
a_on_root = (c * t**3 - 2 * t**2 + b * t) / 2
r = 3 * c * t**2 - 4 * t + b
finite_source = (
    2 / r,
    t - r / 2,
    5 * r**2 / 4 - 3 * t * r / 2 - c * r**3 / 8,
)
finite_image = [
    sp.factor(expr.subs(dict(zip((x, y, z), finite_source))) - target)
    for expr, target in zip(F, (a_on_root, b, c))
]
assert finite_image == [0, 0, 0]

# The projective-root chart s=1/t is regular at s=0.
D = 1 - b * s + 3 * a * s**2
c_projective = 2 * s - b * s**2 + 2 * a * s**3
z_projective = sp.cancel(
    (
        5 * D**2
        - 3 * D
        - (2 - b * s + 2 * a * s**2) * D**3
    )
    / s**2
)
assert sp.denom(z_projective) == 1
assert sp.expand(z_projective.subs(s, 0)) == a - 4 * b**2
projective_source = (s / D, b - 3 * a * s, z_projective)
projective_image = [
    sp.factor(
        (expr.subs(dict(zip((x, y, z), projective_source))) - target)
        .subs(c, c_projective)
    )
    for expr, target in zip(F, (a, b, c))
]
assert projective_image == [0, 0, 0]

# Gamma is the triple-root (hence omitted) curve.
gamma = {a: 4 / (27 * c**2), b: 4 / (3 * c)}
assert sp.factor(P.subs(gamma) - c * (T - 2 / (3 * c)) ** 3) == 0
assert sp.factor(Q.subs(gamma)) == 0

# The exact rational collision lies in the integral unit ball at every odd p.
collision_target = (sp.Rational(-1, 4), 0, 0)
collision_sources = (
    (0, 0, sp.Rational(-1, 4)),
    (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
)
for source in collision_sources:
    assert tuple(sp.expand(expr.subs(dict(zip((x, y, z), source)))) for expr in F) == (
        collision_target
    )

# The origin branch has an open-unit implicit equation, while Gamma supplies
# a target on the boundary of the centered closed unit polydisc for p >= 5.
assert sp.diff(c_projective, s).subs(s, 0) == 2
gamma_boundary = (sp.Rational(4, 27), sp.Rational(4, 3), 1)
assert gamma_boundary[1] != 0 and gamma_boundary[2] != 0
assert tuple(
    sp.factor(expr)
    for expr in (
        3 * gamma_boundary[1] * gamma_boundary[2] - 4,
        12 * gamma_boundary[0] - gamma_boundary[1] ** 2,
    )
) == (0, 0)


def foundational_map_mod_p(point: tuple[int, int, int], p: int) -> tuple[int, int, int]:
    xx, yy, zz = point
    uu = (1 + xx * yy) % p
    return (
        (uu**3 * zz + yy**2 * uu * (4 + 3 * xx * yy)) % p,
        (yy + 3 * xx * uu**2 * zz + 3 * xx * yy**2 * (4 + 3 * xx * yy)) % p,
        (2 * xx - 3 * xx**2 * yy - xx**3 * zz) % p,
    )


def q_mod_p(point: tuple[int, int, int], p: int) -> int:
    aa, bb, cc = point
    return (
        27 * aa**2 * cc**2
        - 18 * aa * bb * cc
        + 16 * aa
        + bb**3 * cc
        - bb**2
    ) % p


def on_gamma_mod_p(point: tuple[int, int, int], p: int) -> bool:
    aa, bb, cc = point
    return (3 * bb * cc - 4) % p == 0 and (12 * aa - bb**2) % p == 0


for prime in (5, 7, 11):
    fiber_sizes: Counter[tuple[int, int, int]] = Counter()
    for xx in range(prime):
        for yy in range(prime):
            for zz in range(prime):
                fiber_sizes[foundational_map_mod_p((xx, yy, zz), prime)] += 1

    strata = Counter()
    for aa in range(prime):
        for bb in range(prime):
            for cc in range(prime):
                target = (aa, bb, cc)
                size = fiber_sizes[target]
                simple_roots = 0
                for tt in range(prime):
                    value = (cc * tt**3 - 2 * tt**2 + bb * tt - 2 * aa) % prime
                    derivative = (3 * cc * tt**2 - 4 * tt + bb) % prime
                    simple_roots += value == 0 and derivative != 0
                expected_rational_size = simple_roots + (cc == 0)
                assert size == expected_rational_size
                if q_mod_p(target, prime) != 0:
                    strata["off_discriminant"] += 1
                elif on_gamma_mod_p(target, prime):
                    assert size == 0
                    strata["omitted_gamma"] += 1
                else:
                    strata["discriminant_survivor"] += 1
    assert sum(strata.values()) == prime**3

print(
    "PASS p-adic inverse branches: charts, discriminant, Gamma obstruction, "
    "collision, and good-prime rational-root/fiber regressions"
)

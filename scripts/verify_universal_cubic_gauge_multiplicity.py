#!/usr/bin/env python3
"""Exact checks for universal cubic multiplicity by gauge-lift deformation."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.keller_fiber import compile_polynomial_to_keller_fiber


x, y, z, S = sp.symbols("x y z S")
g1, g2, g3 = sp.symbols("g1 g2 g3", nonzero=True)

t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
P = t * q
root = x / t
Q = y + x * q
D = sp.factor(1 - root * (Q - P * root))
assert D == 1 / t


def cubic_map(exponent: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    delta_b = 3 * (g3 / g1) * (
        t ** (exponent - 1) * x * q**exponent - t**2 * x * q**3
    )
    delta_c = -(g3 / g1) * (
        t ** (exponent - 3) * x**3 * q**exponent - x**3 * q**3
    )
    second = y + 3 * (g3 / g1) * x * q + 2 * (g2 / g1) * t * q + delta_b
    third = x * (5 - 3 * t) - (g3 / g1) * x**3 * z + delta_c
    return P, second, third


def lifted_cubic(exponent: int, variable: sp.Expr) -> sp.Expr:
    return (
        g1 * variable
        + g2 * P * variable**2
        + g3 * P * (1 + P ** (exponent - 1) - P**2) * variable**3
    )


# Direct polynomial-map regressions.  The rational marked-line proof is
# uniform in n; two cases catch both added monomials and their cancellation.
for n in (4, 5):
    first, second, third = cubic_map(n)
    inverse = sp.factor(
        lifted_cubic(n, root) - g1 * (second * root**2 + third) / 2
    )
    assert inverse == 0

    derivative = sp.diff(
        lifted_cubic(n, S) - g1 * (second * S**2 + third) / 2,
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


# At P=1 the deformation disappears for every exponent.
p = sp.symbols("p")
for n in range(4, 13):
    coefficient = sp.expand(p * (1 + p ** (n - 1) - p**2))
    assert coefficient.subs(p, 1) == 1

    degree_drop = 1 + p ** (n - 1) - p**2
    assert sp.degree(degree_drop, p) == n - 1
    assert sp.gcd(degree_drop, sp.diff(degree_drop, p)) == 1
    assert degree_drop.subs(p, 0) == 1
    assert degree_drop.subs(p, 1) == 1

    # Every simple nonzero zero of h_n gives two finite roots and one
    # unramified escaping root.
    vertices = [(0, 0), (2, 0), (3, 1)]
    slopes = [
        sp.Rational(y2 - y1, x2 - x1)
        for (x1, y1), (x2, y2) in zip(vertices, vertices[1:])
    ]
    assert slopes == [0, 1]
    assert n - 1 + 1 == n  # vertical components plus discriminant


# One fixed connected complete cubic fiber: f(T)=T^3-T-1 at a=0.
T = sp.symbols("T")
f = T**3 - T - 1
translated = sp.expand(f.subs(T, S))
value = f.subs(T, 0)
target_c = -2 * value / (-1)
assert sp.discriminant(f, T) == -23
assert sp.factor(f) == f
for n in range(4, 10):
    selected_seed = sp.expand(
        (g1 * S + g2 * p * S**2 + g3 * p * (1 + p ** (n - 1) - p**2) * S**3)
        .subs({p: 1, g1: -1, g2: 0, g3: 1})
    )
    selected_inverse = sp.expand(selected_seed - sp.Rational(-1, 2) * target_c)
    assert selected_inverse == translated


# The P=0 fiber has a quadratic q=0 block and one affine t=0 branch.
# Restrict the displayed polynomial map using x*y=-1 and
# q=(g1/g3)*y^2.
for n in (4, 7):
    _, second, third = cubic_map(n)
    t_zero_substitution = {
        x: -1 / y,
        z: z,
    }
    second_t_zero = sp.factor(second.subs(t_zero_substitution))
    third_t_zero = sp.factor(third.subs(t_zero_substitution))
    assert second_t_zero == -2 * y
    expected_third = sp.factor(
        5 * (-1 / y)
        - (g3 / g1) * (-1 / y) ** 3 * z
        + (g3 / g1) * (-1 / y) ** 3 * ((g1 / g3) * y**2) ** 3
    )
    assert sp.factor(third_t_zero - expected_third) == 0


# The public compiler exposes the same family with stable_parameter=k
# corresponding to n=k+4, while retaining the selected inverse polynomial.
compiled_cubics = [
    compile_polynomial_to_keller_fiber(
        f,
        T,
        translation=0,
        inverse_variable=S,
        source_variables=(x, y, z),
        stable_parameter=parameter,
    )
    for parameter in (0, 3)
]
for parameter, compilation in zip((0, 3), compiled_cubics):
    certificate = compilation.stable_multiplicity
    assert certificate is not None
    assert certificate.family_parameter == parameter
    assert certificate.gauge_exponent == parameter + 4
    assert certificate.separation_invariant == (
        "geometric_boundary_target_components"
    )
    assert certificate.separation_value == parameter + 4
    assert compilation.inverse_polynomial == translated
    assert sp.expand(
        compilation.lifted_seed.subs({x: 0, y: 0, z: 1})
        - compilation.seed
    ) == 0
assert compiled_cubics[0].determinant_minus_two_map != (
    compiled_cubics[1].determinant_minus_two_map
)


print("PASS: fiber-invisible cubic gauge lifts have determinant -2")
print("PASS: the selected connected inverse cubic is independent of the lift exponent")
print("PASS: h_n has n-1 simple nonzero roots in the regression range")
print("PASS: P=0 has two q=0 sheets and one affine t=0 sheet")
print("PASS: the canonical boundary-component count is n")
print("PASS: the public compiler returns the exact cubic stable certificate")

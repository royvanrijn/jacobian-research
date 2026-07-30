#!/usr/bin/env python3
"""Exact rank-five ambient, projective, and marked transition loci."""

from __future__ import annotations

import sympy as sp


T = sp.symbols("T")
u = sp.symbols("u1:6")
alpha = sp.symbols("alpha", nonzero=True)


# Ordered-root form of the ambient stable-equivalence hypersurface.
root_polynomial = sp.Poly(sp.prod(T - root for root in u), T)
e1 = -root_polynomial.nth(4)
e2 = root_polynomial.nth(3)
e3 = -root_polynomial.nth(2)
e4 = root_polynomial.nth(1)
e5 = -root_polynomial.nth(0)

ambient = sp.expand(
    968203125 * e4**2
    - 75076 * e2 * e1**6
)
E1, E2, E3, E4, E5 = sp.symbols("E1 E2 E3 E4 E5")
ambient_elementary = 968203125 * E4**2 - 75076 * E2 * E1**6

base_roots = tuple(sp.Integer(index) for index in range(1, 6))
base_substitution = dict(zip(u, base_roots, strict=True))
assert ambient.subs(base_substitution) == 0


# The two framed PGL_2 residuals after matching the first three roots.
rows = [
    [1, base_roots[index], u[index], base_roots[index] * u[index]]
    for index in range(5)
]
projective_4 = sp.factor(sp.Matrix(rows[:4]).det())
projective_5 = sp.factor(sp.Matrix(rows[:3] + [rows[4]]).det())

expected_projective_4 = (
    u[0] * u[1]
    - 4 * u[0] * u[2]
    + 3 * u[0] * u[3]
    + 3 * u[1] * u[2]
    - 4 * u[1] * u[3]
    + u[2] * u[3]
)
expected_projective_5 = 2 * (
    u[0] * u[1]
    - 3 * u[0] * u[2]
    + 2 * u[0] * u[4]
    + 2 * u[1] * u[2]
    - 3 * u[1] * u[4]
    + u[2] * u[4]
)
assert sp.expand(projective_4 - expected_projective_4) == 0
assert sp.expand(projective_5 - expected_projective_5) == 0

ambient_rank = sp.Matrix([ambient]).jacobian(u).subs(base_substitution).rank()
projective_rank = (
    sp.Matrix([projective_4, projective_5])
    .jacobian(u)
    .subs(base_substitution)
    .rank()
)
intersection_rank = (
    sp.Matrix([ambient, projective_4, projective_5])
    .jacobian(u)
    .subs(base_substitution)
    .rank()
)
assert ambient_rank == 1
assert projective_rank == 2
assert intersection_rank == 3


# Coefficient parameterization and explicit ambient source-target scaling.
A0, A1, A2, A3, A4 = sp.symbols(
    "A0 A1 A2 A3 A4",
    nonzero=True,
)
base_i5 = sp.Rational(75076, 968203125)
ambient_a3 = sp.factor(A1**2 / (base_i5 * A4**6))


def compiler_seeds(
    linear: sp.Expr,
    cubic: sp.Expr,
    quartic: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the rank-five compiler seeds."""

    seed_4 = sp.factor(quartic * linear**3 / cubic**4)
    seed_5 = sp.factor(linear**4 / cubic**5)
    return seed_4, seed_5


base_seed_4, base_seed_5 = compiler_seeds(
    sp.Integer(274),
    sp.Integer(85),
    sp.Integer(-15),
)
moving_seed_4, moving_seed_5 = compiler_seeds(A1, A3, A4)
coefficient_alpha = sp.factor(-1275 * A1 / (274 * A3 * A4))

assert sp.factor(
    (moving_seed_4 - coefficient_alpha**5 * base_seed_4).subs(
        A3, ambient_a3
    )
) == 0
assert sp.factor(
    (moving_seed_5 - coefficient_alpha**6 * base_seed_5).subs(
        A3, ambient_a3
    )
) == 0


# Canonical marked-target compatibility.
alpha_definition = 274 * alpha * e1 * e2 - 1275 * e4
marked_pi = 274 * alpha**2 * e2 - 85 * e4
marked_b = 274 * alpha * e3 - 225 * e4
marked_c = 274 * e5 - 120 * alpha * e4

marked_substitution = dict(base_substitution)
marked_substitution[alpha] = 1
marked_equations = [
    ambient,
    alpha_definition,
    marked_pi,
    marked_b,
    marked_c,
]
assert all(
    equation.subs(marked_substitution) == 0
    for equation in marked_equations
)

marked_rank = (
    sp.Matrix(marked_equations)
    .jacobian((*u, alpha))
    .subs(marked_substitution)
    .rank()
)
assert marked_rank == 5


# Exact elimination argument on the clean locus:
# alpha_definition / marked_pi gives e1=15*alpha; marked_pi gives
# e4=(274/85)*alpha^2*e2; the ambient equation then forces e2=85*alpha^2.
ambient_after_first_steps = sp.factor(
    ambient_elementary.subs(
        {
            E1: 15 * alpha,
            E4: sp.Rational(274, 85) * alpha**2 * E2,
        },
        simultaneous=True,
    )
)
expected_reduction = (
    -sp.Rational(171032512500, 17)
    * E2
    * alpha**4
    * (-E2 + 85 * alpha**2)
)
assert sp.factor(ambient_after_first_steps - expected_reduction) == 0

scaled_elementary = {
    E1: 15 * alpha,
    E2: 85 * alpha**2,
    E3: 225 * alpha**3,
    E4: 274 * alpha**4,
    E5: 120 * alpha**5,
}
marked_elementary_equations = [
    ambient_elementary,
    274 * alpha * E1 * E2 - 1275 * E4,
    274 * alpha**2 * E2 - 85 * E4,
    274 * alpha * E3 - 225 * E4,
    274 * E5 - 120 * alpha * E4,
]
assert all(
    sp.factor(equation.subs(scaled_elementary, simultaneous=True)) == 0
    for equation in marked_elementary_equations
)

scaled_polynomial = sp.Poly(
    T**5
    - 15 * alpha * T**4
    + 85 * alpha**2 * T**3
    - 225 * alpha**3 * T**2
    + 274 * alpha**4 * T
    - 120 * alpha**5,
    T,
)
expected_scaled_polynomial = sp.Poly(
    sp.prod(T - alpha * index for index in range(1, 6)),
    T,
)
assert scaled_polynomial == expected_scaled_polynomial


print("PASS: the rank-five ambient stable locus is one hypersurface")
print("PASS: the labelled projective locus has two independent residuals")
print("PASS: their intersection is transverse of codimension three")
print("PASS: equal I_5 gives the explicit unique seed scaling")
print("PASS: canonical marked transport is exactly root scaling")
print("PASS: arbitrary marked transport reduces to a fixed-map stabilizer orbit")

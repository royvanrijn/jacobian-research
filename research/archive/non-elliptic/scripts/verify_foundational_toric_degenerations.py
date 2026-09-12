#!/usr/bin/env python3
"""Exact diagonal one-parameter degenerations of the foundational cubic.

This is a narrow orbit-closure certificate.  It classifies diagonal source
weights in the displayed foundational coordinates, with the unique output
renormalization which keeps a nonzero limit.  Constant linear changes before
the one-parameter subgroup and nonlinear polynomial degenerations are not
covered.
"""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
a, b, c = sp.symbols("a b c")
X, Y, Z = sp.symbols("X Y Z")

u = 1 + x * y
F = (
    sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y)),
    sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)),
    sp.expand(2 * x - 3 * x**2 * y - x**3 * z),
)

assert sp.expand(
    sp.det(sp.Matrix([[sp.diff(fi, v) for v in (x, y, z)] for fi in F]))
) == -2

# The linear terms z,y,2x have weights c,b,a.  If m_i is the least
# monomial weight in F_i, then m_i <= (c,b,a)_i.  Preservation of the
# nonzero constant Jacobian requires
#
#     m_1 + m_2 + m_3 = a + b + c,
#
# so all three inequalities are equalities separately.  It remains only to
# classify when every monomial lies above the corresponding linear term.
linear_exponents = ((0, 0, 1), (0, 1, 0), (1, 0, 0))
differences: list[list[tuple[tuple[int, int, int], sp.Expr, sp.Expr]]] = []
for fi, linear_exp in zip(F, linear_exponents):
    rows = []
    linear_weight = sum(e * w for e, w in zip(linear_exp, (a, b, c)))
    for exponent, coefficient in sp.Poly(fi, x, y, z).terms():
        monomial_weight = sum(e * w for e, w in zip(exponent, (a, b, c)))
        rows.append((exponent, coefficient, sp.expand(monomial_weight - linear_weight)))
    differences.append(rows)

# The entire admissible cone is the quadrant p,q >= 0 times the lineality
# direction (-1,1,2):
#
#     p = 2b-c,  q = a-b+c,  and a+b=p+q.
p = 2 * b - c
q = a - b + c
assert sp.expand(a + b - p - q) == 0

p_vector = sp.Matrix((0, 2, -1))
q_vector = sp.Matrix((1, -1, 1))
coefficient_matrix = sp.Matrix.hstack(p_vector, q_vector)

decompositions: dict[tuple[int, int, int], tuple[int, int]] = {}
for rows in differences:
    for exponent, _, difference in rows:
        diff_vector = sp.Matrix(
            (
                sp.diff(difference, a),
                sp.diff(difference, b),
                sp.diff(difference, c),
            )
        )
        alpha, beta = tuple(
            next(iter(sp.linsolve((coefficient_matrix, diff_vector))))
        )
        assert alpha.is_Integer and beta.is_Integer
        assert int(alpha) >= 0 and int(beta) >= 0
        assert sp.expand(difference - alpha * p - beta * q) == 0
        decompositions[exponent] = (int(alpha), int(beta))

nonzero_difference_pairs = {
    (int(alpha), int(beta))
    for rows in differences
    for _, _, difference in rows
    for alpha, beta in [
        tuple(
            next(
                iter(
                    sp.linsolve(
                        (
                            coefficient_matrix,
                            sp.Matrix(
                                (
                                    sp.diff(difference, a),
                                    sp.diff(difference, b),
                                    sp.diff(difference, c),
                                )
                            ),
                        )
                    )
                )
            )
        )
    ]
    if alpha != 0 or beta != 0
}
assert nonzero_difference_pairs == {
    (0, 1),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 1),
    (2, 2),
    (2, 3),
    (3, 2),
    (3, 3),
}


def initial_map(p_zero: bool, q_zero: bool) -> tuple[sp.Expr, ...]:
    """Return the initial map on one relatively open face of the quadrant."""

    result = []
    for rows in differences:
        terms = []
        for exponent, coefficient, difference in rows:
            diff_vector = sp.Matrix(
                (
                    sp.diff(difference, a),
                    sp.diff(difference, b),
                    sp.diff(difference, c),
                )
            )
            alpha, beta = tuple(
                next(iter(sp.linsolve((coefficient_matrix, diff_vector))))
            )
            survives = (p_zero or alpha == 0) and (q_zero or beta == 0)
            if survives:
                monomial = coefficient
                for variable, power in zip((x, y, z), exponent):
                    monomial *= variable**power
                terms.append(monomial)
        result.append(sp.expand(sum(terms)))
    return tuple(result)


initials = {
    "interior": initial_map(False, False),
    "p_face": initial_map(True, False),
    "q_face": initial_map(False, True),
    "lineality": initial_map(True, True),
}

assert initials["interior"] == (z, y, 2 * x)
assert initials["p_face"] == (4 * y**2 + z, y, 2 * x)
assert initials["q_face"] == (z, 3 * x * z + y, 2 * x)
assert initials["lineality"] == F

inverses = {
    "interior": (Z / 2, Y, X),
    "p_face": (Z / 2, Y, X - 4 * Y**2),
    "q_face": (Z / 2, Y - sp.Rational(3, 2) * X * Z, X),
}

for name, initial in initials.items():
    jacobian = sp.expand(
        sp.det(
            sp.Matrix(
                [[sp.diff(fi, variable) for variable in (x, y, z)] for fi in initial]
            )
        )
    )
    assert jacobian == -2
    if name == "lineality":
        continue
    inverse = inverses[name]
    forward_after_inverse = tuple(
        sp.expand(fi.subs({x: inverse[0], y: inverse[1], z: inverse[2]}))
        for fi in initial
    )
    assert forward_after_inverse == (X, Y, Z)

print("PASS foundational diagonal toric degeneration cone: p>=0, q>=0")
print("PASS its four faces are the foundational map or three automorphisms")
print("SCOPE: displayed diagonal weights only; no oblique or nonlinear orbit limits")

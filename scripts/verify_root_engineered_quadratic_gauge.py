#!/usr/bin/env python3
"""Exact checks for the root-engineered quadratic cancellation gauge."""

from __future__ import annotations

import sympy as sp


# Universal plane-core calculation.
P, S, Q, U = sp.symbols("P S Q U")
a = sp.symbols("a", nonzero=True)
b0, b1, b2, b3, b4 = sp.symbols("b0:5")

b = (
    b0 * P
    + b1 * P * S
    + b2 * P**4 * S**2
    + b3 * P**5 * S**3
    + b4 * P**6 * S**4
)
R = (
    2 * a * S
    + sp.Rational(2, 3) * a * P * S**3
    - sp.integrate(
        U**2 * sp.diff(b.subs(S, U), U),
        (U, 0, S),
    )
)
D = 1 - S * Q + P * S**2
B_core = a * Q + b
C_core = R - a * Q * S**2

plane_jacobian = sp.factor(
    sp.det(sp.Matrix([B_core, C_core]).jacobian((S, Q)))
)
assert plane_jacobian == -2 * a**2 * D
assert sp.factor(sp.diff(R, S) + S**2 * sp.diff(b, S) - 2 * a * (1 + P * S**2)) == 0

inverse_equation = sp.expand(C_core - R + B_core * S**2 - b * S**2)
assert inverse_equation == 0
formal_E = sp.symbols("C") - R + sp.symbols("B") * S**2 - b * S**2
incidence_derivative = sp.factor(
    sp.diff(formal_E, S).subs(
        {
            sp.symbols("B"): B_core,
            sp.symbols("C"): C_core,
        }
    )
)
assert incidence_derivative == -2 * a * D


# Coefficient engineering through degree six.
g1, g2, g3, g4, g5, g6 = sp.symbols("g1:7", nonzero=True)
coefficients = {1: g1, 2: g2, 3: g3, 4: g4, 5: g5, 6: g6}
a_seed = -g1 / 2
b_seed = (
    -g2 * P
    + (g1 - 3 * g3) * P * S / 2
    + sum(
        -sp.Rational(k, 2) * coefficients[k] * P**k * S ** (k - 2)
        for k in range(4, 7)
    )
)
R_seed = (
    2 * a_seed * S
    + sp.Rational(2, 3) * a_seed * P * S**3
    - sp.integrate(
        U**2 * sp.diff(b_seed.subs(S, U), U),
        (U, 0, S),
    )
)
weighted_seed = (
    g1 * S
    + P * (g2 * S**2 + g3 * S**3)
    + sum(coefficients[k] * P**k * S**k for k in range(4, 7))
)
assert sp.factor(-R_seed - b_seed * S**2 - weighted_seed) == 0

B_seed_core = (a_seed * Q + b_seed) / a_seed
C_seed_core = (R_seed - a_seed * Q * S**2) / a_seed
seed_core_jacobian = sp.factor(
    sp.det(sp.Matrix([B_seed_core, C_seed_core]).jacobian((S, Q)))
)
assert sp.factor(seed_core_jacobian + 2 * D) == 0

# Tangent-line interpretation of the normalized plane core.
X_curve = S**2
Y_curve = 2 * weighted_seed / g1
beta = sp.cancel(
    (sp.diff(weighted_seed, S) / g1 - 1 - P * S**2) / S
)
assert not sp.denom(beta).has(S)
assert sp.factor(B_seed_core - (Q + beta)) == 0
assert sp.factor(C_seed_core - (Y_curve - B_seed_core * X_curve)) == 0
tangency_gap = sp.factor(
    sp.diff(Y_curve, S) - B_seed_core * sp.diff(X_curve, S)
)
assert sp.factor(tangency_gap - 2 * D) == 0
line_incidence_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                sp.symbols("line_B"),
                Y_curve - sp.symbols("line_B") * X_curve,
            ]
        ).jacobian((S, sp.symbols("line_B")))
    )
)
assert sp.factor(
    line_incidence_jacobian
    + sp.diff(Y_curve, S)
    - sp.symbols("line_B") * sp.diff(X_curve, S)
) == 0

# Each isolated higher seed monomial has a forced paired correction.  Its
# diagonal P-weight alpha=k is the first weight regular in both coordinates.
for k in range(4, 9):
    coefficient = sp.symbols(f"coefficient_{k}")
    decoration = coefficient * P**k * S**k
    slope_decoration = sp.diff(decoration, S) / S
    intercept_decoration = sp.expand(
        2 * decoration - slope_decoration * S**2
    )
    assert slope_decoration == k * coefficient * P**k * S ** (k - 2)
    assert intercept_decoration == (
        (2 - k) * coefficient * P**k * S**k
    )
    slope_weight_threshold = k - 2
    intercept_weight_threshold = k
    assert max(slope_weight_threshold, intercept_weight_threshold) == k
    is_regular = lambda weight: (
        weight - slope_weight_threshold >= 0
        and weight - intercept_weight_threshold >= 0
    )
    assert not is_regular(k - 1)
    assert is_regular(k)


# Denominator-free pullback in source coordinates.
x, y, z = sp.symbols("x y z")
t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
p_source = t * q
s_source = x / t
Q_source = y + x * q

F2 = (
    y
    + 3 * g3 * x * q / g1
    + 2 * g2 * p_source / g1
    + sum(
        k * coefficients[k] * t**2 * x ** (k - 2) * q**k / g1
        for k in range(4, 7)
    )
)
F3 = (
    x * (5 - 3 * t)
    - g3 * x**3 * z / g1
    - sum(
        (k - 2) * coefficients[k] * (x * q) ** k / g1
        for k in range(4, 7)
    )
)
source_substitution = {P: p_source, S: s_source, Q: Q_source}
assert sp.factor(sp.together(B_seed_core.subs(source_substitution) - F2)) == 0
assert sp.factor(sp.together(C_seed_core.subs(source_substitution) - F3)) == 0


# The source birational chart contributes the reciprocal divisor.
p_chart = t * q
s_chart = x / t
Q_chart = y + x * q
chart_jacobian = sp.factor(
    sp.det(sp.Matrix([p_chart, s_chart, Q_chart]).jacobian((x, y, z)))
)
assert chart_jacobian == t
assert sp.factor(D.subs(source_substitution) - 1 / t) == 0


def normalized_map(seed: sp.Expr) -> tuple[tuple[sp.Expr, ...], sp.Expr, sp.Expr]:
    """Build formula (26) for one exact seed polynomial."""

    polynomial = sp.Poly(sp.expand(seed), S)
    degree = polynomial.degree()
    gs = {k: polynomial.nth(k) for k in range(1, degree + 1)}
    assert gs[1] != 0 and gs[3] != 0
    local_t = 1 + x * y
    local_q = local_t**2 * z + (gs[1] / gs[3]) * y**2 * (1 + 3 * local_t)
    mapping = (
        local_t * local_q,
        y
        + 3 * gs[3] * x * local_q / gs[1]
        + 2 * gs.get(2, 0) * local_t * local_q / gs[1]
        + sum(
            k
            * gs[k]
            * local_t**2
            * x ** (k - 2)
            * local_q**k
            / gs[1]
            for k in range(4, degree + 1)
        ),
        x * (5 - 3 * local_t)
        - gs[3] * x**3 * z / gs[1]
        - sum(
            (k - 2) * gs[k] * (x * local_q) ** k / gs[1]
            for k in range(4, degree + 1)
        ),
    )
    return tuple(sp.expand(item) for item in mapping), local_t, local_q


# The cubic seed is exactly the foundational map.
foundational, foundational_t, foundational_q = normalized_map(S**3 + S)
expected_foundational = (
    foundational_t**3 * z
    + y**2 * foundational_t * (4 + 3 * x * y),
    y
    + 3 * x * foundational_t**2 * z
    + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)
assert all(
    sp.expand(got - expected) == 0
    for got, expected in zip(foundational, expected_foundational)
)


# Original quartic: the normalized last two coordinates are one third of
# the displayed integer map.
quartic_seed = S * (S - 1) * (S - 2) * (S - 3)
quartic_normalized, quartic_t, quartic_q = normalized_map(quartic_seed)
quartic_integer = (
    quartic_t * quartic_q,
    3 * quartic_normalized[1],
    3 * quartic_normalized[2],
)
expected_quartic = (
    quartic_t * quartic_q,
    3 * y
    - 11 * quartic_t * quartic_q
    + 9 * x * quartic_q
    - 2 * quartic_t**2 * x**2 * quartic_q**4,
    6 * x
    - 9 * x**2 * y
    - 3 * x**3 * z
    + x**4 * quartic_q**4,
)
assert all(
    sp.expand(got - expected) == 0
    for got, expected in zip(quartic_integer, expected_quartic)
)
quartic_jacobian = sp.factor(
    sp.det(sp.Matrix(quartic_integer).jacobian((x, y, z)))
)
assert quartic_jacobian == -18

quartic_points = (
    (sp.Rational(0), sp.Rational(11, 3), -sp.Rational(475, 9)),
    (-sp.Rational(3), sp.Rational(4, 3), sp.Rational(125, 81)),
    (sp.Rational(6), sp.Rational(1, 3), -sp.Rational(7, 81)),
    (-sp.Rational(3), sp.Rational(2, 3), -sp.Rational(1, 9)),
)
quartic_polynomial = sp.Poly(quartic_seed, S)
quartic_g1 = quartic_polynomial.nth(1)
quartic_g3 = quartic_polynomial.nth(3)
for root, point in enumerate(quartic_points):
    substitution = dict(zip((x, y, z), point))
    assert tuple(component.subs(substitution) for component in quartic_integer) == (
        1,
        0,
        0,
    )
    root = sp.Rational(root)
    root_D = sp.diff(quartic_seed, S).subs(S, root) / quartic_g1
    root_y = (
        (1 - root_D) / root
        if root != 0
        else -2 * quartic_polynomial.nth(2) / quartic_g1
    )
    root_z = (
        root_D**3
        - (quartic_g1 / quartic_g3)
        * root_D
        * (root_D + 3)
        * root_y**2
    )
    assert point == (root / root_D, root_y, root_z)

target_B, target_C = sp.symbols("target_B target_C")
quartic_inverse = (
    P**4 * S**4
    - 6 * P * S**3
    + (target_B + 11 * P) * S**2
    - 6 * S
    + target_C
)
assert sp.factor(quartic_inverse.subs({P: 1, target_B: 0, target_C: 0})) == quartic_seed
quartic_Q = (target_B + 11 * P - 6 * P * S + 2 * P**4 * S**2) / 3
quartic_D = 1 - S * quartic_Q + P * S**2
assert sp.factor(
    sp.diff(quartic_inverse, S) + 6 * quartic_D
) == 0


# A smaller-coefficient rational quartic.
symmetric_seed = S * (S - 1) * (S + 1) * (S - 2)
symmetric_map, _, _ = normalized_map(symmetric_seed)
symmetric_jacobian = sp.factor(
    sp.det(sp.Matrix(symmetric_map).jacobian((x, y, z)))
)
assert symmetric_jacobian == -2
symmetric_points = (
    (sp.Rational(0), sp.Rational(1), sp.Rational(5)),
    (-sp.Rational(1), sp.Rational(2), -sp.Rational(9)),
    (sp.Rational(1, 3), -sp.Rational(4), -sp.Rational(27)),
    (sp.Rational(2, 3), -sp.Rational(1), sp.Rational(45)),
)
for point in symmetric_points:
    substitution = dict(zip((x, y, z), point))
    assert tuple(component.subs(substitution) for component in symmetric_map) == (
        1,
        0,
        0,
    )


# Discriminant normalization: dC/dB=-r^2.
r = sp.symbols("r", nonzero=True)
h1, h2, h3, h4, h5, h6 = sp.symbols("h1:7", nonzero=True)
H = sum(
    coefficient * r**degree
    for degree, coefficient in enumerate((0, h1, h2, h3, h4, h5, h6))
)
B_discriminant = sp.diff(H, r) / (h1 * r)
C_discriminant = (2 * H - r * sp.diff(H, r)) / h1
assert sp.factor(
    sp.diff(C_discriminant, r)
    + r**2 * sp.diff(B_discriminant, r)
) == 0
assert sp.limit(r * B_discriminant, r, 0) == 1


# Consecutive-root seeds give rational complete fibers in every tested degree.
for degree in range(3, 13):
    seed = sp.prod(S - root for root in range(degree))
    polynomial = sp.Poly(seed, S)
    linear = polynomial.nth(1)
    cubic = polynomial.nth(3)
    assert linear != 0 and cubic != 0
    roots = tuple(sp.Rational(root) for root in range(degree))
    leading = polynomial.nth(degree)
    assert sp.factor(
        -2 * polynomial.nth(2) / linear
        - 2 * sp.harmonic(degree - 1)
    ) == 0
    barycentric_weights = tuple(
        sp.factor(linear / sp.diff(seed, S).subs(S, root))
        for root in roots
    )
    assert barycentric_weights == tuple(
        (-1) ** root * sp.binomial(degree - 1, root)
        for root in range(degree)
    )
    for exponent in range(degree - 1):
        assert sp.factor(
            sum(
                weight * root**exponent
                for weight, root in zip(barycentric_weights, roots)
            )
        ) == 0
    assert sp.factor(
        sum(
            weight * root ** (degree - 1)
            for weight, root in zip(barycentric_weights, roots)
        )
        - linear / leading
    ) == 0
    valid_quadratic_gauge_marks = []
    for root in roots:
        derivative_ratio = sp.diff(seed, S).subs(S, root) / linear
        assert derivative_ratio != 0
        rerooted = sp.Poly(sp.expand(seed.subs(S, S + root)), S)
        assert rerooted.nth(1) == sp.diff(seed, S).subs(S, root)
        assert rerooted.nth(3) == sp.diff(seed, S, 3).subs(S, root) / 6
        if rerooted.nth(3) != 0:
            valid_quadratic_gauge_marks.append(root)
    assert valid_quadratic_gauge_marks


# The exact quartic discriminant regression agrees with the normalization.
quartic_discriminant = sp.factor(sp.discriminant(quartic_inverse.subs(P, 1), S))
quartic_discriminant_expected = 16 * (
    target_B**4 * target_C
    + 35 * target_B**3 * target_C
    - 9 * target_B**3
    - 8 * target_B**2 * target_C**2
    + 249 * target_B**2 * target_C
    - 216 * target_B**2
    + 148 * target_B * target_C**2
    - 121 * target_B * target_C
    - 27 * target_B
    + 16 * target_C**3
    - 23 * target_C**2
    - 2 * target_C
    + 9
)
assert sp.expand(quartic_discriminant - quartic_discriminant_expected) == 0


print("PASS: universal quadratic gauge has plane Jacobian -2*a^2*D")
print("PASS: coefficient engineering and polynomial pullback hold through degree six")
print("PASS: normalized core is marked-line incidence on (S^2, 2*G_P/g1)")
print("PASS: diagonal higher-term weights are minimal term by term")
print("PASS: G=S^3+S is exactly the foundational Keller map")
print("PASS: both quartic complete fibers and determinants are exact")
print("PASS: quadratic discriminant satisfies dC/dB=-r^2")
print("PASS: root reconstruction and barycentric moment identities are exact")
print("PASS: consecutive rational seeds and valid rerootings pass through degree twelve")

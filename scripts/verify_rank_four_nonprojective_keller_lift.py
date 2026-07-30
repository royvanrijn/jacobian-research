#!/usr/bin/env python3
"""Exact rank-four nonprojective Keller-lift regression.

The checker separates three statements:

1. the rank-four quadratic-gauge orbit has a fifth-power arithmetic twist;
2. an arithmetic-neutral primitive quadratic change reduces to two marked
   fibers of one fixed quartic Keller map over QQ;
3. the resulting marked-fiber motion has an exact polynomial first-order
   lift and a formal lift by the repository's formal-orbit theorem.

It does not assert that the formal lift algebraizes at the endpoint.
"""

from __future__ import annotations

import sympy as sp


S, x, y, z = sp.symbols("S x y z")


def assert_zero(expression: sp.Expr) -> None:
    """Assert an exact rational identity."""

    assert sp.cancel(expression) == 0


def normalized_relation(
    roots: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, ...]:
    """Return coefficients of the quartic relation normalized at S."""

    polynomial = sp.Poly(
        sp.expand(sp.prod(S - root for root in roots)),
        S,
        domain=sp.QQ.frac_field(*sorted(
            set().union(*(root.free_symbols for root in roots)),
            key=lambda symbol: symbol.name,
        )),
    )
    linear = polynomial.coeff_monomial(S)
    assert linear != 0
    normalized = sp.Poly(sp.cancel(polynomial.as_expr() / linear), S)
    return tuple(
        sp.factor(normalized.coeff_monomial(S**degree))
        for degree in range(5)
    )


def seed_and_target(
    coefficients: tuple[sp.Expr, ...],
) -> tuple[sp.Expr, tuple[sp.Expr, sp.Expr, sp.Expr]]:
    """Compress a normalized quartic relation into (U; pi,b,c)."""

    a0, a1, a2, a3, a4 = coefficients
    assert a1 == 1
    seed = sp.factor(a4 / a3**4)
    return seed, (a3, a2, -2 * a0)


def compressed_quartic_map(seed: sp.Expr) -> tuple[sp.Expr, ...]:
    """Jacobian-one quartic quadratic-gauge map F_seed."""

    local_t = 1 + x * y
    local_q = local_t**2 * z + y**2 * (1 + 3 * local_t)
    pi = local_t * local_q
    raw_b = (
        y
        + 3 * x * local_q
        + 4 * seed * local_t**2 * x**2 * local_q**4
    )
    c = (
        x * (5 - 3 * local_t)
        - x**3 * z
        - 2 * seed * (x * local_q) ** 4
    )
    return tuple(sp.expand(component) for component in (pi, -raw_b / 2, c))


def inverse_relation(
    seed: sp.Expr,
    target: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> sp.Expr:
    """Normalized inverse polynomial at a target (pi,b,c)."""

    pi, b, c = target
    return sp.expand(seed * pi**4 * S**4 + pi * S**3 + b * S**2 + S - c / 2)


def rational_valuations(value: sp.Rational) -> dict[int, int]:
    """Prime valuations of a nonzero rational number."""

    value = sp.Rational(value)
    valuations = dict(sp.factorint(abs(int(value.p))))
    for prime, exponent in sp.factorint(abs(int(value.q))).items():
        valuations[prime] = valuations.get(prime, 0) - exponent
    return {prime: exponent for prime, exponent in valuations.items() if exponent}


# ---------------------------------------------------------------------------
# 1. Arithmetic form of the rank-four two-torus orbit.
# ---------------------------------------------------------------------------

alpha, beta = sp.symbols("alpha beta", nonzero=True)
a3, a4, b3, b4 = sp.symbols("a3 a4 b3 b4", nonzero=True)

beta_from_cubic = a3 / (b3 * alpha**2)
transformed_quartic = a4 * alpha**-3 * beta_from_cubic**-4
kummer_ratio = (b4 / b3**4) / (a4 / a3**4)
assert_zero(transformed_quartic / b4 - alpha**5 / kummer_ratio)

# Thus a K-rational two-torus equivalence exists exactly when this ratio is
# a fifth power in K.  Necessity of the two-torus form is imported from RQG3.


# ---------------------------------------------------------------------------
# 2. The smallest witness has a genuine QQ Kummer twist.
# ---------------------------------------------------------------------------

original_roots = tuple(map(sp.Integer, (1, 2, 3, 4)))
original_changed_roots = tuple(root + root**2 for root in original_roots)
original_coefficients = normalized_relation(original_roots)
original_changed_coefficients = normalized_relation(original_changed_roots)

assert original_coefficients == (
    sp.Rational(-12, 25),
    1,
    sp.Rational(-7, 10),
    sp.Rational(1, 5),
    sp.Rational(-1, 50),
)
assert original_changed_coefficients == (
    sp.Rational(-5, 4),
    1,
    sp.Rational(-127, 576),
    sp.Rational(5, 288),
    sp.Rational(-1, 2304),
)

original_seed, original_target = seed_and_target(original_coefficients)
original_changed_seed, original_changed_target = seed_and_target(
    original_changed_coefficients
)
assert original_seed == sp.Rational(-25, 2)
assert original_target == (
    sp.Rational(1, 5),
    sp.Rational(-7, 10),
    sp.Rational(24, 25),
)
assert original_changed_seed == sp.Rational(-2985984, 625)
assert original_changed_target == (
    sp.Rational(5, 288),
    sp.Rational(-127, 576),
    sp.Rational(5, 2),
)

original_seed_ratio = sp.factor(original_changed_seed / original_seed)
assert original_seed_ratio == sp.Rational(5971968, 15625)
assert rational_valuations(original_seed_ratio) == {2: 13, 3: 6, 5: -6}
assert any(
    exponent % 5
    for exponent in rational_valuations(original_seed_ratio).values()
)

original_projective_matrix = sp.Matrix(
    [
        [1, root, root + root**2, root * (root + root**2)]
        for root in original_roots
    ]
)
assert original_projective_matrix.det() != 0


# ---------------------------------------------------------------------------
# 3. An arithmetic-neutral primitive nonprojective witness.
# ---------------------------------------------------------------------------

neutral_roots = tuple(map(sp.Integer, (-6, -3, 0, 3)))
neutral_changed_roots = tuple(root + root**2 - 18 for root in neutral_roots)
assert neutral_changed_roots == (12, -12, -18, -6)
assert len(set(neutral_changed_roots)) == 4

neutral_coefficients = normalized_relation(neutral_roots)
neutral_changed_coefficients = normalized_relation(neutral_changed_roots)
assert neutral_coefficients == (
    0,
    1,
    sp.Rational(1, 6),
    sp.Rational(-1, 9),
    sp.Rational(-1, 54),
)
assert neutral_changed_coefficients == (
    sp.Rational(9, 2),
    1,
    sp.Rational(1, 96),
    sp.Rational(-1, 144),
    sp.Rational(-1, 3456),
)

seed_0, target_0 = seed_and_target(neutral_coefficients)
seed_1, target_1 = seed_and_target(neutral_changed_coefficients)
assert seed_0 == sp.Rational(-243, 2)
assert seed_1 == -124416
assert target_0 == (sp.Rational(-1, 9), sp.Rational(1, 6), 0)
assert target_1 == (sp.Rational(-1, 144), sp.Rational(1, 96), -9)
assert sp.factor(seed_1 / seed_0) == 4**5

neutral_projective_matrix = sp.Matrix(
    [
        [1, root, root + root**2 - 18, root * (root + root**2 - 18)]
        for root in neutral_roots
    ]
)
assert neutral_projective_matrix.rank() == 4


# ---------------------------------------------------------------------------
# 4. Exact QQ left-right normalization to one fixed Keller map.
# ---------------------------------------------------------------------------

symbolic_seed = sp.symbols("symbolic_seed")
F_symbolic = compressed_quartic_map(symbolic_seed)
jacobian_symbolic = sp.Matrix(F_symbolic).jacobian((x, y, z))
assert_zero(jacobian_symbolic.det() - 1)

source_scaling = {
    x: alpha * x,
    y: y / alpha,
    z: z / alpha**2,
}
F_scaled_seed = compressed_quartic_map(alpha**5 * symbolic_seed)
scaled_source_image = tuple(
    sp.factor(component.subs(source_scaling, simultaneous=True))
    for component in F_scaled_seed
)
scaled_target_image = (
    F_symbolic[0] / alpha**2,
    F_symbolic[1] / alpha,
    alpha * F_symbolic[2],
)
assert all(
    sp.factor(left - right) == 0
    for left, right in zip(scaled_source_image, scaled_target_image)
)

linearized_target_0 = (
    target_0[0] / 16,
    target_0[1] / 4,
    4 * target_0[2],
)
assert linearized_target_0 == (
    sp.Rational(-1, 144),
    sp.Rational(1, 24),
    0,
)
target_residual = tuple(
    sp.factor(right - left)
    for left, right in zip(linearized_target_0, target_1)
)
assert target_residual == (0, sp.Rational(-1, 32), -9)

scaled_old_roots = tuple(4 * root for root in neutral_roots)
assert normalized_relation(scaled_old_roots) == tuple(
    sp.Poly(inverse_relation(seed_1, linearized_target_0), S).coeff_monomial(
        S**degree
    )
    for degree in range(5)
)
assert normalized_relation(neutral_changed_roots) == tuple(
    sp.Poly(inverse_relation(seed_1, target_1), S).coeff_monomial(S**degree)
    for degree in range(5)
)
assert tuple(
    scaled_root**2 / 16 + scaled_root / 4 - 18
    for scaled_root in scaled_old_roots
) == neutral_changed_roots


# ---------------------------------------------------------------------------
# 5. A rational path through primitive nonprojective labels.
# ---------------------------------------------------------------------------

path_parameter = sp.symbols("path_parameter")
path_roots = tuple(
    root + path_parameter * (root**2 - 18)
    for root in neutral_roots
)
raw_path_relation = sp.Poly(
    sp.expand(sp.prod(S - root for root in path_roots)),
    S,
)
path_linear_coefficient = sp.factor(
    raw_path_relation.coeff_monomial(S)
)
assert_zero(
    path_linear_coefficient
    + 54
    * (3 * path_parameter - 1)
    * (36 * path_parameter**2 - 3 * path_parameter - 1)
)

path_coefficients = normalized_relation(path_roots)
expected_path_coefficients = (
    18
    * path_parameter
    * (3 * path_parameter - 1)
    * (3 * path_parameter + 1)
    / (36 * path_parameter**2 - 3 * path_parameter - 1),
    1,
    (27 * path_parameter**2 - 24 * path_parameter + 1)
    / (
        6
        * (3 * path_parameter - 1)
        * (36 * path_parameter**2 - 3 * path_parameter - 1)
    ),
    -(3 * path_parameter + 1)
    / (
        9
        * (3 * path_parameter - 1)
        * (36 * path_parameter**2 - 3 * path_parameter - 1)
    ),
    -1
    / (
        54
        * (3 * path_parameter - 1)
        * (36 * path_parameter**2 - 3 * path_parameter - 1)
    ),
)
assert all(
    sp.factor(got - expected) == 0
    for got, expected in zip(path_coefficients, expected_path_coefficients)
)

path_seed, path_target = seed_and_target(path_coefficients)
expected_path_seed = (
    -243
    * (3 * path_parameter - 1) ** 3
    * (36 * path_parameter**2 - 3 * path_parameter - 1) ** 3
    / (2 * (3 * path_parameter + 1) ** 4)
)
assert_zero(path_seed - expected_path_seed)
assert sp.factor(path_seed.subs(path_parameter, 0) - seed_0) == 0
assert sp.factor(path_seed.subs(path_parameter, 1) - seed_1) == 0
assert tuple(
    sp.factor(coordinate.subs(path_parameter, 0))
    for coordinate in path_target
) == target_0
assert tuple(
    sp.factor(coordinate.subs(path_parameter, 1))
    for coordinate in path_target
) == target_1

path_vandermonde = sp.factor(
    sp.prod(
        path_roots[second] - path_roots[first]
        for first in range(4)
        for second in range(first + 1, 4)
    )
)
assert_zero(
    path_vandermonde
    - 8748
    * (3 * path_parameter - 1) ** 2
    * (3 * path_parameter + 1)
    * (6 * path_parameter - 1)
    * (9 * path_parameter - 1)
)
assert path_vandermonde.subs(path_parameter, 0) != 0
assert path_vandermonde.subs(path_parameter, 1) != 0

seed_derivative = sp.factor(sp.diff(path_seed, path_parameter).subs(path_parameter, 0))
target_derivative = tuple(
    sp.factor(sp.diff(coordinate, path_parameter).subs(path_parameter, 0))
    for coordinate in path_target
)
assert seed_derivative == 1458
assert target_derivative == (
    sp.Rational(-1, 3),
    -4,
    -36,
)


# ---------------------------------------------------------------------------
# 6. Exact source reconstruction and first-order marked-fiber transport.
# ---------------------------------------------------------------------------

F_0 = compressed_quartic_map(seed_0)
jacobian_0 = sp.Matrix(F_0).jacobian((x, y, z))
assert_zero(jacobian_0.det() - 1)


def reconstruct_source(
    root: sp.Expr,
    seed: sp.Expr,
    target: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Reconstruct one source point from an inverse root."""

    pi, b, _ = target
    beta_polynomial = 2 * pi * root + 4 * seed * pi**4 * root**2
    q_coordinate = -2 * b - beta_polynomial
    derivative = 1 - root * q_coordinate + pi * root**2
    local_t = 1 / derivative
    source_x = root * local_t
    source_y = q_coordinate - pi * root
    source_q = pi * derivative
    source_z = derivative**2 * (
        source_q - source_y**2 * (1 + 3 * local_t)
    )
    return tuple(
        sp.factor(coordinate)
        for coordinate in (source_x, source_y, source_z)
    )


source_points = tuple(
    reconstruct_source(root, seed_0, target_0)
    for root in neutral_roots
)
assert source_points == (
    (sp.Integer(-2), sp.Rational(1, 3), sp.Integer(-5)),
    (sp.Integer(3), sp.Rational(-2, 3), sp.Integer(1)),
    (sp.Integer(0), sp.Rational(-1, 3), sp.Rational(-5, 9)),
    (sp.Integer(-1), sp.Rational(4, 3), sp.Integer(3)),
)

for root, point in zip(neutral_roots, source_points):
    assert tuple(
        sp.factor(component.subs(dict(zip((x, y, z), point))))
        for component in F_0
    ) == target_0
    recovered_root = sp.factor(
        (x / (1 + x * y)).subs(dict(zip((x, y, z), point)))
    )
    assert recovered_root == root

parameter_derivative = sp.Matrix(
    [
        sp.diff(component, symbolic_seed).subs(symbolic_seed, seed_0)
        for component in F_symbolic
    ]
)
moving_map_derivative = seed_derivative * parameter_derivative
moving_target_derivative = sp.Matrix(target_derivative)

# If p_t is transported inside F_{U(t)}^{-1}(y_t), then p'_0=X(p_0).
transport_field = jacobian_0.adjugate() * (
    moving_target_derivative - moving_map_derivative
)
transport_field = sp.Matrix(
    [sp.expand(component) for component in transport_field]
)
assert all(
    sp.expand(component) == 0
    for component in (
        jacobian_0 * transport_field
        + moving_map_derivative
        - moving_target_derivative
    )
)
assert_zero(
    sum(
        sp.diff(transport_field[index], variable)
        for index, variable in enumerate((x, y, z))
    )
)

transport_polynomials = [
    sp.Poly(component, x, y, z, domain=sp.QQ)
    for component in transport_field
]
assert tuple(polynomial.total_degree() for polynomial in transport_polynomials) == (
    55,
    53,
    55,
)
assert tuple(len(polynomial.terms()) for polynomial in transport_polynomials) == (
    512,
    510,
    612,
)

expected_point_velocities = (
    (-30, -8, 303),
    (45, 6, -25),
    (-18, 2, sp.Rational(31, 3)),
    (3, 4, 9),
)
root_coordinate = x / (1 + x * y)
root_gradient = sp.Matrix(
    [sp.diff(root_coordinate, variable) for variable in (x, y, z)]
)
for root, point, expected_velocity in zip(
    neutral_roots,
    source_points,
    expected_point_velocities,
):
    substitution = dict(zip((x, y, z), point))
    point_velocity = tuple(
        sp.factor(component.subs(substitution))
        for component in transport_field
    )
    assert point_velocity == expected_velocity
    induced_root_velocity = sp.factor(
        root_gradient.dot(transport_field).subs(substitution)
    )
    assert induced_root_velocity == root**2 - 18


# ---------------------------------------------------------------------------
# 7. The straight fixed-map target line is formally liftable but frame-wrong.
# ---------------------------------------------------------------------------

line_parameter = sp.symbols("line_parameter")
line_target = tuple(
    start + line_parameter * delta
    for start, delta in zip(linearized_target_0, target_residual)
)
line_relation = sp.factor(inverse_relation(seed_1, line_target))
assert_zero(
    line_relation
    + (S - 12)
    * (S + 12)
    * (S**2 + 24 * S + 108 * line_parameter)
    / 3456
)
line_discriminant = sp.factor(sp.discriminant(line_relation, S))
assert_zero(
    line_discriminant
    + (line_parameter + 4) ** 2
    * (3 * line_parameter - 4) ** 3
    / sp.Integer(1358954496)
)
assert line_discriminant.subs(line_parameter, 0) != 0
assert line_discriminant.subs(line_parameter, 1) != 0

F_1 = compressed_quartic_map(seed_1)
jacobian_1 = sp.Matrix(F_1).jacobian((x, y, z))
translation_field = jacobian_1.adjugate() * sp.Matrix(target_residual)
translation_field = sp.Matrix(
    [sp.expand(component) for component in translation_field]
)
assert all(
    sp.expand(component) == 0
    for component in (
        jacobian_1 * translation_field - sp.Matrix(target_residual)
    )
)
assert_zero(
    sum(
        sp.diff(translation_field[index], variable)
        for index, variable in enumerate((x, y, z))
    )
)
translation_polynomials = [
    sp.Poly(component, x, y, z, domain=sp.QQ)
    for component in translation_field
]
assert tuple(
    polynomial.total_degree() for polynomial in translation_polynomials
) == (31, 29, 31)
assert tuple(
    len(polynomial.terms()) for polynomial in translation_polynomials
) == (118, 115, 149)

scaled_source_points = tuple(
    (4 * point[0], point[1] / 4, point[2] / 16)
    for point in source_points
)
line_root_velocities = []
for point in scaled_source_points:
    substitution = dict(zip((x, y, z), point))
    line_root_velocities.append(
        sp.factor(root_gradient.dot(translation_field).subs(substitution))
    )
assert tuple(line_root_velocities) == (
    sp.Rational(9, 2),
    0,
    sp.Rational(-9, 2),
    0,
)

# The line family has constant sections S=+/-12 and one quadratic sheet.
# The desired quadratic Tschirnhaus labels cross that decomposition, so an
# algebraization of this based straight-line formal lift would not by itself
# realize the required collision frame.
constant_start = {-12, 12}
constant_end = {-12, 12}
desired_images_of_constant_start = {
    root**2 / 16 + root / 4 - 18
    for root in constant_start
}
assert desired_images_of_constant_start == {-12, -6}
assert desired_images_of_constant_start != constant_end


print("PASS: the quartic orbit twist is the exact fifth-power Kummer class")
print("PASS: the 1,2,3,4 witness is nonprojective and nontrivial over QQ")
print("PASS: the neutral witness has seed ratio 4^5 and one fixed QQ map")
print("PASS: the endpoint fibers differ by target residual (0,-1/32,-9)")
print("PASS: the rational label path is primitive at both endpoints")
print("PASS: det(DF_U)=1 and all four source points reconstruct exactly")
print("PASS: the marked motion has a divergence-free polynomial lift")
print("PASS: lift degrees are (55,53,55), with terms (512,510,612)")
print("PASS: the fixed-map target line is finite etale at both endpoints")
print("PASS: its formal translation lift has degrees (31,29,31)")
print("PASS: the straight-line sheet partition does not realize the q-frame")
print("SCOPE: no endpoint algebraization or global self-equivalence is claimed")

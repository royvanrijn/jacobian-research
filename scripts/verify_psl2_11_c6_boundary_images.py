#!/usr/bin/env python3
"""Replay the exact C6 boundary traces through its elliptic quotients.

This is the expensive residue-field companion to
verify_psl2_11_keller_action_spectrum.py.  Singular separates the normalized
boundary primes and computes their canonical (t,w)-images; SymPy performs
Riemann--Roch interpolation on the two elliptic quotients and identifies the
resulting K-rational traces in the pinned Mordell--Weil basis.
"""

import re
import subprocess

import sympy as sp


K = sp.QQ.algebraic_field(sp.sqrt(-11))
r = K.ext
X, Y, T, W = sp.symbols("X Y T W")


def parse_singular_polynomial(source):
    source = source.rstrip(",")
    source = re.sub(r"(?<=\d)r", "*r", source)
    source = re.sub(r"([XYtw])(\d+)", r"\1**\2", source)
    return sp.sympify(
        source,
        locals={"X": X, "Y": Y, "t": T, "w": W, "r": r},
    )


def quotient_reduce_rational(expression, modulus, variable=W):
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    modulus_poly = sp.Poly(modulus, variable, domain=K)
    numerator_poly = sp.Poly(numerator, variable, domain=K)
    denominator_poly = sp.Poly(denominator, variable, domain=K)
    inverse_denominator = sp.invert(denominator_poly, modulus_poly)
    return sp.rem(numerator_poly * inverse_denominator, modulus_poly).as_expr()


def elliptic_map_expressions(sign):
    znum = 8 * T - 3 + r
    zden = 8 * T + 5 + r
    twist_d = 2 * (1 + 3 * r)
    if sign == "PLUS":
        scale = -22 - 6 * r
        shift = 110 + 330 * r
        ordinate_shift = 7744 - 3168 * r
        x_expression = (
            -121 * twist_d * znum**2 - shift * zden**2
        ) / (zden**2 * scale**2)
        y_expression = (
            3872 * twist_d * W - ordinate_shift * zden**3
        ) / (zden**3 * scale**3)
    else:
        scale = 6 - 2 * r
        shift = 2 + 6 * r
        ordinate_shift = -288 - 64 * r
        x_expression = (
            twist_d * zden**2 - shift * znum**2
        ) / (znum**2 * scale**2)
        y_expression = (
            -32 * twist_d * W - ordinate_shift * znum**3
        ) / (znum**3 * scale**3)
    return x_expression, y_expression


def image_ideal_from_residue(modulus, variable, x_expression, y_expression):
    x_value = quotient_reduce_rational(x_expression, modulus, variable)
    y_value = quotient_reduce_rational(y_expression, modulus, variable)
    coefficient_domain = K.poly_ring(Y)
    modulus_coefficients = sp.Poly(
        modulus,
        variable,
        domain=K,
    ).rep.to_list()
    modulus_for_resultant = sp.Poly.from_list(
        [coefficient_domain.convert(value, K) for value in modulus_coefficients],
        gens=variable,
        domain=coefficient_domain,
    )
    ordinate_coefficients = [
        -coefficient_domain.convert(value, K)
        for value in sp.Poly(y_value, variable, domain=K).rep.to_list()
    ]
    if not ordinate_coefficients:
        ordinate_coefficients = [coefficient_domain.zero]
    ordinate_coefficients[-1] += coefficient_domain.gens[0]
    ordinate_for_resultant = sp.Poly.from_list(
        ordinate_coefficients,
        gens=variable,
        domain=coefficient_domain,
    )
    f_y_raw = sp.Poly(
        modulus_for_resultant.resultant(ordinate_for_resultant),
        Y,
        domain=K,
    )
    f_y = sp.sqf_part(f_y_raw).monic()
    image_degree = f_y.degree()
    modulus_poly = sp.Poly(modulus, variable, domain=K)
    powers = []
    current = sp.Poly(1, variable, domain=K)
    for _ in range(image_degree):
        powers.append(current)
        current = sp.rem(
            current * sp.Poly(y_value, variable, domain=K),
            modulus_poly,
        )
    x_value_poly = sp.Poly(x_value, variable, domain=K)
    source_degree = modulus_poly.degree()
    matrix = sp.Matrix(
        [
            [power.nth(exponent) for power in powers]
            for exponent in range(source_degree)
        ]
    )
    target = sp.Matrix(
        [x_value_poly.nth(exponent) for exponent in range(source_degree)]
    )
    coefficient_tuple = next(iter(sp.linsolve((matrix, target))))
    x_of_y = sum(
        K.to_sympy(K.from_sympy(sp.cancel(coefficient))) * Y**exponent
        for exponent, coefficient in enumerate(coefficient_tuple)
    )
    return f_y.as_expr(), X - x_of_y


def image_ideal_from_canonical(polynomials, sign):
    f_w = next(polynomial for polynomial in polynomials if not polynomial.has(T))
    t_relation = next(polynomial for polynomial in polynomials if polynomial.has(T))
    t_coefficient = sp.Poly(t_relation, T, domain=K.poly_ring(W)).coeff_monomial(T)
    t_expression = sp.cancel(-(t_relation - t_coefficient * T) / t_coefficient)
    x_expression, y_expression = elliptic_map_expressions(sign)
    return image_ideal_from_residue(
        f_w,
        W,
        x_expression.subs(T, t_expression),
        y_expression.subs(T, t_expression),
    )


def reduce_y(expression, f_y):
    return sp.rem(
        sp.Poly(sp.expand(expression), Y, domain=K),
        sp.Poly(f_y, Y, domain=K),
    ).as_expr()


def closed_point_sum(polynomials, coefficients):
    f_y = next(p for p in polynomials if not p.has(X))
    x_relation = next(p for p in polynomials if p.has(X))
    x_coefficient = sp.Poly(x_relation, X, domain=K.poly_ring(Y)).coeff_monomial(X)
    x_expression = sp.cancel(-(x_relation - x_coefficient * X) / x_coefficient)
    degree = sp.degree(f_y, Y)

    for pole_bound in (degree, degree + 1):
        basis = []
        for y_degree in (0, 1):
            for x_degree in range(pole_bound // 2 + 1):
                if 2 * x_degree + 3 * y_degree <= pole_bound:
                    basis.append(X**x_degree * Y**y_degree)
        reductions = [
            sp.Poly(
                reduce_y(term.subs(X, x_expression), f_y),
                Y,
                domain=K,
            )
            for term in basis
        ]
        matrix = sp.Matrix(
            [
                [reduction.nth(power) for reduction in reductions]
                for power in range(degree)
            ]
        )
        kernel = matrix.nullspace()
        if kernel:
            kernel_vector = [
                K.to_sympy(K.from_sympy(sp.cancel(coefficient)))
                for coefficient in kernel[0]
            ]
            function = sp.expand(
                sum(coefficient * term for coefficient, term in zip(kernel_vector, basis))
            )
            break
    else:
        raise AssertionError("Riemann--Roch interpolation failed")

    if pole_bound == degree:
        return degree, None, function

    curve_a1, curve_a2, curve_a3, curve_a4, curve_a6 = coefficients
    assert (curve_a1, curve_a3) == (0, 1)
    curve_rhs = X**3 + curve_a2 * X**2 + curve_a4 * X + curve_a6
    function_poly = sp.Poly(function, Y)
    function_b = function_poly.coeff_monomial(Y)
    function_a = function_poly.coeff_monomial(1)
    norm = sp.Poly(
        sp.expand(function_a * (function_a - function_b) - function_b**2 * curve_rhs),
        X,
        domain=K,
    )
    q_x = sp.Poly(sp.resultant(f_y, x_relation, Y), X, domain=K)
    quotient, remainder = sp.div(norm, q_x)
    if not remainder.is_zero:
        # Resultants and norms can carry powers/constant scalings; remove gcd.
        quotient, remainder = sp.div(norm, sp.gcd(norm, q_x))
    assert remainder.is_zero
    assert quotient.degree() == 1
    residual_x = -quotient.nth(0) / quotient.nth(1)
    a_at_residual = sp.Poly(function_a, X, domain=K).eval(residual_x)
    b_at_residual = sp.Poly(function_b, X, domain=K).eval(residual_x)
    if b_at_residual == 0:
        assert degree == 1
        f_y_poly = sp.Poly(f_y, Y, domain=K)
        point_sum = (
            sp.cancel(residual_x),
            sp.cancel(-f_y_poly.nth(0) / f_y_poly.nth(1)),
        )
        return degree, point_sum, function
    residual_y = -a_at_residual / b_at_residual
    # The closed point sum is -R.
    point_sum = (sp.cancel(residual_x), sp.cancel(-residual_y - 1))
    lhs = point_sum[1] ** 2 + point_sum[1]
    rhs = (
        point_sum[0] ** 3
        + curve_a2 * point_sum[0] ** 2
        + curve_a4 * point_sum[0]
        + curve_a6
    )
    if K.from_sympy(sp.expand(lhs - rhs)) != K.zero:
        print("DEBUG", f_y, x_relation, function, norm, q_x, quotient, point_sum, lhs-rhs)
        raise AssertionError("residual point is off the target curve")
    return degree, point_sum, function


def inverse(point):
    if point is None:
        return None
    return point[0], -point[1] - K.one


def to_field_point(point):
    if point is None:
        return None
    return tuple(K.from_sympy(sp.cancel(value)) for value in point)


def add(left, right, coefficients):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    _, a2_raw, _, a4_raw, a6_raw = coefficients
    a2, a4, a6 = map(K.convert, (a2_raw, a4_raw, a6_raw))
    if x1 == x2:
        if y1 + y2 + K.one == K.zero:
            return None
        denominator = 2 * y1 + K.one
        slope = (3 * x1**2 + 2 * a2 * x1 + a4) / denominator
        intercept = (-x1**3 + a4 * x1 + 2 * a6 - y1) / denominator
    else:
        slope = (y2 - y1) / (x2 - x1)
        intercept = (y1 * x2 - y2 * x1) / (x2 - x1)
    x3 = slope**2 - a2 - x1 - x2
    y3 = -slope * x3 - intercept - K.one
    return x3, y3


def multiply(scalar, point, coefficients):
    if scalar < 0:
        return multiply(-scalar, inverse(point), coefficients)
    answer = None
    summand = point
    while scalar:
        if scalar & 1:
            answer = add(answer, summand, coefficients)
        scalar >>= 1
        if scalar:
            summand = add(summand, summand, coefficients)
    return answer


def normalize(point):
    if point is None:
        return "O"
    return "(" + ", ".join(str(K.to_sympy(value)) for value in point) + ")"


def points_equal(left, right):
    if left is None or right is None:
        return left is None and right is None
    return left == right


result = subprocess.run(
    ["Singular", "-q", "scripts/psl2_11_c6_boundary_images.sing"],
    text=True,
    capture_output=True,
    check=True,
)
lines = result.stdout.splitlines()
images = {}
for index, line in enumerate(lines):
    if line.startswith(("PLUS ", "MINUS ")):
        label = line.replace(" ", ":", 1)
        images[label] = tuple(
            parse_singular_polynomial(lines[index + offset])
            for offset in (1, 2)
        )

node_result = subprocess.run(
    ["Singular", "-q", "scripts/psl2_11_c6_node_tangent.sing"],
    text=True,
    capture_output=True,
    check=True,
)
node_lines = node_result.stdout.splitlines()
for index, line in enumerate(node_lines):
    if line.startswith("CANONICAL "):
        label = line.split()[1]
        canonical_polynomials = tuple(
            parse_singular_polynomial(node_lines[index + offset])
            for offset in (1, 2)
        )
        for sign in ("PLUS", "MINUS"):
            images[f"{sign}:{label}"] = image_ideal_from_canonical(
                canonical_polynomials,
                sign,
            )

canonical_qa = (
    (128 * r + 320) * T**6
    + (64 * r - 736) * T**5
    + (-110 * r - 50) * T**4
    + (-10 * r + 138) * T**3
    + (sp.Rational(79, 8) * r - sp.Rational(181, 8)) * T**2
    + (-sp.Rational(15, 16) * r + sp.Rational(69, 16)) * T
    - sp.Rational(1, 64) * r
    - sp.Rational(9, 32)
)
canonical_qb = (
    (544 * r - 1760) * T**6
    + (-2320 * r - 3344) * T**5
    + (128 * r + 6908) * T**4
    + (508 * r - 1056) * T**3
    + (-sp.Rational(781, 8) * r - sp.Rational(1859, 8)) * T**2
    + (sp.Rational(87, 16) * r + sp.Rational(1353, 16)) * T
    - sp.Rational(97, 128) * r
    - sp.Rational(957, 128)
)
canonical_square = (
    T**3
    + (r / 40 - sp.Rational(7, 40)) * T**2
    + (7 * r / 160 - sp.Rational(19, 160)) * T
    - r / 320
    + sp.Rational(1, 160)
)
_, infinity_factors = sp.factor_list(canonical_qa, T, extension=r)
assert [sp.degree(factor, T) for factor, _ in infinity_factors] == [1, 5]
infinity_w = -2 * canonical_qb / (3 * r * canonical_square)
for infinity_index, (infinity_factor, exponent) in enumerate(
    infinity_factors,
    start=1,
):
    assert exponent == 1
    for sign in ("PLUS", "MINUS"):
        x_expression, y_expression = elliptic_map_expressions(sign)
        images[f"{sign}:infinity_{infinity_index}"] = image_ideal_from_residue(
            infinity_factor,
            T,
            x_expression.subs(W, infinity_w),
            y_expression.subs(W, infinity_w),
        )

plus_coefficients = (0, 2, 1, 1, 0)
minus_coefficients = (0, 2, 1, -6, 3)
source_degrees = {
    "p1_a": 2,
    "p1_b": 2,
    "p2_a": 2,
    "p2_b": 2,
    "p3_a": 2,
    "p3_b": 2,
    "q1_a": 3,
    "q1_b": 3,
    "q2_a": 3,
    "q2_b": 3,
    "q2_c": 6,
    "q3_a": 3,
    "q3_b": 3,
    "q3_c": 6,
    "p1_node_1": 2,
    "p2_node_a": 4,
    "p2_node_b": 4,
    "q2_node": 6,
    "infinity_1": 1,
    "infinity_2": 5,
}

pushforwards = {}
image_degrees = {}
for key, polynomials in images.items():
    sign, label = key.split(":")
    coefficients = plus_coefficients if sign == "PLUS" else minus_coefficients
    degree, point, _ = closed_point_sum(polynomials, coefficients)
    multiplicity = source_degrees[label] // degree
    pushforward = multiply(multiplicity, to_field_point(point), coefficients)
    pushforwards[key] = pushforward
    image_degrees[key] = degree

direct_labels = (
    "p1_node_1",
    "p1_a",
    "p1_b",
    "p2_a",
    "p2_node_a",
    "p2_node_b",
    "p2_b",
    "p3_a",
    "p3_b",
    "q1_a",
    "q1_b",
    "q2_a",
    "q2_b",
    "q2_c",
    "q3_a",
    "q3_b",
    "q3_c",
    "infinity_1",
    "infinity_2",
)
assert set(images) == {
    f"{sign}:{label}"
    for sign in ("PLUS", "MINUS")
    for label in direct_labels
}

# Expected traces in the basis T on E_+ and (G0,G2) on E_-.
plus_classes = {
    "p1_node_1": 0,
    "p1_a": 1,
    "p1_b": 0,
    "p2_a": 1,
    "p2_node_a": 4,
    "p2_node_b": 3,
    "p2_b": 4,
    "p3_a": 0,
    "p3_b": 4,
    "q1_a": 1,
    "q1_b": 0,
    "q2_a": 4,
    "q2_b": 1,
    "q2_c": 2,
    "q3_a": 4,
    "q3_b": 0,
    "q3_c": 0,
    "infinity_1": 2,
    "infinity_2": 4,
}
minus_classes = {
    "p1_node_1": (4, 0),
    "p1_a": (-1, 1),
    "p1_b": (3, -1),
    "p2_a": (-1, -1),
    "p2_node_a": (0, 0),
    "p2_node_b": (12, 0),
    "p2_b": (1, 1),
    "p3_a": (3, 1),
    "p3_b": (1, -1),
    "q1_a": (2, -1),
    "q1_b": (4, 1),
    "q2_a": (0, -1),
    "q2_b": (2, 1),
    "q2_c": (14, 0),
    "q3_a": (0, 1),
    "q3_b": (4, -1),
    "q3_c": (10, 0),
    "infinity_1": (1, 0),
    "infinity_2": (5, 0),
}

torsion_generator = (K.zero, K.zero)
rational_generator = (K.convert(3), K.convert(5))
anti_generator = (K.convert(-8), K.from_sympy((-1 + 11 * r) / 2))
for label in direct_labels:
    expected_plus = multiply(
        plus_classes[label],
        torsion_generator,
        plus_coefficients,
    )
    assert points_equal(pushforwards[f"PLUS:{label}"], expected_plus), (
        label,
        normalize(pushforwards[f"PLUS:{label}"]),
        normalize(expected_plus),
    )
    rational_coefficient, anti_coefficient = minus_classes[label]
    expected_minus = add(
        multiply(
            rational_coefficient,
            rational_generator,
            minus_coefficients,
        ),
        multiply(anti_coefficient, anti_generator, minus_coefficients),
        minus_coefficients,
    )
    assert points_equal(pushforwards[f"MINUS:{label}"], expected_minus), (
        label,
        normalize(pushforwards[f"MINUS:{label}"]),
        normalize(expected_minus),
    )

# The remaining q2-node prime is recovered without another degree-six
# elimination: div(q2(x)) says that its class is three times the infinity
# class minus the other three q2 classes.
plus_q2_sum = None
minus_q2_sum = None
for label in ("q2_a", "q2_b", "q2_c"):
    plus_q2_sum = add(
        plus_q2_sum,
        pushforwards[f"PLUS:{label}"],
        plus_coefficients,
    )
    minus_q2_sum = add(
        minus_q2_sum,
        pushforwards[f"MINUS:{label}"],
        minus_coefficients,
    )
plus_infinity_sum = add(
    pushforwards["PLUS:infinity_1"],
    pushforwards["PLUS:infinity_2"],
    plus_coefficients,
)
minus_infinity_sum = add(
    pushforwards["MINUS:infinity_1"],
    pushforwards["MINUS:infinity_2"],
    minus_coefficients,
)
plus_q2_node = add(
    multiply(3, plus_infinity_sum, plus_coefficients),
    inverse(plus_q2_sum),
    plus_coefficients,
)
minus_q2_node = add(
    multiply(3, minus_infinity_sum, minus_coefficients),
    inverse(minus_q2_sum),
    minus_coefficients,
)
assert points_equal(
    plus_q2_node,
    multiply(1, torsion_generator, plus_coefficients),
)
assert points_equal(
    minus_q2_node,
    multiply(2, rational_generator, minus_coefficients),
)

print(
    "PASS exact C6 boundary images: 19 residue traces and the q2-node relation"
)
print(
    "PASS C6 class rows: Eplus Z/5 and Eminus coordinates in (G0,G2)"
)

#!/usr/bin/env python3
"""Certify the angle signs for consecutive mirrored moments at orders five and six."""

from fractions import Fraction

from flint import arb, ctx


def real(value):
    if isinstance(value, Fraction):
        return arb(value.numerator) / value.denominator
    return arb(value)


def span(lower, upper):
    return real(lower).union(real(upper))


def rational_power(base, numerator, denominator):
    if base.is_zero():
        return arb(0)
    return base.root(denominator)**numerator


def mobius(value, parameter):
    parameter = real(parameter)
    return ((1 + parameter) * value + 1 - parameter) / (
        (1 - parameter) * value + 1 + parameter
    )


def quintic_shape_upper(u):
    return u * mobius(u, Fraction(1, 2))


def sextic_outer_lower(u):
    return rational_power(u, 5, 4)


def sextic_outer_upper(u):
    return 24 * u / (25 - u)


def quadratic_roots(linear, constant):
    discriminant = (linear * linear - 4 * constant).sqrt()
    return (-linear - discriminant) / 2, (-linear + discriminant) / 2


def cubic_at_root(root, quadratic, linear, constant):
    center = root.mid()
    center_value = ((center + quadratic) * center + linear) * center + constant
    derivative = 3 * root * root + 2 * quadratic * root + linear
    return center_value + derivative * (root - center)


def direct_angle_signs(u, w):
    u_lower = u.lower().max(0)
    u_upper = u.upper().min(1)
    v = span(rational_power(u_lower, 3, 2), quintic_shape_upper(u_upper))
    left, right = quadratic_roots(5 * u - 1, 10 * v - 5 * u - 1)
    if not left.is_finite() or not right.is_finite():
        return False

    w_lower = w.lower().max(0)
    w_upper = w.upper().min(1)
    y_corners = (
        span(rational_power(w_lower, 8, 5), rational_power(w_upper, 8, 5)),
        span(w_lower * mobius(w_lower, Fraction(3, 5)), w_upper * mobius(w_upper, Fraction(3, 5))),
    )
    t_corners = (
        span(rational_power(w_lower, 1, 5), rational_power(w_upper, 1, 5)),
        span(mobius(w_lower, Fraction(1, 5)), mobius(w_upper, Fraction(1, 5))),
    )
    for y in y_corners:
        for t in t_corners:
            quadratic = 6 * w
            linear = 15 * y - 3
            constant = 20 * y * t - 12 * w
            if not (
                cubic_at_root(left, quadratic, linear, constant) > 0
                and cubic_at_root(right, quadratic, linear, constant) < 0
            ):
                return False
    return True


def v_lower_remainder(p):
    return ((1 - p).root(2)**3 - 1 + 3 * p / 2) / (p * p)


def y_lower_remainder(q):
    return ((1 - q).root(5)**8 - 1 + 8 * q / 5 - 12 * q * q / 25) / q**3


def t_lower_remainder(q):
    return ((1 - q).root(5) - 1 + q / 5 + 2 * q * q / 25) / q**3


def y_upper_remainder(q):
    return 12 / (25 * (5 - q))


def t_upper_remainder(q):
    return -4 / (25 * (5 - 2 * q))


def scaled_cubic(q, ratio, y_remainder, t_remainder):
    y_quadratic = real(Fraction(12, 25)) + y_remainder * q
    y = 1 - 8 * q / 5 + y_quadratic * q * q
    constant_remainder = -10 * y_remainder + 20 * (
        y * t_remainder
        - y_quadratic / 5
        - real(Fraction(2, 25)) * (-real(Fraction(8, 5)) + y_quadratic * q)
    )
    return -6 * ratio, 15 * ratio * ratio * y_quadratic, ratio**3 * constant_remainder


def scaled_angle_signs(u, w):
    p = 1 - u
    q = 1 - w
    if p.lower() <= 0 or q.lower() <= 0 or w.lower() <= 0:
        return False

    ratio = q / p
    ratio_upper = real(Fraction(5, 4))
    if ratio.lower() > ratio_upper or ratio.upper() < 1:
        return True
    ratio = span(ratio.lower().max(1), ratio.upper().min(ratio_upper))
    v_remainder = span(v_lower_remainder(p.lower()), 3 / (2 * (4 - p.upper())))
    left, right = quadratic_roots(arb(-5), 10 * v_remainder)
    if not left.is_finite() or not right.is_finite():
        return False

    q_lower, q_upper = q.lower(), q.upper()
    y_corners = (
        span(y_lower_remainder(q_lower), y_lower_remainder(q_upper)),
        span(y_upper_remainder(q_lower), y_upper_remainder(q_upper)),
    )
    t_corners = (
        span(t_lower_remainder(q_upper), t_lower_remainder(q_lower)),
        span(t_upper_remainder(q_upper), t_upper_remainder(q_lower)),
    )
    for y_remainder in y_corners:
        for t_remainder in t_corners:
            cubic = scaled_cubic(q, ratio, y_remainder, t_remainder)
            if not (cubic_at_root(left, *cubic) > 0 and cubic_at_root(right, *cubic) < 0):
                return False
    return True


def near_flat_angle_signs(p, ratio):
    q = p * ratio
    w_lower = (1 - q).lower()
    if w_lower <= 0:
        return False
    v_remainder = span(Fraction(3, 8), 3 / (2 * (4 - p.upper())))
    y_remainder = span(Fraction(8, 125), y_upper_remainder(q.upper()))
    t_remainder = span(-6 / (125 * w_lower**3), t_upper_remainder(q.lower()))
    left, right = quadratic_roots(arb(-5), 10 * v_remainder)

    quadratic, linear, constant = scaled_cubic(q, ratio, y_remainder, t_remainder)

    def cubic(root):
        return root**3 + quadratic * root * root + linear * root + constant

    return cubic(left) > 0 and cubic(right) < 0


def outside_cross_order_bounds(u, w):
    return (
        w.upper() < sextic_outer_lower(u.lower().max(0))
        or w.lower() > sextic_outer_upper(u.upper().min(1))
    )


def verify_boxes(ranges, verify_signs, verify_cross_order_bounds):
    stack = [ranges]
    visited = 0
    while stack:
        bounds = stack.pop()
        first, second = (span(*interval) for interval in bounds)
        visited += 1
        if verify_cross_order_bounds and outside_cross_order_bounds(first, second):
            continue
        if verify_signs(first, second):
            continue
        widths = [upper - lower for lower, upper in bounds]
        if visited >= 50_000 or max(widths) < Fraction(1, 10**12):
            raise ArithmeticError(("unresolved angle box", bounds))
        index = max(range(2), key=lambda k: widths[k])
        lower, upper = bounds[index]
        midpoint = (lower + upper) / 2
        stack.append(bounds[:index] + ((lower, midpoint),) + bounds[index + 1:])
        stack.append(bounds[:index] + ((midpoint, upper),) + bounds[index + 1:])
    return visited


def main():
    ctx.prec = 110
    direct_upper = Fraction(13, 20)
    scaled_upper = Fraction(9, 10)
    near_flat_p = Fraction(13, 100)
    if 1 - near_flat_p > scaled_upper:
        raise ArithmeticError("The scaled and near-flat ranges must overlap.")
    counts = (
        verify_boxes(((Fraction(0), direct_upper), (Fraction(0), direct_upper)), direct_angle_signs, True),
        verify_boxes(((direct_upper, scaled_upper), (Fraction(0), scaled_upper)), scaled_angle_signs, True),
        verify_boxes(
            ((Fraction(0), near_flat_p), (25 / (24 + near_flat_p), Fraction(5, 4))),
            near_flat_angle_signs,
            False,
        ),
    )
    print(f"Mirrored terminal-six angle signs passed. Box counts: {counts}.")


if __name__ == "__main__":
    main()

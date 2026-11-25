#!/usr/bin/env python3
"""Rational near-flat angle certificate for fully shared ordinary terminal six."""

from fractions import Fraction as Q

from verify_multivariate_terminal_six_private_angles import multiply_interval, require


ZERO = (Q(0), Q(0))
ONE = (Q(1), Q(1))
SCALE_RANGE = (Q(0), Q(19, 60))
RATIO_RANGE = (Q(517, 500), Q(53, 50))
RAYS = ((Q(527, 500), Q(7, 20)), (Q(43, 50), Q(0)), (Q(53, 100), Q(0)), (Q(39, 100), Q(0)))


def add(*values):
    return sum(value[0] for value in values), sum(value[1] for value in values)


def neg(value):
    return -value[1], -value[0]


def scale(value, coefficient):
    if coefficient >= 0:
        return coefficient * value[0], coefficient * value[1]
    return coefficient * value[1], coefficient * value[0]


def symmetric(value, coefficient):
    radius = coefficient * value[1]
    return -radius, radius


def power(value, exponent):
    return value[0] ** exponent, value[1] ** exponent


def check_constants():
    pmax, qmax = Q(101, 1000), Q(107, 1000)
    require(SCALE_RANGE[1] ** 2 < pmax and RATIO_RANGE[1] * SCALE_RANGE[1] ** 2 < qmax,
            "scaled deficits")
    require((2 - pmax) / (4 * (1 - pmax)**2) < Q(588, 1000), "quintic reciprocal scale")
    require((2 - qmax) / (5 * (1 - qmax)**2) < Q(475, 1000), "sextic reciprocal scale")
    require((Q(588, 1000)**2 < Q(346, 1000)
             and Q(3, 5) * Q(588, 1000)**4 < Q(72, 1000)
             and (Q(346, 1000) + Q(72, 1000) * pmax**2) / 2 < Q(174, 1000)),
            "quintic skew powers")
    require((Q(5, 3) * Q(475, 1000)**2 < Q(377, 1000)
             and 6 * Q(475, 1000)**4 < Q(306, 1000)
             and Q(4, 5) * Q(377, 1000) + Q(306, 5000) * qmax**2 < Q(303, 1000)),
            "sextic skew powers")
    require((Q(501, 1000) * Q(346, 1000)**2 < Q(6, 100)
             and Q(501, 1000) * Q(174, 1000)**2 < Q(16, 1000)
             and Q(501, 1000) * Q(377, 1000)**2 < Q(72, 1000)
             and Q(501, 1000) * Q(303, 1000)**2 < Q(46, 1000)),
            "even skew corrections")
    require(Q(6, 5) * Q(475, 1000)**5 < Q(3, 100), "central sextic link")
    for prefactor, max_gap, cap in ((Q(1, 6), pmax, Q(23, 100)),
                                     (Q(1, 5), qmax, Q(282, 1000))):
        bound = prefactor / (1 - max_gap)**3
        require(bound / (1 - bound * max_gap**3) < cap, "mean defect cap")
    require((Q(1, 16)**2 / (1 - pmax)**3 < Q(74, 1000)**2
             and Q(74, 1000) + Q(23, 100) <= Q(304, 1000)),
            "quintic cubic correction")
    require((Q(14, 625)**5 < Q(3, 100)**5 * (1 - qmax)**12
             and Q(9, 625)**5 < Q(2, 100)**5 * (1 - qmax)**11),
            "sextic fourth-order Taylor remainders")
    require((Q(1, 5)**5 < Q(22, 100)**5 * (1 - qmax)**4
             and Q(1) < Q(118, 100)**5 * (1 - qmax)**7
             and Q(1) < Q(108, 100)**5 * (1 - qmax)**3),
            "central defect derivative bounds")
    require(Q(3, 2) * Q(22, 100) * Q(282, 1000)
            + Q(3, 8) * Q(118, 100) * Q(282, 1000)**2 * qmax**2
            + Q(108, 100) * Q(3, 100) * qmax < Q(98, 1000), "linked central defect")
    require(((1 - pmax) * (2 - pmax) / 48)**2
            > (RATIO_RANGE[0] - 1)**2 * (1 + pmax / 12)**3, "ratio lower bound")
    require(((1 - pmax)**-2 - 1) / 4 < Q(1, 16) and Q(20, 19) < RATIO_RANGE[1],
            "ratio upper bound")


def check_limit_anchor():
    square_root = (Q(3162, 1000), Q(3163, 1000))
    require(square_root[0] ** 2 < 10 < square_root[1] ** 2, "quintic radical")
    quintic_low = ((5 - square_root[1]) / 8, (5 - square_root[0]) / 8)
    quintic_high = ((5 + square_root[0]) / 8, (5 + square_root[1]) / 8)

    def cubic(value):
        return ((value - 6) * value + Q(36, 5)) * value - Q(24, 25)

    roots = [(Q(a, 1000), Q(a + 1, 1000)) for a in (152, 1427, 4420)]
    for lower, upper in roots:
        require(cubic(lower) * cubic(upper) < 0, ("sextic anchor interval", lower, upper))
    sextic_low, sextic_middle, sextic_high = [
        scale(multiply_interval(RATIO_RANGE, interval), Q(1, 4)) for interval in roots
    ]
    intervals = (sextic_high, quintic_high, sextic_middle, quintic_low, sextic_low)
    for (high, low), (ray, _) in zip(zip(intervals, intervals[1:]), RAYS, strict=True):
        require(high[0] > ray**2 > low[1], ("anchor ray order", ray, high, low))


def coefficients(order, lam, ratio):
    p = power(lam, 2)
    q = multiply_interval(ratio, p)
    if order == 5:
        leading = add(ONE, scale(p, Q(-5, 4)), scale(power(p, 2), Q(15, 64)),
                      (Q(0), Q(19, 100) * power(p, 3)[1] + Q(29, 1000) * power(p, 4)[1]))
        remainder = (Q(3, 8), Q(3, 8) + Q(304, 1000) * p[1])
        cubic = add((Q(5, 4), Q(5, 4)), scale(multiply_interval(p, remainder), Q(-5, 4)),
                    (Q(-20, 1000) * power(p, 3)[1], Q(38, 1000) * power(p, 3)[1]))
        linear = add(scale(remainder, Q(5, 8)),
                     (Q(-57, 1000) * power(p, 2)[1], Q(1, 100) * power(p, 2)[1]))
        return [symmetric(lam, Q(58, 1000)), linear, symmetric(lam, Q(436, 1000)),
                cubic, symmetric(power(lam, 3), Q(435, 1000)), leading]
    y_remainder = (Q(64, 1000), Q(64, 1000) + Q(3, 100) * q[1])
    z_remainder = (Q(48, 1000), Q(48, 1000) + Q(2, 100) * q[1])
    first_defect = (Q(0), Q(282, 1000))
    leading = add(
        ONE, scale(q, Q(-3, 2)), scale(power(q, 2), Q(9, 20)),
        multiply_interval(power(q, 3), scale(add(scale(y_remainder, 15),
                                                  scale(z_remainder, 10),
                                                  scale(first_defect, 30)), Q(1, 32))),
        (Q(-31, 1000) * power(q, 4)[1], Q(68, 1000) * power(q, 4)[1]),
    )
    quartic = multiply_interval(ratio, add(
        (Q(3, 2), Q(3, 2)), scale(q, Q(-9, 10)),
        multiply_interval(power(q, 2), scale(add(y_remainder, scale(z_remainder, 2),
                                                   scale(first_defect, 4)), Q(-15, 32))),
        (Q(-115, 1000) * power(q, 3)[1], Q(163, 1000) * power(q, 3)[1]),
    ))
    quadratic = multiply_interval(power(ratio, 2), add(
        (Q(9, 20), Q(9, 20)),
        multiply_interval(q, scale(add(scale(z_remainder, 2), neg(y_remainder),
                                         scale(first_defect, 2)), Q(15, 32))),
        (Q(-182, 1000) * power(q, 2)[1], Q(92, 1000) * power(q, 2)[1]),
    ))
    constant = multiply_interval(power(ratio, 3), scale(add(
        scale(y_remainder, 15), scale(z_remainder, -10),
        (Q(-1412, 1000) * q[1], Q(167, 100) * q[1]),
    ), Q(1, 32)))
    return [
        constant,
        symmetric(multiply_interval(power(ratio, 3), lam), Q(179, 1000)),
        quadratic,
        symmetric(multiply_interval(power(ratio, 2), lam), Q(57, 100)),
        quartic,
        symmetric(multiply_interval(power(ratio, 2), power(lam, 3)), Q(57, 100)),
        leading,
    ]


def ray_imaginary(lam, real, ray):
    # On the physical ray, (Y-b)(1+lam^2 b Y)=-b lam^2 real^2.
    square = power(lam, 2)
    b, slope = ray
    moving = add((b, b), scale(square, slope))
    num = multiply_interval(moving, multiply_interval(square, power(real, 2)))
    imaginary = (Q(0), moving[1])
    for _ in range(2):
        denominator = add(ONE, multiply_interval(moving, multiply_interval(square, imaginary)))
        lower = moving[0] - num[1] / denominator[0]
        upper = moving[1] - num[0] / denominator[1]
        imaginary = max(imaginary[0], lower), min(imaginary[1], upper)
    return imaginary


def verify_projection(order, ray, lam, ratio, real):
    imag = ray_imaginary(lam, real, ray)
    square = add(power(real, 2), power(imag, 2))
    imag_powers = [ZERO]
    real_power = ONE
    for _ in range(order):
        real_power, imag_power = (
            add(multiply_interval(real_power, real), neg(multiply_interval(imag_powers[-1], imag))),
            add(multiply_interval(real_power, imag), multiply_interval(imag_powers[-1], real)),
        )
        imag_powers.append(imag_power)
    modulus_powers = [power(square, k) for k in range(order + 1)]
    coeff = coefficients(order, lam, ratio)
    for rotation in range(order + 1):
        pieces = []
        for k, c in enumerate(coeff):
            sine = imag_powers[k - rotation] if k >= rotation else neg(imag_powers[rotation - k])
            pieces.append(multiply_interval(c, multiply_interval(modulus_powers[min(k, rotation)], sine)))
        low, high = add(*pieces)
        if low > 0 or high < 0:
            return True
    return False


def verify_ray(order, ray):
    # All scaled roots have modulus less than two by leading-term dominance.
    full_ratio = RATIO_RANGE if order == 6 else ONE
    coeff = coefficients(order, SCALE_RANGE, full_ratio)
    require(coeff[-1][0] * 2**order
            > sum(max(map(abs, value)) * 2**k for k, value in enumerate(coeff[:-1])),
            ("root bound", order))
    intervals = (SCALE_RANGE, full_ratio, (Q(0), Q(2)))
    widths = tuple(b - a or Q(1) for a, b in intervals)
    stack = [(intervals, 0)]
    seen = 0
    while stack:
        box, depth = stack.pop()
        seen += 1
        if verify_projection(order, ray, *box):
            continue
        require(depth < 24, ("unresolved scaled ray", order, ray, box))
        k = max(range(3), key=lambda i: (box[i][1] - box[i][0]) / widths[i])
        a, b = box[k]
        midpoint = (a + b) / 2
        stack.append((box[:k] + ((a, midpoint),) + box[k + 1:], depth + 1))
        stack.append((box[:k] + ((midpoint, b),) + box[k + 1:], depth + 1))
    return seen


if __name__ == "__main__":
    check_constants()
    check_limit_anchor()
    for order in (5, 6):
        print(f"Order {order}: {[verify_ray(order, ray) for ray in RAYS]} scaled rectangles.", flush=True)

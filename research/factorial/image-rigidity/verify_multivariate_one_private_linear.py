"""Exact phase certificates for two regimes with private degree one."""

import argparse
from fractions import Fraction as Q
from math import factorial

SCALE = 10**8
PI_UPPER = Q(22, 7)


def upper_rational(value):
    numerator = SCALE * value.numerator
    return Q((numerator + value.denominator - 1) // value.denominator, SCALE)


def upper_root(value, degree):
    # Return the least positive integer c with (c/SCALE)**degree >= value.
    target = value.numerator * SCALE**degree
    denom = value.denominator
    low, high = 0, SCALE
    while high**degree * denom < target:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if middle**degree * denom < target:
            low = middle
        else:
            high = middle
    return high


def endpoint_envelopes(order, shared_curvature):
    floor = order // 2
    right_degree = order - floor - 1
    lower_norm = Q(upper_root(Q(factorial(order)), order) - 1, SCALE)
    normalized = []
    for k in range(1, order):
        factorial_bound = Q(factorial(k)) / lower_norm**k
        if shared_curvature:
            coefficient = factorial_bound * Q(1, 2)**min(k, order - k)
        else:
            coefficient = factorial_bound**2
        normalized.append(upper_rational(coefficient))
    left = normalized[:floor]
    right = normalized[-right_degree:][::-1]
    assert sum(left) < 1 and sum(right) < 1, order
    radius = Q(upper_root((1 + sum(left))/(1 - sum(right)), order - floor), SCALE)
    inverse = Q(upper_root((1 + sum(right))/(1 - sum(left)), order - right_degree), SCALE)
    left = [upper_rational(c * radius**j) for j, c in enumerate(left, 1)]
    right = [upper_rational(c * inverse**j) for j, c in enumerate(right, 1)]
    assert sum(left) < 1 and sum(right) < 1, order
    return left, right


def sine_lower(x):
    return x - x**3/6


def cosine_upper(x):
    return 1 - x*x/2 + x**4/24


def rotated_upper(coeffs, angle):
    result = Q()
    for j, c in enumerate(coeffs, 1):
        ratio = (j - angle) / angle
        if ratio < 1:
            # This branch is used only for 0 < theta <= pi/2.
            ratio /= 1 - (angle * PI_UPPER / 2)**2/6
        result += ratio*c
    return result


def rotated_lower(coeffs, angle, sine_at_pi, sine_at_split):
    near_zero = (
        2 * PI_UPPER / (9 * sine_at_pi)
        * sum((j + angle) * c for j, c in enumerate(coeffs, 1) if j >= 2)
    )
    near_pi = coeffs[0] + sum(coeffs[1:]) / sine_at_split
    return max(near_zero, near_pi)


def small_order_checks(order, left, right, shared_curvature):
    # Columns give the two small-angle positive phase coefficients,
    # and the two absolute bounds for the negative large-angle phases.
    if shared_curvature:
        bounds = {
            4: (".46", ".26", ".28", "0"),
            5: (".35", ".40", ".16", ".19"),
            6: (".40", ".36", ".18", ".15"),
            7: (".31", ".40", ".11", ".18"),
        }
        upper_slope, lower_angle, upper_angle = map(Q, (".76", ".35", "1.17"))
    else:
        bounds = {
            4: (".61", ".29", ".485", "0"),
            5: (".36", ".45", ".18", ".27"),
            6: (".44", ".38", ".20", ".17"),
            7: (".25", ".43", ".09", ".20"),
        }
        upper_slope, lower_angle, upper_angle = map(Q, (".90", ".485", "1.06"))
    left_angle, right_angle, left_negative, right_negative = map(Q, bounds[order])
    assert left_angle + right_angle <= upper_slope
    assert rotated_upper(left, left_angle) < 1, order
    assert rotated_upper(right, right_angle) < 1, order
    negative_slope = Q(2, 9) * sum(
        sum(j*c for j, c in enumerate(coeffs, 1) if j >= 3) / (1 - sum(coeffs))
        for coeffs in (left, right)
    )
    assert negative_slope < Q(9, 100), order

    assert left_negative + right_negative <= lower_angle
    for coeffs, angle in ((left, left_negative), (right, right_negative)):
        tail = sum(coeffs[1:]) / (1 - coeffs[0])
        assert tail == 0 or tail < sine_lower(angle), order
    if order == 4 and not shared_curvature:
        # On pi/2 <= theta <= pi, the degree-two sine term is nonpositive.
        assert left[0] / (1 - sum(left)) < sine_lower(Q(".65")) / cosine_upper(Q(".65"))
        assert right[0] < sine_lower(Q(".41"))
    else:
        # The j=2 term never has positive imaginary part in this band.
        positive_bound = sum(
            (coeffs[0] + sum(coeffs[2:])) / (1 - sum(coeffs))
            for coeffs in (left, right)
        )
        assert positive_bound < upper_angle, (order, positive_bound)


def verify_range(shared_curvature):
    if shared_curvature:
        left_upper, right_upper, angle, sine_at_pi, sine_at_split = map(
            Q, (".31", ".40", ".14", ".42", ".37")
        )
        label = "private degrees 1 and 2 with the shared-curvature bound"
    else:
        left_upper, right_upper, angle, sine_at_pi, sine_at_split = map(
            Q, (".25", ".43", ".15", ".45", ".398")
        )
        label = "private degrees 1 and at least 3"
    pi_lower = Q(157, 50)
    assert sine_lower(angle * pi_lower) > sine_at_pi
    assert sine_lower(angle * pi_lower / (1 + angle)) > sine_at_split
    for order in range(4, 361):
        left, right = endpoint_envelopes(order, shared_curvature)
        if order <= 7:
            small_order_checks(order, left, right, shared_curvature)
        if order >= 7:
            assert rotated_upper(left, left_upper) < 1, (order, "left upper phase")
            assert rotated_upper(right, right_upper) < 1, (order, "right upper phase")
            for coeffs in (left, right):
                assert rotated_lower(coeffs, angle, sine_at_pi, sine_at_split) < 1, order
    print(f"Certified all 357 phase envelopes at orders 4 through 360 for {label}.")


def main():
    if not __debug__:
        raise SystemExit("Verification requires assertions. Run Python without -O.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--shared-curvature", action="store_true", help="run only the shared-curvature checks"
    )
    args = parser.parse_args()
    if not args.shared_curvature:
        verify_range(shared_curvature=False)
    verify_range(shared_curvature=True)


if __name__ == "__main__":
    main()

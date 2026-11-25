#!/usr/bin/env python3
"""Exact ray exclusions for terminal-six private-degree and curvature bounds."""

from fractions import Fraction
from math import comb, factorial


RAYS = (Fraction(81, 100), Fraction(2, 5), Fraction(-7, 40), Fraction(-17, 25))
CONCENTRATED_RAYS = (Fraction(41, 50),) + RAYS[1:]
THREE_PRIVATE_RAYS = (Fraction(2, 5), Fraction(-1, 5), Fraction(-7, 10))
THREE_PRIVATE_UPPER_RAY = Fraction(39, 50)
DENOMINATOR = 10000


def require(condition, details):
    if not condition:
        raise ArithmeticError(details)


def exponential_lower(value):
    total = term = Fraction(1)
    for index in range(1, 65):
        term *= value / index
        total += term
    return total


def rational_cap(power_bound, order):
    lower, upper = 0, DENOMINATOR * comb(order, order // 2)
    require(Fraction(upper, DENOMINATOR) ** order >= power_bound, ("cap bracket", order))
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if Fraction(middle, DENOMINATOR) ** order >= power_bound:
            upper = middle
        else:
            lower = middle
    return Fraction(upper, DENOMINATOR)


def coefficient_bounds(order, private_rows, other_curvature, reversed_rows=()):
    coefficients = [Fraction(1)]
    for index in range(1, order):
        power_bound = Fraction(comb(order, index) ** order)
        for endpoint_index, gaps in ((index, private_rows), (order - index, reversed_rows)):
            for gap in gaps:
                power_bound *= Fraction(
                    factorial(gap * endpoint_index) ** order,
                    factorial(gap * order) ** endpoint_index,
                )
        exponential_argument = other_curvature * index * (order - index) / 2
        power_bound /= exponential_lower(exponential_argument)
        coefficients.append(rational_cap(power_bound, order))
    return tuple(coefficients) + (Fraction(1),)


def verify_root_annulus(coefficients):
    order = len(coefficients) - 1
    require(sum(c / 2**k for k, c in enumerate(coefficients) if k) < 1,
            ("inner-circle domination", coefficients))
    require(sum(c * 2**k for k, c in enumerate(coefficients[:-1])) < 2**order,
            ("outer-circle domination", coefficients))


def sine_ratios(order, cosine):
    ratios = [Fraction(0), Fraction(1)]
    for index in range(2, order + 1):
        ratios.append(2 * cosine * ratios[index - 1] - ratios[index - 2])
    return ratios


def projection_bounds(coefficients, cosine):
    order = len(coefficients) - 1
    ratios = sine_ratios(order, cosine)
    projections = []
    for rotation in range(order + 1):
        signs = [ratios[k - rotation] if k >= rotation else -ratios[rotation - k]
                 for k in range(order + 1)]
        lower = [signs[0]]
        upper = [signs[0]]
        for factor, coefficient in zip(signs[1:-1], coefficients[1:-1], strict=True):
            product = factor * coefficient
            lower.append(min(Fraction(0), product))
            upper.append(max(Fraction(0), product))
        lower.append(signs[-1])
        upper.append(signs[-1])
        projections.append((lower, upper))
    return projections


def verify_ray(coefficients, cosine):
    order = len(coefficients) - 1
    projections = projection_bounds(coefficients, cosine)
    stack = [(Fraction(1, 2), Fraction(2), 0)]
    boxes = 0
    while stack:
        left, right, depth = stack.pop()
        boxes += 1
        left_powers = [left**k for k in range(order + 1)]
        right_powers = [right**k for k in range(order + 1)]
        for lower, upper in projections:
            low = sum(c * (left_powers[k] if c >= 0 else right_powers[k])
                      for k, c in enumerate(lower))
            high = sum(c * (right_powers[k] if c >= 0 else left_powers[k])
                       for k, c in enumerate(upper))
            if low > 0 or high < 0:
                break
        else:
            require(depth < 18, ("unresolved ray interval", cosine, left, right))
            midpoint = (left + right) / 2
            stack.append((left, midpoint, depth + 1))
            stack.append((midpoint, right, depth + 1))
    return boxes


def verify_case(private_rows, other_curvature, rays, extra_sextic_ray=None):
    total = 0
    for order in (5, 6):
        coefficients = coefficient_bounds(order, private_rows, other_curvature)
        verify_root_annulus(coefficients)
        total += sum(verify_ray(coefficients, cosine) for cosine in rays)
        if order == 6 and extra_sextic_ray is not None:
            total += verify_ray(coefficients, extra_sextic_ray)
    return total


def multiply_interval(first, second):
    products = tuple(x * y for x in first for y in second)
    return min(products), max(products)


def sine_ratio_intervals(order, cosines):
    ratios = [(Fraction(0), Fraction(0)), (Fraction(1), Fraction(1))]
    twice_cosines = (2 * cosines[0], 2 * cosines[1])
    for index in range(2, order + 1):
        lower, upper = multiply_interval(twice_cosines, ratios[index - 1])
        previous = ratios[index - 2]
        ratios.append((lower - previous[1], upper - previous[0]))
    return ratios


def verify_angle_radius_rectangle(coefficients, cosines, radii):
    order = len(coefficients) - 1
    cosine_width = cosines[1] - cosines[0]
    radius_width = radii[1] - radii[0]
    stack = [(cosines, radii, 0)]
    boxes = 0
    while stack:
        x, r, depth = stack.pop()
        boxes += 1
        ratios = sine_ratio_intervals(order, x)
        powers = [(r[0] ** k, r[1] ** k) for k in range(order + 1)]
        for rotation in range(order + 1):
            lower = upper = Fraction(0)
            for k, coefficient in enumerate(coefficients):
                if k >= rotation:
                    sign = ratios[k - rotation]
                else:
                    lo, hi = ratios[rotation - k]
                    sign = (-hi, -lo)
                value = multiply_interval(sign, powers[k])
                if 0 < k < order:
                    value = multiply_interval(value, (Fraction(0), coefficient))
                lower += value[0]
                upper += value[1]
            if lower > 0 or upper < 0:
                break
        else:
            require(depth < 22, ("unresolved angle-radius interval", order, x, r))
            if (x[1] - x[0]) * radius_width > (r[1] - r[0]) * cosine_width:
                midpoint = (x[0] + x[1]) / 2
                stack.append(((x[0], midpoint), r, depth + 1))
                stack.append(((midpoint, x[1]), r, depth + 1))
            else:
                midpoint = (r[0] + r[1]) / 2
                stack.append((x, (r[0], midpoint), depth + 1))
                stack.append((x, (midpoint, r[1]), depth + 1))
    return boxes


def verify_private_degree_three_inequalities():
    cosines = (Fraction(39, 50), Fraction(81, 100))
    boxes = sum((
        verify_angle_radius_rectangle(
            coefficient_bounds(5, (1, 1, 1), Fraction(0)),
            cosines, (Fraction(11, 10), Fraction(2)),
        ),
        verify_angle_radius_rectangle(
            coefficient_bounds(6, (1, 1, 1), Fraction(0)),
            cosines, (Fraction(1, 2), Fraction(17, 20)),
        ),
    ))
    constants = (
        (1, 1, Fraction(4, 5)), (2, 2, Fraction(2, 3)),
        (3, 3, Fraction(3, 5)), (4, 3, Fraction(11, 10)),
    )
    for private_rows in ((1, 1, 1), (2, 1), (3,)):
        for upper, lower, limit in constants:
            ratio_power = Fraction(comb(6, upper), comb(5, lower)) ** 30
            for gap in private_rows:
                ratio_power *= Fraction(factorial(gap * upper), factorial(gap * lower)) ** 30
                ratio_power *= Fraction(factorial(5 * gap) ** (6 * lower),
                                        factorial(6 * gap) ** (5 * upper))
            require(ratio_power < limit**30, ("private link constant", private_rows, upper))
    radius_ratio = Fraction(22, 17)
    reciprocal_sextic_radius = Fraction(20, 17)
    multiplier = Fraction(31, 10)
    require(max(
        Fraction(4, 5) * radius_ratio**4 * reciprocal_sextic_radius,
        Fraction(2, 3) * radius_ratio**3 * reciprocal_sextic_radius,
        radius_ratio**2 * (Fraction(3, 5) * reciprocal_sextic_radius + Fraction(11, 10)),
    ) < multiplier, "projected imaginary contribution")
    # Verify -sin(6 theta) + (31/10) sin(5 theta) > 0 after dividing by sin(theta).
    for index in range(100):
        left = cosines[0] + (cosines[1] - cosines[0]) * index / 100
        right = cosines[0] + (cosines[1] - cosines[0]) * (index + 1) / 100
        ratios = sine_ratio_intervals(6, (left, right))
        require(-ratios[6][1] + multiplier * ratios[5][0] > 0,
                ("angular sine inequality", left, right))
    return boxes


def verify_two_sided_private_degrees():
    both_two = 0
    for order in (5, 6):
        bounds = coefficient_bounds(order, (1, 1), Fraction(0), (1, 1))
        verify_root_annulus(bounds)
        both_two += sum(verify_ray(bounds, ray) for ray in RAYS)

    quintic = coefficient_bounds(5, (1, 1), Fraction(0), (1,))
    sextic = coefficient_bounds(6, (1, 1), Fraction(0), (1,))
    for bounds in (quintic, sextic):
        verify_root_annulus(bounds)
    two_and_one = sum(verify_ray(bounds, ray)
                      for bounds in (quintic, sextic) for ray in RAYS[1:])
    band = (Fraction(157, 200), Fraction(81, 100))
    two_and_one += verify_ray(sextic, band[0])
    two_and_one += verify_angle_radius_rectangle(
        quintic, band, (Fraction(109, 100), Fraction(2)))
    two_and_one += verify_angle_radius_rectangle(
        sextic, band, (Fraction(1, 2), Fraction(9, 10)))

    ratios = ((1, 1, Fraction(7, 8)), (2, 2, Fraction(4, 5)),
              (3, 3, Fraction(2, 3)), (4, 3, Fraction(9, 10)))
    for partition in ((1, 1), (2,)):
        for upper, lower, limit in ratios:
            ratio_power = Fraction(comb(6, upper), comb(5, lower)) ** 30
            for gaps, numerator_index, denominator_index in (
                (partition, upper, lower), ((1,), 6 - upper, 5 - lower)
            ):
                for gap in gaps:
                    ratio_power *= Fraction(
                        factorial(gap * numerator_index), factorial(gap * denominator_index)
                    ) ** 30
                    ratio_power *= Fraction(
                        factorial(5 * gap) ** (6 * denominator_index),
                        factorial(6 * gap) ** (5 * numerator_index)
                    )
            require(ratio_power < limit**30, ("two-sided link", partition, upper))

    radius_ratio = Fraction(109, 90)
    inverse_radius = Fraction(10, 9)
    multiplier = Fraction(5, 2)
    require(max(
        Fraction(7, 8) * radius_ratio**4 * inverse_radius,
        Fraction(4, 5) * radius_ratio**3 * inverse_radius,
        radius_ratio**2 * (Fraction(2, 3) * inverse_radius + Fraction(9, 10)),
    ) < multiplier, "two-sided imaginary contribution")
    for index in range(100):
        left = band[0] + (band[1] - band[0]) * index / 100
        right = band[0] + (band[1] - band[0]) * (index + 1) / 100
        values = sine_ratio_intervals(6, (left, right))
        require(-values[6][1] + multiplier * values[5][0] > 0,
                ("two-sided sine inequality", left, right))
    return both_two, two_and_one


def verify_total_curvature_six():
    """The middle sextic coefficient is bounded from both endpoints."""
    bounds = {j: coefficient_bounds(j, (), Fraction(6)) for j in (5, 6)}
    rays = Fraction(2, 5), Fraction(-11, 50), Fraction(-17, 25)
    boxes = 0
    for coefficients in bounds.values():
        verify_root_annulus(coefficients)
        boxes += sum(verify_ray(coefficients, ray) for ray in rays)

    x_lo, x_hi = Fraction(3907, 5000), Fraction(8091, 10000)
    require(4 * x_hi - 1 > 0 and (4 * x_hi - 1) ** 2 > 5,
            "cosine of pi/5 is below the upper endpoint")
    r, s = Fraction(537, 500), Fraction(459, 500)
    boxes += verify_ray(bounds[6], x_lo)
    boxes += verify_angle_radius_rectangle(bounds[5], (x_lo, x_hi), (r, Fraction(2)))
    boxes += verify_angle_radius_rectangle(bounds[6], (x_lo, x_hi), (Fraction(1, 2), s))

    sine_ratio = Fraction(621, 1000)
    # sine_ratio * U_2(x) - U_3(x) is positive and decreasing on the band.
    require(8 * sine_ratio * x_hi - 24 * x_lo**2 + 4 < 0 and
            sine_ratio * (4 * x_hi**2 - 1) - (8 * x_hi**3 - 4 * x_hi) > 0,
            "sine ratio bound")
    radius_ratio = r / s
    outer_one = Fraction(6, 5) * radius_ratio**4 / s
    outer_two = Fraction(3, 2) * radius_ratio**3 / s
    inner_four = Fraction(3, 2) * radius_ratio**2 * sine_ratio
    middle_from_two = 2 * radius_ratio**3
    middle_from_three = 2 * radius_ratio**2 / s
    # The remaining coefficient is bounded by a minimum of two linear functions.
    at_intersection = (
        outer_two * middle_from_three
        + inner_four * middle_from_two
        + middle_from_two * middle_from_three
    ) / (middle_from_two + middle_from_three)
    multiplier = Fraction(3469, 1000)
    require(max(outer_one, outer_two, inner_four, at_intersection) < multiplier,
            "linked contribution threshold")

    # D(x) = -U_5(x) + multiplier * U_4(x) is positive and increasing.
    ratios = sine_ratios(6, x_lo)
    require(-ratios[6] + multiplier * ratios[5] > 0,
            "positive sine separation at the lower endpoint")
    derivative_lower = (-160 * x_hi**4 + 96 * x_lo**2 - 6
                        + multiplier * (64 * x_lo**3 - 24 * x_hi))
    require(derivative_lower > 0, "sine separation derivative")
    return boxes


if __name__ == "__main__":
    certificate = {
        f"weighted curvature, {degree} unit rows": verify_case(
            (1,) * degree, max(Fraction(0), Fraction(34, 5) - 2 * degree), RAYS
        )
        for degree in range(5)
    }
    certificate["one private gap three"] = verify_case((3,), Fraction(0), CONCENTRATED_RAYS)
    certificate["three private degrees, angular localization"] = verify_case(
        (1, 1, 1), Fraction(0), THREE_PRIVATE_RAYS, THREE_PRIVATE_UPPER_RAY
    )
    certificate["three private degrees, linked radius bounds"] = verify_private_degree_three_inequalities()
    both_two, two_and_one = verify_two_sided_private_degrees()
    certificate["private degrees two and two"] = both_two
    certificate["private degrees two and one"] = two_and_one
    certificate["total curvature at least six"] = verify_total_curvature_six()
    for name, boxes in certificate.items():
        print(f"{name}: {boxes} rational radius intervals")

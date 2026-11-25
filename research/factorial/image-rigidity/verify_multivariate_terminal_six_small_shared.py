#!/usr/bin/env python3
"""Exact terminal-six angles for one or two private degrees and low shared curvature."""

from fractions import Fraction
from math import comb, factorial, isqrt

from verify_multivariate_terminal_six_private_angles import (
    exponential_lower, rational_cap, require, sine_ratios,
)


CASES = {
    (1,): (
        Fraction(9, 100), (146, -196, -640, -771),
        ((24361, 317791), (-295067, 303015), (-458588, 0)),
        ((56123, 258241), (-180097, 304099), (-376026, 133474)),
    ),
    (1, 1): (
        Fraction(7, 8), (616, 160, -376, -720),
        ((66587, 112161), (-83645, 132698), (-165883, 0)),
        ((61553, 74189), (-31046, 111015), (-113841, 51593)),
    ),
    (2,): (
        Fraction(13, 10), (688, 200, -322, -710),
        ((25365, 34902), (-25664, 45111), (-54957, 0)),
        ((21653, 22110), (-8195, 35965), (-36186, 17107)),
    ),
}


def multiply_complex(first, second):
    x, y = first
    u, v = second
    return x * u - y * v, x * v + y * u


def absolute_upper(value):
    return sum(map(abs, value))


def check_anchor(order, private_rows, cosines, given_centers):
    raw = [
        comb(order, index) * product_of_factorials(private_rows, index)
        for index in range(order + 1)
    ]
    centers = [(Fraction(u, 10**6), Fraction(v, 10**6)) for u, v in given_centers]
    centers += [(u, -v) for u, v in centers if v > 0]
    radius = Fraction(1, 10000)
    require(len(centers) == order, ("anchor count", order, private_rows))
    for i, center in enumerate(centers):
        for other in centers[i + 1:]:
            require(max(abs(x - y) for x, y in zip(center, other, strict=True)) > 2 * radius,
                    ("anchor disk overlap", center, other))
        powers = [(Fraction(1), Fraction(0))]
        for k in range(1, order + 1):
            powers.append(multiply_complex(powers[-1], center))
        coefficients = [
            tuple(sum(raw[k] * comb(k, m) * powers[k - m][coordinate]
                      for k in range(m, order + 1))
                  for coordinate in range(2))
            for m in range(order + 1)
        ]
        require(
            max(map(abs, coefficients[1])) * radius
            > absolute_upper(coefficients[0]) + sum(
                absolute_upper(coefficients[m]) * radius**m for m in range(2, order + 1)
            ),
            ("anchor Rouché bound", order, private_rows, center),
        )
        u, v = center
        if v < 0:
            continue
        if v == 0:
            require(u < -radius, ("negative real anchor", center))
            continue
        require(v > radius, ("upper anchor disk", center))
        after = 0
        for cosine in cosines:
            square = 1 - cosine**2
            denominator = 10**7
            floor = isqrt((square.numerator * denominator**2) // square.denominator)
            sines = (Fraction(floor, denominator), Fraction(floor + 1, denominator))
            lo, hi = sorted(v * cosine - u * sine for sine in sines)
            if lo > radius:
                after += 1
            else:
                require(hi < -radius, ("anchor disk intersects separator", center, cosine))
        expected = (1, 3) if order == 5 else (0, 2, 4)
        require(after == expected[i],
                ("wrong anchor sector", order, center, after))


def product_of_factorials(private_rows, index):
    value = 1
    for gap in private_rows:
        value *= factorial(gap * index)
    return value


def bounds(order, private_rows, variance):
    lower = [Fraction(1)]
    upper = [Fraction(1)]
    for k in range(1, order):
        power = Fraction(comb(order, k) ** order)
        for gap in private_rows:
            power *= Fraction(factorial(gap * k) ** order, factorial(gap * order) ** k)
        high = rational_cap(power, order)
        upper.append(high)
        u = variance * k * (order - k) / 2
        require(u < 66, ("exponential tail bound", u))
        power /= exponential_lower(u) + u**65 / factorial(65) / (1 - u / 66)
        low = rational_cap(power, order) - Fraction(1, 10000)
        lower.append(max(Fraction(0), low))
    return lower + [Fraction(1)], upper + [Fraction(1)]


def verify_ray(order, lower, upper, cosine):
    ratios = sine_ratios(order, cosine)
    projections = []
    for rotation in range(order + 1):
        signs = [ratios[k - rotation] if k >= rotation else -ratios[rotation - k]
                 for k in range(order + 1)]
        projections.append((
            [min(x * sign, y * sign) for x, y, sign in zip(lower, upper, signs, strict=True)],
            [max(x * sign, y * sign) for x, y, sign in zip(lower, upper, signs, strict=True)],
        ))
    stack = [(Fraction(1, 5), Fraction(5), 0)]
    boxes = 0
    while stack:
        a, b, depth = stack.pop()
        boxes += 1
        a_powers = [a**k for k in range(order + 1)]
        b_powers = [b**k for k in range(order + 1)]
        for low, high in projections:
            if sum(c * (a_powers[k] if c >= 0 else b_powers[k]) for k, c in enumerate(low)) > 0:
                break
            if sum(c * (b_powers[k] if c >= 0 else a_powers[k]) for k, c in enumerate(high)) < 0:
                break
        else:
            require(depth < 20, ("unresolved ray interval", order, cosine, a, b))
            midpoint = (a + b) / 2
            stack.append((a, midpoint, depth + 1))
            stack.append((midpoint, b, depth + 1))
    return boxes


if __name__ == "__main__":
    for private_rows, (variance, cosines, fifth_centers, sixth_centers) in CASES.items():
        cosines = tuple(Fraction(cosine, 1000) for cosine in cosines)
        boxes = 0
        for order, centers in ((5, fifth_centers), (6, sixth_centers)):
            check_anchor(order, private_rows, cosines, centers)
            lower, upper = bounds(order, private_rows, variance)
            require(sum(c / 5**k for k, c in enumerate(upper) if k) < 1,
                    ("inner-circle domination", order, private_rows))
            require(sum(c * 5**k for k, c in enumerate(upper[:-1])) < 5**order,
                    ("outer-circle domination", order, private_rows))
            boxes += sum(verify_ray(order, lower, upper, x) for x in cosines)
        print(f"Private rows {private_rows}: endpoint curvature at most {variance}, {boxes} intervals.")

"""Certify the rational and radical thresholds for the analytic 287/500 bound."""

from fractions import Fraction as Q
from itertools import combinations
from math import isqrt


DENOMINATOR = 10**12


def root_lower(value):
    """Return the floor of the square root on the rational grid."""
    scaled_square = value.numerator * DENOMINATOR**2 // value.denominator
    return Q(isqrt(scaled_square), DENOMINATOR)


def root_difference_squared(a, b):
    """Upper bound the squared positive part of sqrt(a) minus sqrt(b)."""
    a, b = Q(a), Q(b)
    if a <= b:
        return Q(0)
    upper = root_lower(a) + Q(1, DENOMINATOR)
    lower = root_lower(b)
    return (upper - lower)**2


def check(name, bound, cap):
    cap = Q(cap)
    if bound >= cap:
        raise ValueError(f"{name}: bound {bound} exceeds strict threshold {cap}")
    print(f"{name}: < {cap}")


def polygon_vertices(i, j, k, low, high):
    """Enumerate the vertices of the square cut by ia+jb <= 1-k*low."""
    limit = 1 - k * low
    endpoints = (low, high)
    candidates = {(a, b) for a in endpoints for b in endpoints}
    candidates.update((a, (limit - i * a) / j) for a in endpoints)
    candidates.update(((limit - j * b) / i, b) for b in endpoints)
    return (
        (a, b) for a, b in candidates
        if low <= a <= high and low <= b <= high and i * a + j * b <= limit
    )


def main():
    target = Q(287, 500)
    x = 1 - target
    z = 2 * target - 1
    large = Q("0.1592")
    dominant = Q("0.283")
    discrepancy = Q("0.097")
    p = root_difference_squared

    check("large non-degree-one", max((3 + x) / 4 - dominant, 1 - z - dominant), target)
    check("tall degree-one", max((1 + x) / 2 - large + p(2 * x, "0.75"),
                                (1 + x) / 2 - Q("0.25") + p(2 * x, 1 - x)), "0.55706")
    check("tall degree-two or three, lower endpoint",
          (1 + x) / 2 - large + p(2 * x, 1 - 3 * large) / 2, "0.57386")
    check("tall degree-two or three, upper endpoint",
          (1 + x) / 2 - dominant + p(2 * x, 1 - 3 * dominant) / 2, "0.57283")

    if not (large > 8 * x / 25 and x - Q("0.25") > large and
            3 * x + large > 2 - target and 9 * large > 1 and x > Q("0.25")):
        raise ValueError("large tall coordinates: mass or geometry thresholds")
    check("two tall large singletons", x + p(2 * x, 4 * x - 1), "0.43306")
    check("third-largest mass at a high degree",
          Q("0.5") + large - 3 * discrepancy / 2 + p(2 * x, 5 * discrepancy), "0.56507")
    degree_caps = ("0.52557", "0.52557", "0.53881", "0.56676")
    for (i, j, k), cap in zip(combinations(range(1, 5), 3), degree_caps, strict=True):
        bound = max(
            (1 + x) / 2 - a - b + p(2 * x, 1 - i * a - j * b)
            for a, b in polygon_vertices(i, j, k, discrepancy, large)
        )
        check(f"third-largest mass for degrees {i}, {j}, {k}", bound, cap)

    g2 = p(3, 1 / x)
    g3 = p(4, Q("1.5") / x)
    short_balanced = max(
        x * g2,
        Q("0.06") * (1 - Q("0.28") / x) + (x - Q("0.25")) * g2,
        (x - Q("0.25"))**2 / (2 * x) + (x - Q("0.28")) * g3,
    )
    tall_balanced = max(x * g3, x * (3 - 1 / x)**2 / 24 + (x - Q(1, 6)) * g3)
    check("balanced short penalty", short_balanced, "0.03859")
    check("balanced tall penalty", tall_balanced, "0.01152")
    check("balanced exponent", Q("0.5") + discrepancy / 2 +
          (short_balanced + tall_balanced) / 2, target)

    if not 2 * large - discrepancy < target / 2:
        raise ValueError("need a tall letter with a large coordinate if the bins are unbalanced")
    check("unbalanced bins", Q("0.75") - (target - 2 * large + discrepancy) / 2, target)


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""Exact factor-shape bounds for alternative Eisenstein least-prime budgets."""

from bisect import bisect_right
from fractions import Fraction
from math import factorial, isqrt

from verify_eisenstein_budget_intersection import BETA, logarithm, show_interval, total_mass


def rectangle_mass(eta, cutoff, mesh, partner_bound):
    """Bound the mass of full rectangles admitted by a decreasing partner cap."""
    nodes = [eta] + [
        i * mesh for i in range(1, int(cutoff / mesh) + 1) if i * mesh > eta
    ]
    totals = [total_mass(node) for node in nodes]
    total_lo, total_hi = totals[0]

    lo = [Fraction(0)]
    for _, right in totals[1:]:
        lo.append(max(lo[-1], total_lo - right))
    hi = [total_hi - left for left, _ in totals]
    for i in range(len(hi) - 2, -1, -1):
        hi[i] = min(hi[i], hi[i + 1])
    hi[0] = Fraction(0)

    partner = [
        max(0, bisect_right(nodes, partner_bound(node)) - 1) for node in nodes[1:]
    ] + [0]

    def abel(cdf):
        return sum(
            cdf[i] * (cdf[partner[i - 1]] - cdf[partner[i]])
            for i in range(1, len(nodes))
        )

    return abel(lo), abel(hi)


def normalized_margin(eta, marked_lo, marked_hi):
    lo, hi = total_mass(eta)
    return 2 / hi - 1 + marked_lo / hi**2, 2 / lo - 1 + marked_hi / lo**2


ROUNDING_SCALE = 10**12


def sqrt_ceil(value):
    numerator = value.numerator * ROUNDING_SCALE**2
    denominator = value.denominator
    root = isqrt(numerator // denominator)
    return root + (root**2 * denominator < numerator)


def ceil_scaled(value):
    return (value.numerator * ROUNDING_SCALE + value.denominator - 1) // value.denominator


def lower_sieve_upper(ratio):
    if ratio <= 1:
        return 0
    if ratio >= 3:
        return ceil_scaled(Fraction(2, 3))
    argument = Fraction(sqrt_ceil(1 - 1 / ratio), ROUNDING_SCALE)
    series = sum(argument**(2*j + 1) / (2*j + 1) for j in range(16))
    tail = argument**33 / (33 * (1 - argument**2))
    inverse_root = Fraction(sqrt_ceil(1 / ratio), ROUNDING_SCALE)
    return ceil_scaled((series + tail) * inverse_root)


def unequal_vector_upper():
    cap = Fraction(19, 40)
    budget = 1 - cap
    start = Fraction(2, 25)
    x_mesh = Fraction(1, 400)
    upper_ratios = [Fraction(i, 100) for i in range(164)]
    roots = [sqrt_ceil(ratio) for ratio in upper_ratios[1:]]
    worst = 0
    for i in range(int((cap / 2 - start) / x_mesh) + 1):
        x_lo = start + i * x_mesh
        y_lo = cap - min(cap / 2, x_lo + x_mesh)
        lower_a = [lower_sieve_upper((budget - y_lo * u) / x_lo) for u in upper_ratios[:-1]]
        lower_b = [lower_sieve_upper((budget - x_lo * u) / y_lo) for u in upper_ratios[:-1]]
        local_max = max((
            root_a * lower_bound_a + root_b * lower_bound_b
            for root_a, lower_bound_b in zip(roots, lower_b) if lower_bound_b
            for root_b, lower_bound_a in zip(roots, lower_a) if lower_bound_a
        ), default=0)
        worst = max(worst, local_max)
    return Fraction(worst, ROUNDING_SCALE**2)


def main():
    eta = Fraction(1, 2) / BETA
    mesh = Fraction(1, 1000)

    def polygon_partner(x, slack):
        return min(
            Fraction(1, 2) + 2 * slack - x,
            1 + 4 * slack - 3 * x,
            (1 + 4 * slack - x) / 3,
        )

    inner = rectangle_mass(eta, Fraction(1, 2), mesh, lambda x: polygon_partner(x, 0))[0]
    outer = rectangle_mass(eta, Fraction(1, 2), mesh, lambda x: polygon_partner(x, mesh))[1]
    margin = normalized_margin(eta, inner, outer)
    assert margin[1] < -Fraction(261, 10000)
    print("common leastness cutoff", eta, show_interval(margin))

    cap = Fraction(47, 100)
    inner = rectangle_mass(eta, Fraction(1, 2), mesh, lambda x: cap - x)[0]
    outer = rectangle_mass(eta, Fraction(1, 2), mesh, lambda x: cap + 2 * mesh - x)[1]
    margin = normalized_margin(eta, inner, outer)
    assert margin[1] < -Fraction(89, 10000)
    print("normalized vector-sieve cap", eta, show_interval(margin))

    gamma_floor = Fraction(572, 1000)
    harmonic = sum(Fraction(1, i) for i in range(1, 101))
    log101_hi = 6 * logarithm(Fraction(2))[1] + logarithm(Fraction(101, 64))[1]
    assert harmonic - log101_hi > gamma_floor
    assert sum(gamma_floor**j / factorial(j) for j in range(5)) > Fraction(9*22, 16*7)
    ratio_ceiling = unequal_vector_upper()
    assert ratio_ceiling < Fraction(983, 1000)
    print("unequal upper/lower vector ratio ceiling", show_interval((ratio_ceiling, ratio_ceiling)))

    cap = Fraction(19, 40)
    inner = rectangle_mass(eta, Fraction(1, 2), mesh, lambda x: cap - x)[0]
    outer = rectangle_mass(eta, Fraction(1, 2), mesh, lambda x: cap + 2 * mesh - x)[1]
    margin = normalized_margin(eta, inner, outer)
    assert margin[1] < -Fraction(54, 10000)
    print("unequal vector-sieve cap", eta, show_interval(margin))

    eta, cap, mesh = Fraction(117, 1000), Fraction(499, 1000), Fraction(1, 2000)
    marked = rectangle_mass(eta, Fraction(1, 3), mesh, lambda x: cap - x)
    margin = normalized_margin(eta, *marked)
    assert margin[0] > Fraction(16, 10000)
    print("cross mark below one third", show_interval(margin))

    eta = Fraction(29, 250)
    raw_lo, raw_hi = total_mass(Fraction(1, 3))
    total_lo, total_hi = total_mass(eta)
    band = (raw_lo - 1) / total_hi, (raw_hi - 1) / total_lo
    assert band[0] > Fraction(143, 1000)
    print("first-side semiprime band above one third", show_interval(band))


if __name__ == "__main__":
    main()
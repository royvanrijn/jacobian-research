#!/usr/bin/env python3
"""Rational source-mass bounds for the rough Eisenstein divisor exchanges."""

from bisect import bisect_right
from fractions import Fraction

from verify_eisenstein_budget_intersection import BETA, logarithm, show_interval, total_mass


def balanced_triple_lower(eta, cap, grid):
    ceiling = 1 - cap
    step = (ceiling - eta) / grid
    minimum_sum = (1 - ceiling - 2 * eta) / step
    maximum_sum = (1 - 3 * eta) / step
    first_index_sum = -(-minimum_sum.numerator // minimum_sum.denominator)
    last_index_sum = maximum_sum.numerator // maximum_sum.denominator - 2
    precision = 10**14
    lower = 0
    for i in range(grid):
        x = eta + (i + Fraction(1, 2)) * step
        for j in range(max(0, first_index_sum - i), min(grid, last_index_sum - i + 1)):
            y = eta + (j + Fraction(1, 2)) * step
            middle = precision * step**2 / (24 * x * y * (1 - x - y))
            lower += middle.numerator // middle.denominator
    return Fraction(lower, precision)


def missed_composite_lower(eta, cap, mesh):
    nodes = [eta] + [
        i * mesh for i in range(1, int(Fraction(1, 2) / mesh) + 1)
        if i * mesh > eta
    ]
    totals = [total_mass(node) for node in nodes]
    upper_cdf = [totals[0][1] - low for low, _ in totals]
    for i in range(len(upper_cdf) - 2, -1, -1):
        upper_cdf[i] = min(upper_cdf[i], upper_cdf[i + 1])
    upper_cdf[0] = Fraction(0)

    lower_endpoints = nodes[:-1]
    partner = [
        bisect_right(lower_endpoints, cap - left) for left in lower_endpoints
    ] + [0]
    retained_upper = sum(
        upper_cdf[i] * (upper_cdf[partner[i - 1]] - upper_cdf[partner[i]])
        for i in range(1, len(nodes))
    )
    return (totals[0][0] - 1)**2 - retained_upper


def main():
    eta = Fraction(1, 2) / BETA
    mesh = Fraction(1, 4000)
    for cap, minimum in (
        (Fraction(101, 200), Fraction(1445, 1000)),
        (Fraction(277, 500), Fraction(100186, 100000)),
    ):
        missed = missed_composite_lower(eta, cap, mesh)
        triple = balanced_triple_lower(eta, cap, grid=768)
        capacity = missed + 4 * triple
        assert capacity > minimum
        print("cap", cap, "missed mass lower bound", show_interval((missed, missed)),
              "triple mass lower bound", show_interval((triple, triple)),
              "exchange capacity lower bound", show_interval((capacity, capacity)))


    eta, cap = Fraction(29, 250), Fraction(253, 500)
    companion_cap = logarithm((1 - cap) / eta)[1] / 2
    minimum = cap - eta
    restricted_mass = logarithm((1 - minimum) / minimum)[0] / 2
    scale = Fraction(3, 5)
    assert scale * companion_cap / 2 < restricted_mass
    assert 2 * scale * logarithm(1 / eta)[1]**2 / 8 < 1
    removed = 2 * scale * balanced_triple_lower(eta, cap, grid=768)
    assert removed > Fraction(3, 20)
    print("prime mass removed with both weight banks", show_interval((removed, removed)))


if __name__ == "__main__":
    main()
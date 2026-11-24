#!/usr/bin/env python3
"""Certified simplex bounds for Eisenstein factor counts and a finite check."""

from fractions import Fraction
from itertools import accumulate
from math import factorial, gcd

from experiment_rough_goldbach import (
    coefficients,
    integer_root,
    least_prime_factors,
    squarefree_data,
)


GRID = 2048
SCALE = 10**12


def simplex_sum(weights, dimensions, budget):
    """Sum products over nonnegative indices with total at most budget."""
    prefix = [0, *accumulate(weights)]

    def recurse(remaining, cap):
        if remaining == 1:
            return prefix[min(cap + 1, len(weights))]
        return sum(weights[index] * recurse(remaining - 1, cap - index)
                   for index in range(min(cap + 1, len(weights))))

    return recurse(dimensions, budget)


def factor_constant(eta, degree, grid=GRID):
    if degree == 1:
        return Fraction(1), Fraction(1)
    numerator, denominator = eta.numerator, eta.denominator
    width = denominator - degree * numerator
    if width <= 0:
        return Fraction(0), Fraction(0)

    dimension = degree - 1
    low_weights = []
    high_weights = []
    for index in range(grid):
        left = grid * numerator + index * width
        right = left + width
        low_weights.append(SCALE * width // right)
        high_weights.append((SCALE * width + left - 1) // left)

    normalization = factorial(dimension) * SCALE**dimension
    low = Fraction(simplex_sum(low_weights, dimension, grid - dimension), normalization)
    high = Fraction(simplex_sum(high_weights, dimension, grid - 1), normalization)
    if degree == 2:
        ratio = (1 - eta) / eta
        x = (ratio - 1) / (ratio + 1)
        # log(ratio) = 2 atanh(x). The nonnegative series gives an
        # independent exact enclosure and a bound on the omitted tail.
        terms = 40
        logarithm_low = 2 * sum(x**(2*k + 1) / (2*k + 1) for k in range(terms))
        logarithm_high = logarithm_low + 2 * x**(2 * terms + 1) / ((2 * terms + 1) * (1 - x*x))
        assert low <= logarithm_low <= logarithm_high <= high
        low, high = logarithm_low, logarithm_high
    return low, high


def linear_interval(bounds, coefficients):
    low = high = Fraction(0)
    for (left, right), scale in zip(bounds, coefficients):
        extremes = left * scale, right * scale
        low += min(extremes)
        high += max(extremes)
    return low, high


def compute_certificates(eta):
    constants = [factor_constant(eta, j) for j in range(1, 5)]
    results = {}
    for weight, label in ((Fraction(1, 2), "values"), (1, "lattice")):
        total = linear_interval(constants, [weight**i for i in range(4)])
        derivative = linear_interval(constants, [(1 - 2*i)*weight**i for i in range(4)])
        derivative_without_four = linear_interval(constants, [(1-2*i)*weight**i for i in range(3)])
        proper = linear_interval(constants, [1, *[-i*weight**i for i in range(1, 4)]])
        composite_fraction = (1 - 1/total[0], 1 - 1/total[1])
        proper_ratios = [value / mass for value in proper for mass in total]
        normalized_proper = min(proper_ratios), max(proper_ratios)
        results[label] = (derivative, derivative_without_four, proper,
                          composite_fraction, normalized_proper)
    return constants, results


def show_interval(bounds):
    # Display an enclosure rounded outwards to six decimal places.
    lower, upper = bounds
    units = 10**6
    floor = (units * lower.numerator) // lower.denominator
    ceiling = -((-units * upper.numerator) // upper.denominator)
    return f"[{floor / units:.6f}, {ceiling / units:.6f}]"


def arithmetic_samples():
    limit = 999_998
    least = least_prime_factors(limit)
    omega, squarefree = squarefree_data(least)
    split = bytearray(limit + 1)
    split[1] = 1
    for n in range(2, limit + 1):
        p = least[n]
        split[n] = int(p % 6 == 1 and split[n // p])
    for target in (998, 9_998, 99_998, limit):
        cutoff = integer_root(target, 5)
        cells = [[0]*5 for _ in range(5)]
        for a in range(7, target // 2, 6):
            b = target - a
            if (least[a] > cutoff and least[b] > cutoff and split[a] and split[b]
                    and squarefree[a] and squarefree[b] and gcd(a, b) == 1):
                cells[omega[a]][omega[b]] += 1
                cells[omega[b]][omega[a]] += 1
        unrestricted = coefficients(target, least, omega, squarefree, (1, 5))
        assert all(cells[i][j] <= unrestricted[i][j] for i in range(1, 5) for j in range(1, 5))
        print("sample", target, "prime pairs", cells[1][1])
        for k, label in ((1, "values"), (2, "lattice")):
            rows = [sum(cells[i][j]*k**(i+j-2) for j in range(1, 5)) for i in range(1, 5)]
            total = sum(rows)
            derivative = sum((1-2*i)*count for i, count in enumerate(rows))
            proper = rows[0] - sum(i*count for i, count in enumerate(rows[1:], 1))
            composite_fraction = Fraction(sum(rows[1:]), total)
            print(f"  {label}: derivative={float(Fraction(derivative, total)):+.5f}, "
                  f"proper={float(Fraction(proper, total)):+.5f}, "
                  f"composite={float(composite_fraction):.5f}")


def main():
    eta = Fraction(1, 5)
    constants, results = compute_certificates(eta)
    print("at", eta, "factor constants", *(show_interval(item) for item in constants))
    for label, bounds in results.items():
        print(label, "derivative, omit degree four, proper, composite, proper/total",
              *(show_interval(item) for item in bounds))
    value_derivative, omit_four, value_proper, value_composite, value_ratio = results["values"]
    lattice_derivative, _, lattice_proper, lattice_composite, _ = results["lattice"]
    assert value_derivative[1] < 0 < omit_four[0]
    assert value_proper[0] > 0 and value_composite[1] < Fraction(1, 2)
    assert Fraction(545, 10000) < value_ratio[0] <= value_ratio[1] < Fraction(55, 1000)
    assert lattice_derivative[1] < 0 and lattice_proper[1] < 0
    assert lattice_composite[0] > Fraction(1, 2)
    for endpoint, label, sign in ((Fraction(101, 500), "values", 1),
                                  (Fraction(29, 100), "lattice", -1),
                                  (Fraction(73, 250), "lattice", 1)):
        _, data = compute_certificates(endpoint)
        bounds = data[label][0]
        assert bounds[0] > 0 if sign > 0 else bounds[1] < 0
        print("derivative sign", endpoint, label, show_interval(bounds))
    arithmetic_samples()


if __name__ == "__main__":
    main()

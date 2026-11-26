#!/usr/bin/env python3
from fractions import Fraction
from itertools import combinations
from math import comb, factorial


def valuation(value):
    if not value:
        return float("inf")
    value = Fraction(value)
    n, d = abs(value.numerator), value.denominator
    return (n & -n).bit_length() - (d & -d).bit_length()


def falling(n, s):
    return factorial(n) // factorial(n - s) if n >= s else 0


def coefficient(j, d, z):
    return Fraction(factorial(j - d), factorial(j) * factorial(d)) / z**d


def check_derivatives():
    count = 0
    for z, c in ((Fraction(1, 8), Fraction(1, 2)), (Fraction(-3, 16), 1)):
        for j in range(3, 15):
            tail = [sum(comb(2 * k, k) * (c / z)**k for k in range(m + 1))
                    for m in range(j + 1)]
            for s in range(3):
                scale = falling(j, s) * z**s
                for d in range(s, j + 1):
                    left = coefficient(j - s, d - s, z)
                    assert left == scale * falling(d, s) * coefficient(j, d, z)
                    count += 1
                left = sum(
                    coefficient(j - s, d, z) * tail[j - s - d]
                    for d in range(j - s + 1)
                )
                right = scale * sum(
                    falling(d, s) * coefficient(j, d, z) * tail[j - d]
                    for d in range(j + 1)
                )
                assert left == right
                count += 1
    return count


def block(j, length, z, x):
    relative = [Fraction(1)]
    for b in range(1, length):
        relative.append(relative[-1] / (z * (j - x - b + 1) * (x + b)))
    return relative


def divided(first_x, first_y, second_x, second_y):
    return (first_y - second_y) / (first_x - second_x)


def check_blocks():
    valuation_count = 0
    interpolation_count = 0
    for tau, y, unit in ((4, 12, 1), (4, 36, -3), (5, 24, 5)):
        length = 2**tau
        j = length * (y + 1) - 1
        z = Fraction(unit, 2 ** (tau - 1))
        nodes = [length * q for q in (0, 1, 2, 3, 5, 8, 12)]
        weights = [block(j, length, z, x) for x in nodes]
        for weight in weights:
            for b, term in enumerate(weight):
                assert valuation(term) == (tau - 3) * b + 2 * b.bit_count()
                valuation_count += 1
        for i, m, n in combinations(range(len(nodes)), 3):
            for b in range(1, length):
                x_i, x_m, x_n = nodes[i], nodes[m], nodes[n]
                first = divided(x_i, weights[i][b], x_m, weights[m][b])
                second = divided(x_i, weights[i][b], x_n, weights[n][b])
                assert valuation(first) >= 2
                assert valuation((first - second) / (x_m - x_n)) >= 2
                valuation_count += 2
        linear = []
        quadratic = []
        for x, weight in zip(nodes, weights):
            total = sum(weight)
            first = sum(b * term for b, term in enumerate(weight)) / total
            second = sum(b * b * term for b, term in enumerate(weight)) / total
            p = x + first
            linear.append(p)
            quadratic.append(p * p - p + second - first * first)
        for i, m in combinations(range(len(nodes)), 2):
            assert valuation(linear[i] - linear[m]) == valuation(nodes[i] - nodes[m])
            interpolation_count += 1
            slope = divided(linear[i], quadratic[i], linear[m], quadratic[m])
            assert valuation(slope) >= 0
            for n in range(len(nodes)):
                if n in (i, m):
                    continue
                actual = quadratic[n] - quadratic[i] - slope * (linear[n] - linear[i])
                expected = valuation(nodes[n] - nodes[i]) + valuation(nodes[n] - nodes[m])
                assert valuation(actual) == expected
                interpolation_count += 1
    return valuation_count, interpolation_count


def borrow_score(weight, y, q):
    return weight * q + q.bit_count() + y.bit_count() - (y - q).bit_count()


def check_anchors():
    count = 0
    for tau in range(4, 10):
        weight = (tau - 3) * 2**tau
        y = 2 ** (weight + 2)
        bound = y.bit_length()
        base = [borrow_score(weight, y, q) for q in range(bound + 1)]
        assert [q for q, cost in enumerate(base) if cost == min(base)] == [0, 1]
        anchored = [(base[q] + valuation(q), q) for q in range(1, bound + 1)]
        best = min(cost for cost, _ in anchored)
        assert [q for cost, q in anchored if cost == best] == [1]
        count += 2
    return count


def main():
    derivative = check_derivatives()
    block, interpolation = check_blocks()
    anchors = check_anchors()
    print(
        f"{derivative} derivative checks, {block:,} block checks, "
        f"{interpolation} interpolation checks, {anchors} borrow checks"
    )


if __name__ == "__main__":
    main()

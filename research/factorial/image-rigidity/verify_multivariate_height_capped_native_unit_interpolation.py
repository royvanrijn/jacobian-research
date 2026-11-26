#!/usr/bin/env python3
"""Exact residue witnesses for height-capped native unit interpolation."""

from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, prod


def factorint(value):
    factors = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            value //= divisor
        divisor += 1
    if value > 1:
        factors[value] = 1
    return factors


def isprime(value):
    return value > 1 and factorint(value) == {value: 1}


def primitive_root(prime):
    factors = factorint(prime - 1)
    return next(candidate for candidate in range(2, prime) if all(
        pow(candidate, (prime - 1) // factor, prime) != 1 for factor in factors
    ))


def crt(moduli, values):
    modulus = prod(moduli)
    return sum(
        value * (modulus // factor) * pow(modulus // factor, -1, factor)
        for factor, value in zip(moduli, values, strict=True)
    ) % modulus


def unit_at(value, prime, modulus):
    valuation = 0
    while value % prime == 0:
        valuation += 1
        value //= prime
    return valuation, value % modulus


def coefficient_data(rows, order, prime, modulus):
    """Stripped coefficient valuations and units, with the index reversed."""
    result = []
    for index in range(order + 1):
        valuation, unit = unit_at(comb(order, index), prime, modulus)
        for (left, right), count in rows.items():
            base = order * min(left, right)
            length = abs(left - right) * (index if left > right else order - index)
            divisor, row_unit = unit_at(prod(range(base + 1, base + length + 1)), prime, modulus)
            valuation += divisor * count
            unit = unit * pow(row_unit, count, modulus) % modulus
        result.append((valuation, unit))
    return result


def graph_labels(ratios, modulus):
    graph = defaultdict(list)
    for (order, index), ratio in ratios.items():
        left = 0, Fraction(index, order)
        right = 1, Fraction(order - index + 1, order)
        graph[left].append((right, pow(ratio, -1, modulus)))
        graph[right].append((left, ratio))
    labels = {}
    for vertex in graph:
        if vertex in labels:
            continue
        labels[vertex] = 1
        pending = [vertex]
        while pending:
            source = pending.pop()
            for target, multiplier in graph[source]:
                value = labels[source] * multiplier % modulus
                if target in labels:
                    assert labels[target] == value
                else:
                    labels[target] = value
                    pending.append(target)
    return labels


def solve_with_unit_pivots(matrix, target, modulus):
    rows = [row + [value] for row, value in zip(matrix, target, strict=True)]
    width = len(matrix[0])
    pivots = []
    prime = next(iter(factorint(modulus)))
    for top in range(len(rows)):
        column = next(
            k for k in range(width)
            if any(value[k] % prime for value in rows[top:])
        )
        pivot = next(
            k for k in range(top, len(rows))
            if rows[k][column] % prime
        )
        rows[top], rows[pivot] = rows[pivot], rows[top]
        row = rows[top]
        inverse = pow(row[column], -1, modulus)
        rows[top] = [entry * inverse % modulus for entry in row]
        for k, other in enumerate(rows):
            if k != top:
                rows[k] = [
                    (old - other[column] * entry) % modulus
                    for old, entry in zip(other, rows[top], strict=True)
                ]
        pivots.append(column)
    solution = [0] * width
    for column, row in zip(pivots, rows, strict=True):
        solution[column] = row[-1]
    assert all(
        sum(x * y for x, y in zip(row, solution, strict=True)) % modulus == value % modulus
        for row, value in zip(matrix, target, strict=True)
    )
    return solution


def realize(labels, order, height, prime):
    fractions = sorted({fraction for _, fraction in labels})
    roots = [f.numerator * pow(f.denominator, -1, prime) % prime for f in fractions]
    bases = range((height + 3) // 4, height // 3 + 1)
    generator = primitive_root(prime)
    logs = {pow(generator, exponent, prime): exponent for exponent in range(prime - 1)}
    matrix = [[logs[(base + root) % prime] for base in bases] for root in roots]
    matrix.append([1] * len(bases))
    moduli = [factor**power for factor, power in factorint(prime - 1).items()]
    counts = []
    for side in range(2):
        target = [logs[labels[side, fraction]] for fraction in fractions] + [0]
        solutions = [solve_with_unit_pivots(matrix, target, modulus) for modulus in moduli]
        counts.append([crt(moduli, column) for column in zip(*solutions, strict=True)])
    smaller, larger = sorted(counts, key=sum)
    difference = sum(larger) - sum(smaller)
    assert difference % (prime - 1) == 0
    smaller[0] += difference
    assert sum(counts[0]) == sum(counts[1])
    for side, vector in enumerate(counts):
        for fraction, root in zip(fractions, roots, strict=True):
            value = prod(pow(base + root, count, prime) for base, count in zip(bases, vector))
            assert value % prime == labels[side, fraction]
    rows = Counter({(base + 1, base): count for base, count in zip(bases, counts[0]) if count})
    rows.update({(base, base + 1): count for base, count in zip(bases, counts[1]) if count})
    for (left, right) in rows:
        for j in (order, order + 1):
            assert all((j * min(left, right) + k) % prime for k in range(1, j + 1))
    return rows



def realize_prime_power(labels, order, height, prime, depth):
    if depth == 1:
        return realize(labels, order, height, prime)
    modulus = prime**depth
    fractions = sorted({fraction for _, fraction in labels})
    roots = [f.numerator * pow(f.denominator, -1, modulus) % modulus for f in fractions]
    markers = [a for a in range(prime) if all((a + root) % prime for root in roots)]
    markers = markers[:len(fractions)]
    copies = 15 * order
    low, high = (height + 3) // 4, height // 3
    marker_product = [
        prod(pow(base + root, copies, modulus) for base in markers) % modulus
        for root in roots
    ]
    adjusted = {
        (side, fraction): labels[side, fraction] * pow(group, -1, prime) % prime
        for side in range(2)
        for fraction, group in zip(fractions, marker_product, strict=True)
    }
    seed_height = (prime - 1) // order
    seed = realize(adjusted, order, seed_height, prime)
    offset = prime * ((low - (seed_height + 3) // 4 + prime - 1) // prime)
    rows = Counter({(a + offset, b + offset): count for (a, b), count in seed.items()})
    for side in range(2):
        fixed_bases = {min(a, b): count for (a, b), count in rows.items() if (a < b) == side}
        target = [
            labels[side, fraction]
            * pow(prod(pow(base + root, count, modulus)
                       for base, count in fixed_bases.items()), -1, modulus) % modulus
            for fraction, root in zip(fractions, roots, strict=True)
        ]
        cauchy = [[pow((base + root) % prime, -1, prime) for base in markers] for root in roots]
        bases = markers.copy()
        for exponent in range(1, depth):
            step = prime**exponent
            residual = []
            for value, root in zip(target, roots, strict=True):
                current = prod(pow(base + root, copies, modulus) for base in bases) % modulus
                ratio = value * pow(current, -1, modulus) % modulus
                assert ratio % step == 1
                residual.append((ratio - 1) // step % prime)
            if exponent < depth - 1:
                repeated = [[copies * value for value in row] for row in cauchy]
                digits = solve_with_unit_pivots(repeated, residual, prime)
                bases = [base + step * digit for base, digit in zip(bases, digits, strict=True)]
                continue
            digit_sums = solve_with_unit_pivots(cauchy, residual, prime)
            for base, digit_sum in zip(bases, digit_sums, strict=True):
                first = (low - base + step - 1) // step
                last = (high - base) // step
                remaining = (digit_sum - copies * first) % prime
                assert remaining <= copies * (last - first)
                for _ in range(copies):
                    increment = min(remaining, last - first)
                    remaining -= increment
                    actual = base + step * (first + increment)
                    row = (actual, actual + 1) if side else (actual + 1, actual)
                    rows[row] += 1
                assert remaining == 0
    assert all(low <= min(row) <= high for row in rows)
    assert sum(count for (a, b), count in rows.items() if a < b) * 2 == sum(rows.values())
    for (left, right) in rows:
        for j in (order, order + 1):
            assert all((j * min(left, right) + k) % prime for k in range(1, j + 1))
    return rows


def check(order, prime, depth=1):
    modulus = prime**depth
    height = (modulus - 1) // order
    assert isprime(prime) and order * height + 1 == modulus
    rows = Counter([(2, 5), (4, 1), (13, 14), (height - 1, height),
                    (height - 1 + order // 2, height - 1)])
    original = {j: coefficient_data(rows, j, prime, modulus) for j in (order, order + 1)}
    assert [v for v, _ in original[order]] == [0, 0, 0] + [depth] * (order - 2)
    assert [v for v, _ in original[order + 1]] == [0] * (order + 2)
    for compatible in [True, False]:
        if compatible:
            binomial = [comb(order - 1, k) for k in range(order)] + [0, 0]
            odd = [sum(c * binomial[k - i] for i, c in enumerate([1, 3, 2]) if k >= i)
                   for k in range(order + 2)]
        else:
            odd = [comb(order + 1, k) * 3**k for k in range(order + 2)]
        even = [1, 3, 2] + [0] * (order - 2)
        targets = {order: even, order + 1: odd}
        ratios = {}
        for j, coefficients in targets.items():
            for k in range(1, j + 1):
                old = original[j]
                ratios[j, k] = (
                    coefficients[k] * pow(coefficients[k - 1], -1, modulus)
                    * old[k - 1][1] * pow(old[k][1], -1, modulus) % modulus
                    if coefficients[k] else 1
                )
        added = realize_prime_power(graph_labels(ratios, modulus), order, height, prime, depth)
        augmented = rows + added
        for j, target in targets.items():
            exact = coefficient_data(augmented, j, prime, modulus)
            assert [v for v, _ in exact] == [v for v, _ in original[j]]
            assert [u * pow(exact[0][1], -1, modulus) % modulus if v == 0 else 0
                    for v, u in exact] == [c % modulus for c in target]
        for root in (-1, pow(-2, -1, modulus)):
            residue = sum(c * pow(root, k, modulus) for k, c in enumerate(odd)) % modulus
            assert (residue == 0) == compatible
        print(f"order {order}, native modulus {prime}^{depth}, {sum(added.values())} added unit rows, "
              f"quadratic {'divides' if compatible else 'coprime to'} the odd reduction")


if __name__ == "__main__":
    for order, height in [(6, 400), (8, 450), (10, 420)]:
        while not isprime(order * height + 1):
            height += 1
        check(order, order * height + 1)

    for depth in (2, 3):
        check(6, 2437, depth)

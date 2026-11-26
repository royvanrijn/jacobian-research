#!/usr/bin/env python3
"""Replay full linked-anchor residue calibration over a product modulus."""

from collections import Counter
from fractions import Fraction
from math import comb, gcd, lcm, prod

from verify_multivariate_height_capped_native_unit_interpolation import (
    coefficient_data,
    crt,
    factorint,
    graph_labels,
    primitive_root,
    solve_with_unit_pivots,
)


def anchor_primes(height):
    anchors = {}
    for order, h_limit in ((6, 3), (7, 2)):
        for h in range(1, h_limit + 1):
            for prime, exponent in factorint((order * height + h) // gcd(order, h)).items():
                if prime < 151:
                    assert (order, h) not in ((6, 1), (6, 2))
                    continue
                assert exponent == 1 and prime not in anchors
                anchors[prime] = order, h
    return anchors


def solve_mixed_logs(matrix, target, equation_moduli):
    total_modulus = lcm(*equation_moduli)
    width = len(matrix[0])
    solutions, moduli = [], []
    for prime, exponent in factorint(total_modulus).items():
        chosen = [0] * width
        digit_place = 1
        for _ in range(exponent):
            active_modulus = digit_place * prime
            equations = [
                (row, value) for row, value, modulus in
                zip(matrix, target, equation_moduli, strict=True) if modulus % active_modulus == 0
            ]
            active_rows = [row for row, _ in equations]
            discrepancies = [
                value - sum(x * y for x, y in zip(row, chosen, strict=True))
                for row, value in equations
            ]
            assert all(discrepancy % digit_place == 0 for discrepancy in discrepancies)
            digits = solve_with_unit_pivots(
                active_rows, [error // digit_place for error in discrepancies], prime
            )
            chosen = [value + digit_place * digit for value, digit in zip(chosen, digits, strict=True)]
            digit_place = active_modulus
        solutions.append(chosen)
        moduli.append(digit_place)
    return [crt(moduli, column) for column in zip(*solutions, strict=True)]


def complementary_target(marked, compatible):
    order = 13 - marked
    if not compatible:
        return [comb(order, k) * 3**k for k in range(order + 1)]
    binomial = [comb(order - 2, k) for k in range(order - 1)] + [0, 0]
    return [
        sum(c * binomial[k - i] for i, c in enumerate((1, 3, 2)) if k >= i)
        for k in range(order + 1)
    ]


def check(height, anchors, compatible):
    original_rows = Counter([(2, 5), (4, 1), (13, 14), (height - 1, height), (height + 2, height - 1)])
    fractions = sorted({Fraction(k, j) for j in (6, 7) for k in range(1, j + 1)})
    roots = {
        prime: [r.numerator * pow(r.denominator, -1, prime) % prime for r in fractions]
        for prime in anchors
    }
    low, high = (height + 3) // 4, height // 3
    bases = [
        a for a in range(low, high + 1)
        if all(all((a + root) % prime for root in root_list)
               for prime, root_list in roots.items())
    ]
    matrix, equation_moduli = [], []
    target_logs = [[], []]
    target_arrays = {}
    for prime, (marked, _) in anchors.items():
        unmarked = 13 - marked
        coefficient_targets = {
            marked: [1, 3, 2] + [0] * (marked - 2),
            unmarked: complementary_target(marked, compatible),
        }
        target_arrays[prime] = coefficient_targets
        ratios = {}
        for j, coefficients in coefficient_targets.items():
            original = coefficient_data(original_rows, j, prime, prime)
            assert [v == 0 for v, _ in original] == [c != 0 for c in coefficients]
            for k in range(1, j + 1):
                ratios[j, k] = (
                    coefficients[k] * pow(coefficients[k - 1], -1, prime)
                    * original[k - 1][1] * pow(original[k][1], -1, prime) % prime
                    if coefficients[k] else 1
                )
        labels = graph_labels(ratios, prime)
        generator = primitive_root(prime)
        log = {pow(generator, exponent, prime): exponent for exponent in range(prime - 1)}
        for f, root in zip(fractions, roots[prime], strict=True):
            matrix.append([log[(a + root) % prime] for a in bases])
            equation_moduli.append(prime - 1)
            for side in range(2):
                target_logs[side].append(log[labels[side, f]])
    total_modulus = lcm(*equation_moduli)
    matrix.append([1] * len(bases))
    equation_moduli.append(total_modulus)
    for side in range(2):
        target_logs[side].append(0)
    counts = [solve_mixed_logs(matrix, target, equation_moduli) for target in target_logs]
    smaller, larger = sorted(counts, key=sum)
    difference = sum(larger) - sum(smaller)
    assert difference % total_modulus == 0
    smaller[0] += difference
    assert sum(counts[0]) == sum(counts[1])
    added = Counter({
        (a + 1, a) if side == 0 else (a, a + 1): count
        for side, vector in enumerate(counts)
        for a, count in zip(bases, vector, strict=True) if count
    })
    assert all(low <= min(row) <= high for row in added)
    for prime, target in target_arrays.items():
        marked, _ = anchors[prime]
        for j, coefficients in target.items():
            before = coefficient_data(original_rows, j, prime, prime)
            after = coefficient_data(original_rows + added, j, prime, prime)
            assert [v for v, _ in before] == [v for v, _ in after]
            assert [unit * pow(after[0][1], -1, prime) % prime if v == 0 else 0
                    for v, unit in after] == [c % prime for c in coefficients]
        for root in (-1, pow(-2, -1, prime)):
            residue = sum(c * pow(root, k, prime) for k, c in enumerate(target[13 - marked]))
            assert (residue % prime == 0) == compatible
    print(f"B={height}, selected moduli {6 * height + 1} and {3 * height + 1}, "
          f"{len(anchors)} joint primes with {len(bases)} allowed bases, "
          f"{'compatible' if compatible else 'coprime'} reductions, "
          f"{sum(added.values())} appended rows")


if __name__ == "__main__":
    height = 8656
    anchors = anchor_primes(height)
    assert anchors == {167: (6, 1), 311: (6, 1), 25969: (6, 2), 199: (6, 3), 10099: (7, 2)}
    assert prod(prime for prime, anchor in anchors.items() if anchor == (6, 1)) == 6 * height + 1
    assert prod(prime for prime, anchor in anchors.items() if anchor == (6, 2)) == 3 * height + 1
    for compatible in (True, False):
        check(height, anchors, compatible)

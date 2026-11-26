#!/usr/bin/env python3
from fractions import Fraction
from math import comb, factorial, gcd, lcm


def carries(left, right):
    return left.bit_count() + right.bit_count() - (left + right).bit_count()


def valuation(integer):
    return (integer & -integer).bit_length() - 1


def least_levels(order, r, h):
    scale = lcm(r.denominator, h.denominator)
    axis = int(scale * r)
    penalty = int(scale * h)
    levels = {}
    for total in range(order + 1):
        eta = total * axis + scale * (
            (order - total).bit_count() - order.bit_count()
        )
        for cross in range(total + 1):
            half = cross // 2
            remaining = total - cross
            carry_min = min(
                carries(right, half) + carries(remaining - right, half)
                for right in range(remaining + 1)
            )
            value = eta + cross * penalty + scale * (
                carry_min - cross.bit_count()
            )
            levels.setdefault(value, set()).add((total, cross))
    return [
        (Fraction(value, scale), levels[value]) for value in sorted(levels)[:2]
    ]


def residue_exponents(total, cross):
    half = cross // 2
    remaining = total - cross
    return [
        right
        for right in range(remaining + 1)
        if not (right & half or (remaining - right) & half)
    ]


def parity(exponents, period):
    result = set()
    for exponent in exponents:
        result.symmetric_difference_update({exponent % period})
    return result


def geometric_identity(total, cross, step, length, period, shift=()):
    exponents = residue_exponents(total, cross)
    left = parity(exponents + [x + step for x in exponents], period)
    base = [0, step * length]
    right = parity(
        base + [x + extra for x in base for extra in shift], period
    )
    assert left == right


def check_coefficients():
    count = 0
    small_factorial = [factorial(i) for i in range(21)]
    for order in range(1, 11):
        for total in range(order + 1):
            for cross in range(total + 1):
                for right in range(total - cross + 1):
                    left = total - cross - right
                    coefficient = (
                        comb(order, total)
                        * comb(total, cross)
                        * comb(total - cross, right)
                        * small_factorial[2 * left + cross]
                        * small_factorial[2 * right + cross]
                    )
                    expected = (
                        2 * total - cross.bit_count()
                        + (order - total).bit_count() - order.bit_count()
                        + carries(left, cross // 2)
                        + carries(right, cross // 2)
                    )
                    assert valuation(coefficient) == expected
                    count += 1
    return count


def check_dyadic_terms():
    axis_checks = 0
    pair_checks = 0
    for e in range(3, 6):
        d = 2**e
        for high in (1, 3, 5):
            terminal = d * high
            for h in (
                Fraction(3, 4),
                Fraction(1, 2),
                Fraction(3, 8),
                Fraction(1, 4),
                Fraction(1, d),
            ):
                r = (2 - h) / d
                values = [
                    i * r + (terminal - i).bit_count()
                    - terminal.bit_count()
                    for i in range(1, terminal + 1)
                ]
                assert min(values) == d * r - 1
                assert values.count(min(values)) == 1
                assert values[d - 1] == min(values)
                axis_checks += 1
            fixtures = [
                (Fraction(3, 4), terminal - 3, {(d - 3, 1)}),
                (Fraction(1, 2), terminal - 2, {(d - 2, 1)}),
                (Fraction(3, 8), terminal - 1, {(d - 1, 3)}),
                (Fraction(1, 4), terminal - 3, {(d - 3, 3)}),
                (Fraction(1, 4), terminal - 2, {(d - 2, 1), (d - 2, 5)}),
                (Fraction(1, d), terminal - 1, {(d - 1, d - 1)}),
            ]
            for t in range(2, e + 1):
                fixtures.append(
                    (Fraction(3, 2**(t + 1)), terminal - 1, {(d - 1, 2**t - 1)})
                )
            for t in range(2, e):
                fixtures.append(
                    (Fraction(1, 2**t), terminal - 1,
                     {(d - 1, 2**t - 1), (d - 1, 2**(t + 1) - 1)})
                )
            for h, order, expected in fixtures:
                r = (2 - h) / d
                assert least_levels(order, r, h)[0][1] == expected, (d, order, h)
                pair_checks += 1
    return axis_checks, pair_checks


def check_residues():
    count = 0
    for e in range(3, 11):
        d = 2**e
        period = d - 1
        geometric_identity(d - 3, 1, 1, d - 3, period)
        geometric_identity(d - 2, 1, 1, d - 2, period)
        geometric_identity(d - 1, 3, 2, (d - 2) // 2, period)
        geometric_identity(d - 3, 3, 2, (d - 4) // 2, period)
        geometric_identity(d - 2, 5, 4, (d - 4) // 4, period, (1,))
        count += 5
        if period % 3 == 0:
            assert not parity(residue_exponents(d - 2, 5), 3)
            count += 1
        for t in range(2, e + 1):
            step = 2**(t - 1)
            geometric_identity(d - 1, 2**t - 1, step, d // step - 1, period)
            assert gcd(period, step - 1) == 2**gcd(e, t - 1) - 1
            count += 2
    return count


def check_subfield_layer():
    count = 0
    for e, t in ((4, 3), (5, 3), (6, 3), (6, 4), (6, 5)):
        d = 2**e
        step = 2**(t - 1)
        add_tie = Fraction(2, 4 * step - 1)
        drop_tie = Fraction(2, 2 * step + 1)
        strengths = (
            (Fraction(1, 2 * step) + add_tie) / 2,
            add_tie,
            (add_tie + drop_tie) / 2,
            drop_tie,
            (drop_tie + Fraction(1, step)) / 2,
        )
        base_cross = 2**t - 1
        for h in strengths:
            r = (2 - h) / d
            candidates = {
                (d // 2 - 1, base_cross): h / 2,
                (d - 1, base_cross - step): 1 - step * h,
                (d - 1, base_cross + 2 * step): 2 * step * h - 1,
            }
            gap = min(candidates.values())
            expected = {pair for pair, cost in candidates.items() if cost == gap}
            first, second = least_levels(d - 1, r, h)
            assert first[1] == {(d - 1, base_cross)}
            assert second == (first[0] + gap, expected)
            count += 1
        geometric_identity(
            d // 2 - 1, base_cross, step, d // (2 * step) - 1, d - 1
        )
        residue_order = gcd(d - 1, step - 1)
        for exponent in (d // 2 - step, d - step // 2, d - 2 * step):
            assert gcd(residue_order, exponent) == 1
            count += 1
    return count


def main():
    coefficients = check_coefficients()
    axes, pairs = check_dyadic_terms()
    residues = check_residues()
    subfield = check_subfield_layer()
    print(
        f"{coefficients:,} coefficient valuations, {axes} exact-bit axis minima, "
        f"{pairs} least-pair fixtures, {residues} residue checks, "
        f"{subfield} next-layer checks"
    )


if __name__ == "__main__":
    main()

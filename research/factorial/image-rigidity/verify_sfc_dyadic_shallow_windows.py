#!/usr/bin/env python3
"""Replay the finite digit and derivative tables in SFC_DYADIC_SHALLOW_WINDOWS.md."""

from fractions import Fraction


def intercepts(residue, high_valuation, depth):
    modulus = 2**depth
    return tuple(
        m.bit_count() + residue.bit_count() - (residue - m).bit_count()
        if m <= residue
        else m.bit_count() + residue.bit_count()
        - (modulus + residue - m).bit_count() - high_valuation + 1
        for m in range(modulus)
    )


def adjacent_classes(depth):
    modulus = 2**depth
    cap = 2 * depth
    for residue in range(modulus - 1):
        for high_valuation in range(cap + 1):
            yield residue, high_valuation, residue + 1, high_valuation
    for positive_valuation in range(1, cap + 1):
        yield modulus - 1, 0, 0, positive_valuation
        yield modulus - 1, positive_valuation, 0, 0


def maximizers(row, gap):
    peak = max(c - gap * m for m, c in enumerate(row))
    return tuple(m for m, c in enumerate(row) if c - gap * m == peak)


def intersections(row):
    return {
        Fraction(c - prior, m - k)
        for m, c in enumerate(row)
        for k, prior in enumerate(row[:m])
    }


def residual_mask(indices):
    return sum(1 << (m - indices[0]) for m in indices)


def binary_polynomial_gcd(a, b):
    while b:
        while a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def residual_failures(depth, lower_gap, upper_gap):
    result = set()
    for r, s, following_r, following_s in adjacent_classes(depth):
        a = intercepts(r, s, depth)
        b = intercepts(following_r, following_s, depth)
        for gap in intersections(a) & intersections(b):
            if not lower_gap < gap <= upper_gap:
                continue
            left, right = maximizers(a, gap), maximizers(b, gap)
            if binary_polynomial_gcd(residual_mask(left), residual_mask(right)) > 1:
                result.add((depth, 2 - gap, following_r, following_s))
    return result


def ordinary_layers(weights, precision):
    base = min(weights.values())
    layers = tuple(
        tuple(m for m, value in weights.items() if value == base + k)
        for k in range(precision)
    )
    return base, layers


def valuation_2(n):
    return (abs(n) & -abs(n)).bit_length() - 1


def check_derivative_table(failures):
    # Each pair of tuples lists the exponents in layers 0, 1, ... .
    expected = {
        (Fraction(3, 2), 7): (((3, 7), ()), ((2, 6), (1, 5, 7))),
        (Fraction(7, 4), 0): (((7, 15),), ((0,),)),
        (Fraction(7, 4), 11): (((3, 7, 11), ()), ((2, 6, 10), (7,))),
        (Fraction(7, 4), 12): (((3, 7, 11),), ((4, 12),)),
        (Fraction(7, 4), 14): (((6, 14), (), ()), ((5, 13), (), (3, 7, 11))),
        (Fraction(7, 4), 15): (
            ((7, 15), (), (), (), (3, 11), ()),
            ((6, 14), (), (), (5, 13), (2, 10), (3, 11, 15)),
        ),
        (Fraction(11, 6), 3): (((3, 15), ()), ((2, 14), (15,))),
        (Fraction(11, 6), 4): (((3, 15),), ((4,),)),
    }
    found = set()
    for depth, beta, residue, high_valuation in failures:
        found.add((beta, residue))
        gap = 2 - beta
        q = beta.denominator
        row = intercepts(residue, high_valuation, depth)
        expected_moment, expected_derivative = expected[beta, residue]
        weights = {m: q * (gap * m - c) for m, c in enumerate(row)}
        minimum, layers = ordinary_layers(weights, len(expected_moment))
        assert layers == expected_moment
        g = layers[0][0]
        derivative_weights = {m: value + q * valuation_2(m - g)
                              for m, value in weights.items() if m != g}
        derivative_minimum, derivative = ordinary_layers(
            derivative_weights, len(expected_derivative)
        )
        assert derivative == expected_derivative
        if high_valuation == 2 * depth:
            assert max(m for layer in layers + derivative for m in layer) <= residue

        # All higher blocks are separated from the layers used here.
        tail_margin = q * (gap * 2**depth - 2)
        assert tail_margin >= len(expected_moment)
        assert tail_margin >= len(expected_derivative)
        assert minimum + tail_margin + q * depth >= (
            derivative_minimum + len(expected_derivative)
        )
    assert found == set(expected)


def layer_supports(depth, beta, residue, high_valuation, count_s, count_d, shift=0):
    q = beta.denominator
    gap = 2 - beta
    row = intercepts(residue, high_valuation, depth)
    weights = {m: q * (gap * m - c) for m, c in enumerate(row)}
    minimum, layers = ordinary_layers(weights, count_s)
    g = layers[0][0]
    derivative_weights = {
        m: value + q * valuation_2(m - g)
        for m, value in weights.items() if m != g
    }
    derivative_minimum, derivative = ordinary_layers(derivative_weights, count_d)
    if high_valuation == 2 * depth:
        assert max(m for layer in layers + derivative for m in layer) <= residue
    margin = q * (gap * 2**depth - 2)
    assert margin >= max(count_s, count_d)
    assert minimum + margin + q * depth >= derivative_minimum + count_d
    return (
        {k: layer for k, layer in enumerate(layers) if layer},
        {k: tuple(m + shift for m in layer) for k, layer in enumerate(derivative) if layer},
    )


def check_larger_slopes(failures):
    beta15 = Fraction(15, 8)
    beta53 = Fraction(53, 28)
    beta23 = Fraction(23, 12)
    beta27 = Fraction(27, 14)
    expected_residues = {
        beta15: {0, 1, 22, 23, 24, 25, 27, 28, 29, 30, 31},
        beta53: {3, 4},
        beta23: {6, 7, 8, 9},
        beta27: {3, 4},
    }
    found = {beta: set() for beta in expected_residues}
    for depth, beta, r, s in failures:
        found[beta].add(r)
        if beta == beta15:
            if r == 0:
                if s == 0:
                    args = (3, 3, 1)
                    expected = ({0: (15, 31)}, {0: (15, 31), 2: (1,)})
                else:
                    args = (1, 1, 0)
                    expected = ({0: (15, 31)}, {0: (0,)})
            elif r == 1:
                args = (1, 1, 0)
                expected = ({0: (15, 31)}, {0: (1,)})
            elif r == 22:
                args = (3, 3, 1)
                expected = ({0: (6, 14, 22)}, {0: (6, 14, 22), 2: (16,)})
            elif r == 23:
                args = (10, 10, 1)
                extra = {8: (31,)} if s == 0 else {}
                expected = ({0: (7, 15, 23)} | extra,
                            {0: (7, 15, 23), 7: (6, 14, 22), 9: (16,)} | extra)
            elif r == 24:
                args = (3, 3, 1)
                expected = ({0: (7, 15, 23)}, {0: (7, 15, 23), 2: (9, 25)})
            elif r == 25:
                args = (1, 1, 0)
                expected = ({0: (7, 15, 23)}, {0: (9, 25)})
            elif r == 27:
                args = (6, 6, 1)
                expected = ({0: (11, 27), 4: (7, 15, 23)},
                            {0: (11, 27), 4: (7, 15, 23), 5: (8, 16, 24)})
            elif r == 28:
                args = (5, 6, -1)
                expected = ({0: (11, 27), 1: (12, 28), 4: (7, 15, 23)}, {0: (11, 27)})
            elif r == 29:
                args = (6, 5, 2)
                expected = ({0: (13, 29)}, {0: (13, 29), 1: (14, 30), 4: (9, 17, 25)})
            elif r == 30:
                args = (11, 11, 1)
                expected = ({0: (14, 30), 7: (13, 29), 8: (6, 22)},
                            {0: (14, 30), 6: (12, 28), 8: (6, 22), 10: (8, 16, 24)})
            else:
                assert r == 31
                for lifted in (31, 63):
                    for upper_s in range(13):
                        first = {0: (15, 31), 8: (7, 23), 12: (11, 27),
                                 14: (13, 29), 15: (14, 30)}
                        second = {0: (15, 31), 7: (14, 30), 8: (7, 23),
                                  12: (11, 27), 13: (12, 28), 14: (13, 29),
                                  15: (6, 22), 17: (8, 24, 32)}
                        if lifted == 63:
                            first[16] = second[16] = (47, 63)
                        assert layer_supports(6, beta, lifted, upper_s, 18, 18, 1) == (first, second)
                continue
        elif beta in (beta53, beta27):
            if r == 3:
                k = 3 if beta == beta53 else 1
                args = (k + 1, k + 1, 1)
                expected = ({0: (3, 31)}, {0: (3, 31), k: (32,)})
            else:
                assert r == 4
                args = (1, 1, 0)
                expected = ({0: (3, 31)}, {0: (4,)})
        else:
            assert beta == beta23
            if r == 6:
                args = (3, 3, 1)
                expected = ({0: (6, 30)}, {0: (6, 30), 2: (32,)})
            elif r == 8:
                args = (3, 3, 1)
                expected = ({0: (7, 31)}, {0: (7, 31), 2: (9,)})
            elif r == 9:
                args = (1, 1, 0)
                expected = ({0: (7, 31)}, {0: (9,)})
            else:
                assert r == 7 and s == 0
                for upper_s in range(13):
                    shared = {0: (7, 31), 4: (23,), 8: (15, 39)}
                    assert layer_supports(6, beta, 39, upper_s, 14, 14, 1) == (
                        shared, shared | {11: (6, 30), 13: (32,)}
                    )
                continue
        actual = layer_supports(depth, beta, r, s, *args)
        assert actual == expected, (beta, r, s, actual, expected)
    assert found == expected_residues


def main():
    failures = residual_failures(3, Fraction(1, 4), Fraction(1))
    failures |= residual_failures(4, Fraction(1, 8), Fraction(1, 4))
    check_derivative_table(failures)
    check_larger_slopes(residual_failures(5, Fraction(1, 16), Fraction(1, 8)))
    print("PASS: seven shallow dyadic slopes and their derivative layer types")


if __name__ == "__main__":
    main()

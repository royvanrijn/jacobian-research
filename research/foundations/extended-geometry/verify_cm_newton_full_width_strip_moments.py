#!/usr/bin/env python3
"""Replay full-strip residues and the conductor-18 reconstruction over Q."""

from collections import defaultdict
from fractions import Fraction
from math import comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_cm_newton_two_input_strip_moments import normalized_moment_residue, power


def build_phase(n, root_multiplicity, prime, lower_terms):
    angular = power({(2, 0): 1, (0, 0): -1}, root_multiplicity, prime)
    phase = {
        (index - 1, n - 2 * (index - 1)): coefficient
        for (index, _), coefficient in angular.items()
    }
    phase[-1, n] = -angular[0, 0]
    phase.update(lower_terms)
    for angle, height in phase:
        assert height >= 0 and height + 2 * angle <= n
        assert (n - 1) * height + 2 * n * angle >= n * (n - 3)
    return phase


def moment_pair(phase, order, prime):
    return tuple(
        normalized_moment_residue(power(phase, exponent, prime), prime)
        for exponent in (order, order + prime)
    )


def check_wider_strip():
    n, a, prime, order = 10, 5, 17, 5
    lower_terms = {(0, 8): 4, (4, 1): 9}
    for angle, height in lower_terms:
        assert (a + 1) * height + 2 * a * angle < (a + 1) * n - 2 * a
    phase = build_phase(n, 3, prime, lower_terms)
    assert moment_pair(phase, order, prime) == (0, 16)
    print("Full-width control: residues (0, 16) modulo 17")


def check_upper_half_pole():
    n, s, prime = 18, 13, 43
    order = (s * prime - 1) // n
    phase = build_phase(n, 4, prime, {(0, n - 1): 6})
    powers = power(phase, order, prime)
    top_height = s * prime - 1
    inputs = (
        powers.get((0, top_height), 0),
        powers.get((prime, top_height - 2 * prime), 0),
        powers.get((0, top_height - prime), 0),
    )
    assert inputs == (0, 5, 13)
    assert moment_pair(phase, order, prime) == (0, 0)
    print("Upper-half correction: inputs (0, 5, 13), residues (0, 0) modulo 43")

    extra_residues = []
    for exponent in (order + 2 * prime, order + 3 * prime):
        # Every possible angular-constant denominator is in this interval,
        # independently of which coefficients vanish modulo prime.
        assert prime**2 < (n - 2) * exponent + 1 <= n * exponent + 1 < 2 * prime**2
        literal_power = power(phase, exponent, prime)
        residue = sum(
            coefficient * pow((height + 1) // prime, -1, prime)
            for (angle, height), coefficient in literal_power.items()
            if angle == 0 and height % 2 == 0 and (height + 1) % prime == 0
        )
        extra_residues.append(residue % prime)
    assert extra_residues == [0, 16]
    print("Upper-half continuation: orders 117 and 160 give (0, 16) modulo 43")

    prime = 61
    assert s * prime % (2 * n) == 1
    order = (s * prime - 1) // n
    phase = build_phase(n, 4, prime, {(0, n - 1): 6})
    assert normalized_moment_residue(power(phase, order, prime), prime) == 5
    print("Refined progression: order 44 gives 5 modulo 61")


def reconstruction_rows():
    # Monomials in z, t, u, v, h, after b, c, d become u, v, u*h.
    phase = {
        (2 * index - 1, 18, 0, 0, 0): (-1) ** index * comb(4, index)
        for index in range(5)
    }
    phase[-1, 16, 0, 0, 0] = -1
    phase[0, 17, 1, 0, 0] = 1
    phase[0, 16, 0, 1, 0] = 1
    phase[2, 17, 1, 0, 1] = 1
    coherent_weights = {-1: 1, -3: -2, -5: 1}
    current = {(0, 0, 0, 0, 0): 1}
    rows = []
    for _ in range(7):
        product = defaultdict(int)
        for left, left_coefficient in current.items():
            for right, right_coefficient in phase.items():
                exponent = tuple(a + b for a, b in zip(left, right))
                product[exponent] += left_coefficient * right_coefficient
        current = {key: value for key, value in product.items() if value}

        # Columns in (r, v, h), where r=u^2. Absorb one u in each C column.
        coherent, first_pole, second_pole = (defaultdict(Fraction) for _ in range(3))
        for (angle, height, u_degree, v_degree, h_degree), value in current.items():
            assert height % 2 == u_degree % 2
            exponent = (u_degree // 2, v_degree, h_degree)
            if height % 2 == 0 and angle in coherent_weights:
                coherent[exponent] += Fraction(value * coherent_weights[angle], 13 + height)
            if height % 2 and angle == 0:
                first_pole[exponent] += Fraction(value, 12 + height)
            if height % 2 and angle == -1:
                second_pole[exponent] += Fraction(value, 12 + height)
        rows.append((coherent, first_pole, second_pole))
    return rows


def nonzero_terms(polynomial):
    return {exponent: value for exponent, value in polynomial.items() if value}


def eliminate_poles(row):
    coherent, first_pole, second_pole = row
    result = defaultdict(Fraction, coherent)
    for exponent, value in first_pole.items():
        result[exponent] += Fraction(2, 31) * value
    for (r, v, h), value in second_pole.items():
        result[r, v + 1, h] += Fraction(16, 31) * value
    return nonzero_terms(result)


def specialize(polynomial, index, value):
    result = defaultdict(Fraction)
    for exponent, coefficient in polynomial.items():
        remaining = exponent[:index] + (0,) + exponent[index + 1:]
        result[remaining] += coefficient * value ** exponent[index]
    return nonzero_terms(result)


def eliminate_angular_ratio(polynomial):
    # Substitute h=(427/6)r-4454/67.
    slope, intercept = Fraction(427, 6), Fraction(-4454, 67)
    result = defaultdict(Fraction)
    for (r, v, h), coefficient in polynomial.items():
        for index in range(h + 1):
            result[r + index, v, 0] += (
                coefficient * comb(h, index) * slope ** index * intercept ** (h - index)
            )
    return nonzero_terms(result)


def assert_proportional(actual, expected):
    first = next(iter(expected))
    scale = actual[first] / expected[first]
    assert scale and actual == {exponent: scale * value for exponent, value in expected.items()}


def check_radial_and_angular_parameters():
    rows = reconstruction_rows()
    one, r, v, h = (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)
    assert tuple(map(nonzero_terms, rows[0])) == (
        {one: Fraction(-2, 899)}, {one: Fraction(1, 29)}, {},
    )
    assert tuple(map(nonzero_terms, rows[1])) == (
        {v: Fraction(-4, 2115)}, {v: Fraction(2, 45)}, {one: Fraction(-4, 2115)},
    )
    residuals = [eliminate_poles(row) for row in rows[2:]]
    assert_proportional(specialize(residuals[0], 2, 0), {one: 26724, r: -28609})
    assert_proportional(specialize(residuals[1], 2, 0), {v: 33636, (1, 1, 0): -44737})
    assert_proportional(specialize(residuals[2], 2, 0), {
        (0, 2, 0): 4535708, (1, 2, 0): -7209279,
        (2, 0, 0): -1383599, r: 351436, one: 652208,
    })
    assert 33636 * 28609 - 44737 * 26724 != 0
    assert (
        28609**2 * -351436 + 28609 * 1383599 * 26724 - 26921067152 * 28609 == 0
    )
    assert 28609**2 * -652208 + 26921067152 * 26724 == 185622733382800

    assert_proportional(residuals[0], {h: 1, r: Fraction(-427, 6), one: Fraction(4454, 67)})
    fourth, fifth, sixth, seventh = map(eliminate_angular_ratio, residuals[1:])
    assert_proportional(fourth, {v: -3, (1, 1, 0): 7})
    assert_proportional(specialize(fifth, 0, Fraction(3, 7)), {
        (0, 2, 0): -624512896, one: -649463671,
    })
    # In this branch v is nonzero, so cancel it from the sixth row.
    assert_proportional(specialize(sixth, 0, Fraction(3, 7)), {
        (0, 3, 0): 199006080, v: 222720503,
    })
    assert -624512896 * 222720503 + 649463671 * 199006080 != 0

    q5 = {one: -85804, r: -9164116, (2, 0, 0): 5555305}
    q7 = {
        one: -11534389824, r: -1292362283932,
        (2, 0, 0): -1910019598132, (3, 0, 0): 1417868646401,
    }
    assert_proportional(specialize(fifth, 1, 0), q5)
    assert_proportional(specialize(seventh, 1, 0), q7)
    assert [value % 3 for value in q5.values()] == [2, 2, 1]
    assert [value % 3 for value in q7.values()] == [0, 2, 2, 2]
    print("Conductor 18: exact scalar, radial-pair, and three-parameter row identities")


if __name__ == "__main__":
    check_wider_strip()
    check_upper_half_pole()
    check_radial_and_angular_parameters()

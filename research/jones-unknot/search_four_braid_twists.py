#!/usr/bin/env python3
"""Search the seven-block B4 family using exact Temperley–Lieb traces."""

import argparse
from collections import Counter
from functools import cache
import json
from math import isqrt


STRANDS = 4
GENERATORS = (0, 1, 2)
IDENTITY = tuple((i + STRANDS) % (2 * STRANDS) for i in range(2 * STRANDS))
IDENTITY_JET = (1, 0, 0)
ZERO_JET = (0, 0, 0)


def append_smoothing(pairing, generator):
    """Stack E_i below a pairing, returning the new pairing and loop count."""
    left = STRANDS + generator
    right = left + 1
    left_partner, right_partner = pairing[left], pairing[right]
    if left_partner == right:
        return pairing, 1
    result = list(pairing)
    result[left], result[right] = right, left
    result[left_partner], result[right_partner] = right_partner, left_partner
    return tuple(result), 0


def closure_loops(pairing):
    remaining = set(range(2 * STRANDS))
    loops = 0
    while remaining:
        todo = [remaining.pop()]
        loops += 1
        while todo:
            point = todo.pop()
            for neighbor in (pairing[point], IDENTITY[point]):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    todo.append(neighbor)
    return loops


def diagram_basis():
    diagrams = [IDENTITY]
    for diagram in diagrams:
        for generator in GENERATORS:
            next_diagram, _ = append_smoothing(diagram, generator)
            if next_diagram not in diagrams:
                diagrams.append(next_diagram)
    return diagrams


BASIS = diagram_basis()
assert len(BASIS) == 14
PRODUCTS = {
    (index, generator): (BASIS.index(new_diagram), extra_loop)
    for index, diagram in enumerate(BASIS)
    for generator in GENERATORS
    for new_diagram, extra_loop in (append_smoothing(diagram, generator),)
}
CLOSURE_POWERS = [closure_loops(diagram) - 1 for diagram in BASIS]


def add_jets(left, right):
    return tuple(a + b for a, b in zip(left, right))


def multiply_jets(left, right):
    a, b, c = left
    d, e, f = right
    return a * d, a * e + b * d, a * f + b * e + c * d


def monomial_jet(coefficient, exponent):
    return coefficient, coefficient * exponent, coefficient * (exponent * (exponent - 1) // 2)


@cache
def delta_jet(power):
    # delta = -2 - 4(A - 1)^2 mod (A - 1)^3.
    constant = (-2) ** power
    return constant, 0, 2 * power * constant


@cache
def twist_jets(exponent):
    straight = monomial_jet(1 if exponent % 2 == 0 else -1, -2 * exponent)
    difference = add_jets(monomial_jet(1, -6 * exponent), tuple(-c for c in straight))
    # The E coefficient is difference / delta. Division here is exact.
    constant, linear, quadratic = difference
    assert constant % 2 == linear % 2 == quadratic % 2 == 0
    smoothing = (-constant // 2, -linear // 2, constant - quadratic // 2)
    return straight, (smoothing, difference)


def append_twist_jet(state, generator, exponent):
    straight, smoothings = twist_jets(exponent)
    result = [multiply_jets(coefficient, straight) for coefficient in state]
    for index, coefficient in enumerate(state):
        if coefficient == ZERO_JET:
            continue
        target, loops = PRODUCTS[index, generator]
        contribution = multiply_jets(coefficient, smoothings[loops])
        result[target] = add_jets(result[target], contribution)
    return result


def quadratic_screen(state, generator):
    linear = constant = 0
    for index, coefficient in enumerate(state):
        constant += multiply_jets(coefficient, delta_jet(CLOSURE_POWERS[index]))[2]
        target, loops = PRODUCTS[index, generator]
        linear += coefficient[1] * delta_jet(CLOSURE_POWERS[target] + loops)[0]
    linear *= 4
    # A Jones polynomial of 1 requires 24 k^2 - linear k - constant = 0.
    discriminant = linear * linear + 96 * constant
    if discriminant < 0:
        return ()
    root = isqrt(discriminant)
    if root * root != discriminant:
        return ()
    return tuple(sorted({value // 48 for value in (linear + root, linear - root) if value % 48 == 0}))


def add_shifted(result, polynomial, coefficient, shift):
    for exponent, source_coefficient in polynomial.items():
        target = exponent + shift
        result[target] += coefficient * source_coefficient


def remove_zeros(polynomial):
    return {exponent: coefficient for exponent, coefficient in polynomial.items() if coefficient}


@cache
def delta_laurent(power):
    if power == 0:
        return {0: 1}
    result = Counter()
    for shift in (-2, 2):
        add_shifted(result, delta_laurent(power - 1), -1, shift)
    return result


def normalized_bracket(blocks):
    state = [{0: 1}] + [{} for _ in BASIS[1:]]
    for generator, exponent in blocks:
        sign = 1 if exponent > 0 else -1
        for _ in range(abs(exponent)):
            updated = [Counter() for _ in BASIS]
            for index, coefficient in enumerate(state):
                add_shifted(updated[index], coefficient, -1, -2 * sign)
                target, loops = PRODUCTS[index, generator]
                for shift, scalar in delta_laurent(loops).items():
                    add_shifted(updated[target], coefficient, -scalar, shift - 4 * sign)
            state = [remove_zeros(polynomial) for polynomial in updated]
    result = Counter()
    for coefficient, power in zip(state, CLOSURE_POWERS):
        for shift, scalar in delta_laurent(power).items():
            add_shifted(result, coefficient, scalar, shift)
    return remove_zeros(result)


def second_coefficient(bracket):
    return sum(c * e * (e - 1) // 2 for e, c in bracket.items())


def bigelow_power_coefficients():
    # Bigelow's Corollary 6.2, with [x,y] = x y x^-1 y^-1.
    a, b, c = GENERATORS
    psi = [
        (a, -3), (b, -2), (c, -1), (b, 1), (c, 4), (b, -1), (c, 1),
        (b, 1), (a, 1), (b, 1), (c, 2), (b, 1), (a, -1), (b, -1), (c, -2),
    ]

    def inverse(word):
        return [(generator, -exponent) for generator, exponent in reversed(word)]

    x = [(b, 1), (a, 1)] * 3
    y = inverse(psi) + [(b, 1)] + psi
    beta = x + y + inverse(x) + inverse(y)
    tail = [(a, 1), (b, 1), (c, 1)]
    samples = {
        k: second_coefficient(normalized_bracket(prefix + tail))
        for k, prefix in ((-1, inverse(beta)), (0, []), (1, beta))
    }
    assert samples == {-1: 384, 0: 0, 1: -384}
    return samples


def cyclically_cancel(word):
    for start, first in enumerate(word):
        for offset in range(1, len(word)):
            end = (start + offset) % len(word)
            second = word[end]
            if second == -first:
                return [letter for pos, letter in enumerate(word) if pos not in (start, end)]
            if abs(abs(first) - abs(second)) <= 1:
                break
    return word


def shortened_length(blocks):
    """Shorten by cyclic commutations/cancellations and outer destabilizations."""
    word = [
        (generator + 1) if exponent > 0 else -(generator + 1)
        for generator, exponent in blocks
        for _ in range(abs(exponent))
    ]
    left, right = 1, STRANDS - 1
    while word:
        reduced = cyclically_cancel(word)
        if len(reduced) < len(word):
            word = reduced
        elif sum(abs(letter) == left for letter in word) == 1:
            word = [letter for letter in word if abs(letter) != left]
            left += 1
        elif sum(abs(letter) == right for letter in word) == 1:
            word = [letter for letter in word if abs(letter) != right]
            right -= 1
        else:
            break
    return len(word)


def search(odd_max, even_max):
    odd = range(-odd_max, odd_max + 1, 2)
    even = range(-even_max, even_max + 1, 2)
    ranges = ((0, odd), (1, odd), (2, odd), (1, even), (0, even), (2, even))
    final_generator = 1
    prefixes = root_count = unit_count = 0
    max_root_abs = max_unit_length = long_unit_count = 0
    max_shortened_length = 0
    unresolved_units = []

    def visit(blocks, state):
        nonlocal prefixes, root_count, unit_count, max_root_abs, max_unit_length
        nonlocal long_unit_count, max_shortened_length
        if len(blocks) < len(ranges):
            generator, values = ranges[len(blocks)]
            for exponent in values:
                visit(blocks + [(generator, exponent)], append_twist_jet(state, generator, exponent))
            return
        prefixes += 1
        for k in quadratic_screen(state, final_generator):
            root_count += 1
            max_root_abs = max(max_root_abs, abs(k))
            full_word = blocks + [(final_generator, 2 * k)]
            bracket = normalized_bracket(full_word)
            assert second_coefficient(bracket) == 0
            if bracket == {0: 1}:
                unit_count += 1
                length = sum(abs(exponent) for _, exponent in full_word)
                max_unit_length = max(max_unit_length, length)
                if length > 24:
                    long_unit_count += 1
                    reduced_length = shortened_length(full_word)
                    max_shortened_length = max(max_shortened_length, reduced_length)
                    if reduced_length > 24:
                        unresolved_units.append([exponent for _, exponent in full_word])

    visit([], [IDENTITY_JET] + [ZERO_JET] * (len(BASIS) - 1))
    return {
        "odd_max": odd_max,
        "even_max": even_max,
        "prefixes": prefixes,
        "integer_roots": root_count,
        "max_abs_k_before_full_check": max_root_abs,
        "unit_jones_words": unit_count,
        "max_unit_word_length": max_unit_length,
        "unit_jones_words_above_24_before_shortening": long_unit_count,
        "max_shortened_length_of_these_words": max_shortened_length,
        "unit_word_exponents_above_24_after_shortening": unresolved_units,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--odd-max", type=int, default=7)
    parser.add_argument("--even-max", type=int, default=6)
    args = parser.parse_args()
    if args.odd_max < 1 or args.odd_max % 2 == 0 or args.even_max < 0 or args.even_max % 2:
        parser.error("the first bound must be positive odd and the second nonnegative even")
    assert normalized_bracket([(0, 1), (1, -1), (2, 1)]) == {0: 1}
    assert normalized_bracket([(0, 3), (1, 1), (2, 1)]) == {-16: -1, -12: 1, -4: 1}
    print(json.dumps(
        {
            "bigelow_power_second_coefficients": bigelow_power_coefficients(),
            "seven_block_search": search(args.odd_max, args.even_max),
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()

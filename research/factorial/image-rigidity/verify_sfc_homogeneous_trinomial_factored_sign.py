#!/usr/bin/env python3
"""Verify the factored sign and leading part of the initial SFC remainder."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from math import comb, gcd, prod

from flint import fmpz_mpoly_ctx

RING = fmpz_mpoly_ctx.get(("L", "H"))
LEFT, RIGHT = VARIABLES = RING.gens()


def count_terms(order, y_degree, middle, span):
    # i and j count the X and Y factors using the uneliminated columns.
    x_degree = order - y_degree
    for i in range(x_degree + 1):
        for j in range(y_degree + 1):
            right_offset = i * span + j * (span - middle)
            multiplier = (
                comb(order, y_degree) * comb(x_degree, i) * comb(y_degree, j)
                * (-1) ** (order - i - j)
            )
            yield multiplier, right_offset, i, j


@lru_cache(maxsize=4096)
def rising_factors(axis, step, length):
    factors = Counter()
    scalar = 1
    for shift in range(1, length + 1):
        divisor = gcd(step, shift)
        factors[axis, step // divisor, shift // divisor] += 1
        scalar *= divisor
    return factors, scalar


def expand_factors(factors):
    polynomials = [
        (step * VARIABLES[axis] + shift) ** multiplicity
        for (axis, step, shift), multiplicity in sorted(factors.items())
    ]
    if not polynomials:
        return RING.constant(1)
    while len(polynomials) > 1:
        last_pair = len(polynomials) // 2 * 2
        polynomials = [
            polynomials[i] * polynomials[i + 1] for i in range(0, last_pair, 2)
        ] + polynomials[last_pair:]
    return polynomials[0]


def coefficient(order, y_degree, middle, span):
    y_gap = span - middle
    terms = []
    for multiplier, right_offset, i, j in count_terms(order, y_degree, middle, span):
        factors = Counter()
        for axis, step, length, multiplicity in (
            (0, order, order * span - right_offset, 1),
            (1, order, right_offset, 1),
            (0, 1, span, i + j),
            (0, 1, middle, y_degree - j),
            (1, 1, span, order - y_degree - i),
            (1, 1, y_gap, y_degree - j),
        ):
            if multiplicity:
                part, scalar = rising_factors(axis, step, length)
                factors.update({key: count * multiplicity for key, count in part.items()})
                multiplier *= scalar ** multiplicity
        terms.append((multiplier, factors))

    common = terms[0][1].copy()
    for _, factors in terms[1:]:
        common &= factors
    expected = Counter({
        (0, 1, r): 1 + y_degree * (r <= middle) for r in range(1, span + 1)
    })
    for shift in range(1, y_degree * middle + 1):
        if shift % order:
            divisor = gcd(order, shift)
            expected[0, order // divisor, shift // divisor] += 1
    right_length = max(y_gap, (order - y_degree) * span // order)
    expected.update((1, 1, r) for r in range(1, right_length + 1))
    if common != expected:
        raise ArithmeticError(f"Common factors disagree at {(order, y_degree, middle, span)}.")
    quotient = sum(multiplier * expand_factors(factors - common) for multiplier, factors in terms)
    return common, quotient


def factored_negative_remainder(middle, span):
    (a_factors, a), (b_factors, b), (c_factors, c) = (
        coefficient(2, k, middle, span) for k in range(3)
    )
    (d_factors, d), (f_factors, f), (g_factors, g) = (
        coefficient(3, k, middle, span) for k in (0, 2, 3)
    )
    terms = (
        (c_factors + c_factors + d_factors, c, c, d),
        (a_factors + b_factors + g_factors, a, b, g),
        (a_factors + c_factors + f_factors, -a, c, f),
    )
    common = terms[0][0] & terms[1][0] & terms[2][0]
    quotient = sum(p * q * r * expand_factors(factors - common) for factors, p, q, r in terms)
    return common, quotient


def verify_profile(profile):
    middle, span = profile
    common, quotient = factored_negative_remainder(middle, span)
    terms = quotient.to_dict()
    if terms.get((0, 0), 0) <= 0 or min(terms.values()) < 0:
        full = (expand_factors(common) * quotient).to_dict()
        witness = next(((key, value) for key, value in full.items() if value < 0), None)
        if witness is not None or full.get((0, 0), 0) <= 0:
            raise ArithmeticError(
                f"Sign test for -R failed at {profile}. "
                f"Constant: {full.get((0, 0), 0)}. Negative coefficient: {witness}."
            )

    x_degree = sum(n for (axis, _, _), n in common.items() if axis == 0)
    y_degree = sum(n for (axis, _, _), n in common.items() if axis == 1)
    scalar = prod(step ** n for (_, step, _), n in common.items())
    degree = max(map(sum, terms))
    top = RING.from_dict({key: value for key, value in terms.items() if sum(key) == degree})
    gap = span - middle
    scale = (2 ** (4 * span) * 3 ** (3 * span) * middle**2 * span**3 * gap**4) // 72
    expected = (
        scale
        * LEFT ** (7 * span + 4 * middle - 5)
        * RIGHT ** (7 * span - 4 * middle - 5)
        * (LEFT + RIGHT)**4
        * ((span - 2) * LEFT + (span + 2) * RIGHT)
    )
    if scalar * LEFT**x_degree * RIGHT**y_degree * top != expected:
        raise ArithmeticError(f"Leading part disagrees at {profile}.")


def verify_unfactored_comparison():
    from verify_sfc_homogeneous_trinomial_short_spans import initial_remainder

    for span in range(2, 7):
        for middle in range(1, span):
            common, quotient = factored_negative_remainder(middle, span)
            expanded = expand_factors(common) * quotient
            legacy = initial_remainder(middle, span)
            if expanded.to_dict() != {key: -int(value) for key, value in legacy.items()}:
                raise ArithmeticError(f"Unfactored comparison failed at {(middle, span)}.")


def verify_symbolic_leading_part():
    from sympy import expand, symbols

    middle, span, x, y = symbols("a d x y")

    def power_sum(size, power):
        if power == 1:
            return size * (size + 1) / 2
        if power == 2:
            return size * (size + 1) * (2 * size + 1) / 6
        return (size * (size + 1) / 2)**2

    def normalized_coefficients(order, y_degree):
        gap = span - middle
        result = [0, 0, 0, 0]
        for multiplier, right_offset, i, j in count_terms(order, y_degree, middle, span):
            log_coefficients = []
            for power in (1, 2, 3):
                left = (
                    power_sum(order * span - right_offset, power) / order**power
                    + (i + j) * power_sum(span, power)
                    + (y_degree - j) * power_sum(middle, power)
                )
                right = (
                    power_sum(right_offset, power) / order**power
                    + (order - y_degree - i) * power_sum(span, power)
                    + (y_degree - j) * power_sum(gap, power)
                )
                log_coefficients.append(
                    (-1)**(power - 1) * (left * x**power + right * y**power) / power
                )
            u, v, w = log_coefficients
            for index, value in enumerate((1, u, v + u**2 / 2, w + u*v + u**3 / 6)):
                result[index] += multiplier * value
        return list(map(expand, result))

    a, b, c = (normalized_coefficients(2, k) for k in range(3))
    d, f, g = (normalized_coefficients(3, k) for k in (0, 2, 3))
    if any(poly[0] for poly in (a, b, c)) or any(poly[:2] != [0, 0] for poly in (d, f, g)):
        raise ArithmeticError("Unexpected terms of order below four in -R.")
    fourth = c[1]**2 * d[2] + a[1]*b[1]*g[2] - a[1]*c[1]*f[2]
    if expand(fourth) != 0:
        raise ArithmeticError("The order-four term does not vanish.")
    quadratic_correction = (
        2*c[1]*c[2]*d[2] + (a[1]*b[2] + a[2]*b[1])*g[2]
        - (a[1]*c[2] + a[2]*c[1])*f[2]
    )
    cubic_correction = c[1]**2*d[3] + a[1]*b[1]*g[3] - a[1]*c[1]*f[3]
    factor = middle**2 * span**3 * (span - middle)**4 * (x + y)**4 / 72
    if expand(quadratic_correction + 3 * factor * (x - y)) != 0:
        raise ArithmeticError("The quadratic correction disagrees.")
    if expand(cubic_correction - factor * (span*(x + y) + 5*(x - y))) != 0:
        raise ArithmeticError("The cubic correction disagrees.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-span", type=int, default=2)
    parser.add_argument("--max-span", type=int, default=64)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    if not 2 <= args.min_span <= args.max_span or args.jobs < 1:
        parser.error("Use 2 <= min-span <= max-span and jobs >= 1.")
    verify_symbolic_leading_part()
    verify_unfactored_comparison()
    print("Symbolic leading identity and unfactored comparisons passed.", flush=True)
    checked = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for span in range(args.min_span, args.max_span + 1):
            profiles = [(middle, span) for middle in range(1, span)]
            for _ in pool.map(verify_profile, profiles):
                checked += 1
            print(f"Span {span}: all {span - 1} oriented profiles passed.", flush=True)
    print(f"Verified {checked} profiles at all nonnegative real heights.")


if __name__ == "__main__":
    main()

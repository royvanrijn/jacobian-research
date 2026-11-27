#!/usr/bin/env python3
"""Verify odd quadratic Gaussian windows at primes dividing a tested order."""

import argparse
from fractions import Fraction


def prime_divisors(number, minimum):
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            if divisor >= minimum:
                yield divisor
            while number % divisor == 0:
                number //= divisor
        divisor += 1
    if number >= minimum:
        yield number


def remove_prime(number, prime):
    exponent = 0
    while number % prime == 0:
        number //= prime
        exponent += 1
    return exponent, number


def factorial_tables(limit, prime):
    valuations = [0]
    units = [1]
    for number in range(1, limit + 1):
        exponent, unit = remove_prime(number, prime)
        valuations.append(valuations[-1] + exponent)
        units.append(units[-1] * unit % prime)
    return valuations, units


def positive_root_valuations(coefficients):
    vertices = []
    for degree, valuation in enumerate(coefficients):
        while len(vertices) >= 2:
            left, previous = vertices[-2:]
            if (previous[1] - left[1]) * (degree - previous[0]) < (
                valuation - previous[1]
            ) * (previous[0] - left[0]):
                break
            vertices.pop()
        vertices.append((degree, valuation))
    # These are polynomials in x**2.
    return {
        Fraction(left[1] - right[1], 2 * (right[0] - left[0]))
        for left, right in zip(vertices, vertices[1:])
        if left[1] > right[1]
    }


def initial_polynomial(valuations, units, slope):
    weights = [value + 2 * index * slope for index, value in enumerate(valuations)]
    minimum = min(weights)
    support = [index for index, value in enumerate(weights) if value == minimum]
    return [
        units[index] if weights[index] == minimum else 0
        for index in range(support[0], support[-1] + 1)
    ]


def gcd_degree(first, second, prime):
    first = first[:]
    second = second[:]
    while second:
        inverse = pow(second[-1], -1, prime)
        while len(first) >= len(second):
            offset = len(first) - len(second)
            leading = first[-1] * inverse
            for index, value in enumerate(second, offset):
                first[index] = (first[index] - leading * value) % prime
            while first and first[-1] == 0:
                first.pop()
        first, second = second, first
    return len(first) - 1


def quadratic_bound(order, prime, double_factorials):
    exponent, _ = remove_prime(order, prime)
    return min(
        Fraction(
            exponent * removed + double_factorials[order - removed]
            - double_factorials[order],
            2 * removed,
        )
        for removed in range(1, order + 1)
    )


def pure_hermite(order, prime, valuations, units):
    parity = order % 2
    half = (order - parity) // 2
    return (
        [
            valuations[half] - valuations[half - i] - valuations[2 * i + parity]
            for i in range(half + 1)
        ],
        [
            pow(2 * order, i, prime)
            * units[half]
            * pow(units[half - i] * units[2 * i + parity], -1, prime)
            % prime
            for i in range(half + 1)
        ],
    )


def disjoint_initial_forms(first, second, prime):
    common_slopes = (
        positive_root_valuations(first[0]) & positive_root_valuations(second[0])
    )
    for slope in common_slopes:
        initial_first = initial_polynomial(*first, slope)
        initial_second = initial_polynomial(*second, slope)
        if gcd_degree(initial_first, initial_second, prime) > 0:
            return False
    return True


def verify_middle(middle, prime):
    valuations, units = factorial_tables(2 * middle + 2, prime)
    double_factorials = [valuations[2 * i] - valuations[i] for i in range(middle + 2)]
    odd_orders = (middle - 1, middle + 1)
    lower = quadratic_bound(middle, prime, double_factorials)
    higher_term_bound = max(
        Fraction(
            double_factorials[(order - 1) // 2]
            - double_factorials[(order - degree) // 2],
            degree - 3,
        )
        for order in odd_orders
        for degree in range(5, order + 1, 2)
    )
    if lower <= max(0, higher_term_bound):
        return False
    first, second = [pure_hermite(m, prime, valuations, units) for m in odd_orders]
    return disjoint_initial_forms(first, second, prime)


def verify_linear_chart(lower, even_order, odd_order, prime, valuations, units):
    even = pure_hermite(even_order, prime, valuations, units)
    if lower <= max(positive_root_valuations(even[0]), default=0):
        return False
    odd = pure_hermite(odd_order, prime, valuations, units)
    return disjoint_initial_forms(even, odd, prime)


def verify_first(start, prime):
    valuations, units = factorial_tables(2 * start + 4, prime)
    double_factorials = [valuations[2 * i] - valuations[i] for i in range(start + 1)]
    lower = quadratic_bound(start, prime, double_factorials)
    return verify_linear_chart(lower, start + 1, start + 2, prime, valuations, units)


def linear_moment_reduction(order, prime, valuations, units):
    polynomial = [0] * (order + 1)
    for r in range(order // 2 + 1):
        if valuations[order] - valuations[r] - valuations[order - 2 * r] > 0:
            continue
        polynomial[order - 2 * r] = (
            units[order] * pow(units[r] * units[order - 2 * r], -1, prime)
            * pow(2 * order, -r, prime)
        ) % prime
    return polynomial


def verify_last(last, prime):
    valuations, units = factorial_tables(2 * last, prime)
    double_factorials = [valuations[2 * i] - valuations[i] for i in range(last + 1)]
    lower = quadratic_bound(last, prime, double_factorials)
    if lower <= 0:
        return False
    binary = [
        linear_moment_reduction(m, prime, valuations, units) for m in (last - 2, last - 1)
    ]
    powers = [next(i for i, c in enumerate(poly) if c) for poly in binary]
    first, second = [poly[power:] for poly, power in zip(binary, powers, strict=True)]
    if gcd_degree(first, second, prime) != 0:
        return False
    if min(powers) == 0:
        return True
    return verify_linear_chart(lower, last - 1, last - 2, prime, valuations, units)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first", type=int, help="first odd start")
    parser.add_argument("last", type=int, help="last odd start")
    parser.add_argument("--expect-unresolved", type=int, nargs="*", default=[])
    parser.add_argument("--middle-only", action="store_true")
    args = parser.parse_args()
    if args.first < 9 or args.last < args.first or any(
        order % 2 == 0 for order in (args.first, args.last)
    ):
        parser.error("Supply an increasing range of odd starts beginning at nine or later.")

    verified = []
    unresolved = []
    for start in range(args.first, args.last + 1, 2):
        tests = [(start + 1, 5, verify_middle)]
        if not args.middle_only:
            tests.extend([(start, 3, verify_first), (start + 2, 3, verify_last)])
        certified = next(
            (
                (order, prime)
                for order, minimum, verify in tests
                for prime in prime_divisors(order, minimum)
                if verify(order, prime)
            ),
            None,
        )
        if certified is None:
            unresolved.append(start)
        else:
            verified.append((start, *certified))
    print(f"Certified {len(verified)} odd starts from {args.first} through {args.last}.")
    print(f"Unresolved starts under this test: {unresolved}")
    if unresolved != args.expect_unresolved:
        raise ArithmeticError("The unresolved starts differ from the specified list.")


if __name__ == "__main__":
    main()

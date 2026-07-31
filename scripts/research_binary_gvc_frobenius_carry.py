#!/usr/bin/env python3
"""Exact search for the binary jet--carry score and Frobenius gap.

For channels on one ordinary-homogeneous line of degree ``r``, enumerate
selection vectors in moments of order ``p*m``.  The prime is chosen larger
than ``r*m``.  The script verifies the exact valuation identity

    v_p(mult_op * mult_pol * rho_x! * rho_y!)
      = r*m + e_op + e_pol - c_rad

and checks that equality with the Frobenius floor ``r*m`` occurs exactly
when both selection vectors are divisible by ``p``.  It also classifies
the first correction shell, where the excess is one.

This is an exact bounded regression for Lemma 7.4 ter in
``BINARY_GVC_UNIFORM_FACE_TERMINATION.md``.  The written proof is
unbounded.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import factorial


@dataclass(frozen=True)
class SelectionType:
    entries: tuple[int, ...]
    carry: int
    divisible: bool
    full_support: bool


def compositions(total: int, length: int):
    """Yield ordered weak compositions of ``total``."""

    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def is_prime(number: int) -> bool:
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def next_prime(number: int) -> int:
    candidate = number + 1
    while not is_prime(candidate):
        candidate += 1
    return candidate


def valuation(number: int, prime: int) -> int:
    answer = 0
    while number and number % prime == 0:
        answer += 1
        number //= prime
    return answer


def unit_mod_prime(number: int, prime: int, order: int) -> int:
    return (number // prime**order) % prime


def multinomial(total: int, entries: tuple[int, ...]) -> int:
    answer = factorial(total)
    for entry in entries:
        answer //= factorial(entry)
    return answer


def support_subsets(radial_degree: int):
    levels = tuple(range(radial_degree + 1))
    for size in range(1, len(levels) + 1):
        yield from combinations(levels, size)


def verify_primitive_bridge_support(limit: int) -> None:
    """Construct the primitive one-carry support reductions.

    Every triple of distinct levels supports a one-sided carry bridge of
    total ``p``.  Every pair of distinct levels realizes every nonzero
    radial residue, so two arbitrary pairs support a two-sided bridge.
    """

    checked_triples = 0
    checked_pair_pairs = 0
    for radial_degree in range(2, limit + 1):
        prime = next_prime(radial_degree)
        levels = tuple(range(radial_degree + 1))

        for first, second, third in combinations(levels, 3):
            counts = (
                (second - third) % prime,
                (third - first) % prime,
                (first - second) % prime,
            )
            assert all(counts)
            if sum(counts) == 2 * prime:
                counts = tuple(prime - count for count in counts)
            assert sum(counts) == prime
            assert (
                first * counts[0]
                + second * counts[1]
                + third * counts[2]
            ) % prime == 0
            checked_triples += 1

        pairs = tuple(combinations(levels, 2))
        for operator_pair in pairs:
            for polynomial_pair in pairs:
                target_residue = 1
                operator_count = (
                    target_residue
                    * pow(
                        operator_pair[0] - operator_pair[1],
                        -1,
                        prime,
                    )
                ) % prime
                polynomial_count = (
                    target_residue
                    * pow(
                        polynomial_pair[0] - polynomial_pair[1],
                        -1,
                        prime,
                    )
                ) % prime
                assert operator_count
                assert polynomial_count
                operator_residue = (
                    operator_count * operator_pair[0]
                    + (prime - operator_count) * operator_pair[1]
                ) % prime
                polynomial_residue = (
                    polynomial_count * polynomial_pair[0]
                    + (prime - polynomial_count) * polynomial_pair[1]
                ) % prime
                assert operator_residue == target_residue
                assert polynomial_residue == target_residue
                checked_pair_pairs += 1

    print(
        "PASS primitive one-carry support reduction: "
        f"{checked_triples} triples and {checked_pair_pairs} pair-pairs "
        f"through radial degree {limit}"
    )


def verify_nonhomogeneous_jet_score() -> None:
    """Check the general jet--carry formula on a mixed-degree support."""

    channels = ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
    baseline_degree = 1
    largest_degree = 2
    tested = 0

    for base_order in (1, 2):
        prime = next_prime(base_order * largest_degree)
        moment_order = prime * base_order
        grouped: dict[
            tuple[int, int],
            dict[tuple[int, bool], SelectionType],
        ] = defaultdict(dict)

        for entries in compositions(moment_order, len(channels)):
            remainders = tuple(entry % prime for entry in entries)
            carry = sum(remainders) // prime
            divisible = carry == 0
            radial = (
                sum(entry * channel[0] for entry, channel in zip(entries, channels)),
                sum(entry * channel[1] for entry, channel in zip(entries, channels)),
            )
            grouped[radial].setdefault(
                (carry, divisible),
                SelectionType(
                    entries,
                    carry,
                    divisible,
                    all(entries),
                ),
            )

        for radial, type_map in grouped.items():
            types = tuple(type_map.values())
            jet_excess = (
                radial[0]
                + radial[1]
                - baseline_degree * moment_order
            )
            radial_remainders = (
                radial[0] % prime,
                radial[1] % prime,
            )
            carry_radial = (
                sum(radial_remainders) - jet_excess % prime
            ) // prime
            assert carry_radial in (0, 1)
            radial_value = factorial(radial[0]) * factorial(radial[1])
            radial_valuation = valuation(radial_value, prime)

            for operator in types:
                operator_multinomial = multinomial(
                    moment_order, operator.entries
                )
                operator_valuation = valuation(
                    operator_multinomial, prime
                )
                assert operator_valuation == operator.carry
                for polynomial in types:
                    polynomial_multinomial = multinomial(
                        moment_order, polynomial.entries
                    )
                    polynomial_valuation = valuation(
                        polynomial_multinomial, prime
                    )
                    assert polynomial_valuation == polynomial.carry
                    exact = (
                        operator_valuation
                        + polynomial_valuation
                        + radial_valuation
                    )
                    carry_penalty = (
                        operator.carry
                        + polynomial.carry
                        - carry_radial
                    )
                    predicted = (
                        baseline_degree * base_order
                        + jet_excess // prime
                        + carry_penalty
                    )
                    assert exact == predicted
                    assert carry_penalty >= 0
                    assert (carry_penalty == 0) == (
                        operator.divisible
                        and polynomial.divisible
                    )
                    tested += 1

    print(
        "PASS nonhomogeneous jet--carry score: "
        f"{tested} exact mixed-degree return types"
    )


def selection_table(
    support: tuple[int, ...],
    moment_order: int,
    prime: int,
) -> dict[int, tuple[SelectionType, ...]]:
    """Group selection types by their first radial coordinate."""

    grouped: dict[
        int, dict[tuple[int, bool, bool], SelectionType]
    ] = defaultdict(dict)
    for entries in compositions(moment_order, len(support)):
        remainders = tuple(entry % prime for entry in entries)
        remainder_sum = sum(remainders)
        assert remainder_sum % prime == 0
        carry = remainder_sum // prime
        divisible = remainder_sum == 0
        full_support = all(entries)
        radial_x = sum(entry * level for entry, level in zip(entries, support))
        key = (carry, divisible, full_support)
        grouped[radial_x].setdefault(
            key,
            SelectionType(entries, carry, divisible, full_support),
        )
    return {
        radial_x: tuple(types.values())
        for radial_x, types in grouped.items()
    }


def verify_window(radial_limit: int, order_limit: int) -> None:
    tested_pairs = 0
    tested_types = 0
    gap_one = {
        "one_side_carry": 0,
        "two_side_carry_radial": 0,
    }
    full_support_gap_one = 0
    first_full_support_examples: dict[
        str,
        tuple[
            int,
            int,
            int,
            tuple[int, ...],
            tuple[int, ...],
            int,
            tuple[int, ...],
            tuple[int, ...],
        ],
    ] = {}

    for radial_degree in range(1, radial_limit + 1):
        supports = tuple(support_subsets(radial_degree))
        for base_order in range(1, order_limit + 1):
            prime = next_prime(radial_degree * base_order)
            moment_order = prime * base_order
            tables = {
                support: selection_table(support, moment_order, prime)
                for support in supports
            }
            frobenius_floor = radial_degree * base_order

            for operator_support in supports:
                operator_table = tables[operator_support]
                for polynomial_support in supports:
                    polynomial_table = tables[polynomial_support]
                    common_radials = operator_table.keys() & polynomial_table.keys()
                    if not common_radials:
                        continue
                    tested_pairs += 1

                    for radial_x in common_radials:
                        radial_y = (
                            radial_degree * moment_order - radial_x
                        )
                        radial_carry = int(radial_x % prime != 0)
                        radial_valuation = valuation(
                            factorial(radial_x), prime
                        ) + valuation(factorial(radial_y), prime)
                        assert (
                            radial_valuation
                            == frobenius_floor - radial_carry
                        )
                        radial_unit = unit_mod_prime(
                            factorial(radial_x) * factorial(radial_y),
                            prime,
                            radial_valuation,
                        )
                        radial_quotients = (
                            radial_x // prime,
                            radial_y // prime,
                        )
                        radial_remainders = (
                            radial_x % prime,
                            radial_y % prime,
                        )
                        predicted_radial_unit = (
                            (-1) ** sum(radial_quotients)
                            * factorial(radial_quotients[0])
                            * factorial(radial_quotients[1])
                            * factorial(radial_remainders[0])
                            * factorial(radial_remainders[1])
                        ) % prime
                        assert radial_unit == predicted_radial_unit

                        for operator in operator_table[radial_x]:
                            operator_multinomial = multinomial(
                                moment_order, operator.entries
                            )
                            operator_valuation = valuation(
                                operator_multinomial, prime
                            )
                            assert operator_valuation == operator.carry
                            operator_unit = unit_mod_prime(
                                operator_multinomial,
                                prime,
                                operator_valuation,
                            )
                            operator_denominator = 1
                            for entry in operator.entries:
                                operator_denominator *= factorial(
                                    entry // prime
                                )
                                operator_denominator *= factorial(
                                    entry % prime
                                )
                            predicted_operator_unit = (
                                (-1) ** operator.carry
                                * factorial(base_order)
                                * pow(
                                    operator_denominator % prime,
                                    -1,
                                    prime,
                                )
                            ) % prime
                            assert operator_unit == predicted_operator_unit

                            for polynomial in polynomial_table[radial_x]:
                                polynomial_multinomial = multinomial(
                                    moment_order,
                                    polynomial.entries,
                                )
                                polynomial_valuation = valuation(
                                    polynomial_multinomial,
                                    prime,
                                )
                                assert polynomial_valuation == polynomial.carry
                                polynomial_unit = unit_mod_prime(
                                    polynomial_multinomial,
                                    prime,
                                    polynomial_valuation,
                                )
                                polynomial_denominator = 1
                                for entry in polynomial.entries:
                                    polynomial_denominator *= factorial(
                                        entry // prime
                                    )
                                    polynomial_denominator *= factorial(
                                        entry % prime
                                    )
                                predicted_polynomial_unit = (
                                    (-1) ** polynomial.carry
                                    * factorial(base_order)
                                    * pow(
                                        polynomial_denominator % prime,
                                        -1,
                                        prime,
                                    )
                                ) % prime
                                assert (
                                    polynomial_unit
                                    == predicted_polynomial_unit
                                )

                                excess = (
                                    operator.carry
                                    + polynomial.carry
                                    - radial_carry
                                )
                                exact_valuation = (
                                    operator_valuation
                                    + polynomial_valuation
                                    + radial_valuation
                                )
                                assert exact_valuation == frobenius_floor + excess
                                assert excess >= 0
                                assert (excess == 0) == (
                                    operator.divisible
                                    and polynomial.divisible
                                )

                                if excess == 1:
                                    if radial_carry:
                                        assert operator.carry == 1
                                        assert polynomial.carry == 1
                                        category = "two_side_carry_radial"
                                    else:
                                        assert (
                                            operator.carry
                                            + polynomial.carry
                                            == 1
                                        )
                                        category = "one_side_carry"
                                    gap_one[category] += 1
                                    if (
                                        operator.full_support
                                        and polynomial.full_support
                                    ):
                                        full_support_gap_one += 1
                                        first_full_support_examples.setdefault(
                                            category,
                                            (
                                                radial_degree,
                                                base_order,
                                                prime,
                                                operator_support,
                                                polynomial_support,
                                                radial_x,
                                                operator.entries,
                                                polynomial.entries,
                                            ),
                                        )
                                tested_types += 1

    print(
        "PASS binary Frobenius carry gap: "
        f"{tested_types} exact return types across {tested_pairs} "
        "support/order pairs"
    )
    print(f"gap-one shell: {gap_one}")
    print(f"full-support gap-one types: {full_support_gap_one}")
    for category, example in first_full_support_examples.items():
        print(f"first {category} full-support example: {example}")
    print(
        "STATUS: bounded regression; the unbounded result follows from "
        "Legendre--Kummer carry identities"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radial-limit", type=int, default=3)
    parser.add_argument("--order-limit", type=int, default=3)
    parser.add_argument("--bridge-limit", type=int, default=40)
    arguments = parser.parse_args()
    verify_primitive_bridge_support(arguments.bridge_limit)
    verify_nonhomogeneous_jet_score()
    verify_window(arguments.radial_limit, arguments.order_limit)


if __name__ == "__main__":
    main()

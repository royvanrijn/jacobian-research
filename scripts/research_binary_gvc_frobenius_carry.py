#!/usr/bin/env python3
"""Exact search for the binary jet--carry score and affine digit gap.

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
unbounded.  The two-digit checks additionally test the affine extension
at orders ``p*m+r``:

    v_p(multinomial(p*m+r; u))
      = (sum_i (u_i mod p)-r)/p,

with the corresponding factorial-unit formula.  These identities are the
input for the still-open carry-promotion problem; verifying them does not
prove that a partial Hall shell inherits a pure-zero identity.
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


def verify_nonhomogeneous_two_digit_score(residue_limit: int) -> None:
    """Check the affine score on the same mixed-degree support."""

    channels = ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
    baseline_degree = 1
    largest_degree = 2
    tested_selections = 0
    tested_pairs = 0

    for high_digit in (1, 2):
        prime = next_prime(largest_degree * (high_digit + 1))
        for low_digit in range(min(prime - 1, residue_limit) + 1):
            moment_order = prime * high_digit + low_digit
            grouped: dict[
                tuple[int, int],
                dict[tuple[int, int], SelectionType],
            ] = defaultdict(dict)

            for entries in compositions(moment_order, len(channels)):
                remainders = tuple(entry % prime for entry in entries)
                carry = (sum(remainders) - low_digit) // prime
                assert carry >= 0
                value = multinomial(moment_order, entries)
                assert valuation(value, prime) == carry
                denominator = 1
                for entry in entries:
                    denominator *= factorial(entry // prime)
                    denominator *= factorial(entry % prime)
                predicted_unit = (
                    (-1) ** carry
                    * factorial(high_digit)
                    * factorial(low_digit)
                    * pow(denominator % prime, -1, prime)
                ) % prime
                assert unit_mod_prime(value, prime, carry) == predicted_unit
                radial = (
                    sum(
                        entry * channel[0]
                        for entry, channel in zip(entries, channels)
                    ),
                    sum(
                        entry * channel[1]
                        for entry, channel in zip(entries, channels)
                    ),
                )
                grouped[radial].setdefault(
                    (carry, denominator % prime),
                    SelectionType(entries, carry, False, all(entries)),
                )
                tested_selections += 1

            for radial, type_map in grouped.items():
                jet_excess = (
                    radial[0]
                    + radial[1]
                    - baseline_degree * moment_order
                )
                radial_remainders = (
                    radial[0] % prime,
                    radial[1] % prime,
                )
                expected_residue = (
                    baseline_degree * low_digit + jet_excess
                ) % prime
                radial_carry = (
                    sum(radial_remainders) - expected_residue
                ) // prime
                assert radial_carry in (0, 1)
                radial_value = factorial(radial[0]) * factorial(radial[1])
                radial_order = valuation(radial_value, prime)
                predicted_radial_order = (
                    baseline_degree * high_digit
                    + (
                        baseline_degree * low_digit + jet_excess
                    ) // prime
                    - radial_carry
                )
                assert radial_order == predicted_radial_order

                types = tuple(type_map.values())
                for operator in types:
                    for polynomial in types:
                        exact = (
                            operator.carry
                            + polynomial.carry
                            + radial_order
                        )
                        predicted = (
                            baseline_degree * high_digit
                            + (
                                baseline_degree * low_digit
                                + jet_excess
                            ) // prime
                            + operator.carry
                            + polynomial.carry
                            - radial_carry
                        )
                        assert exact == predicted
                        tested_pairs += 1

    print(
        "PASS nonhomogeneous two-digit jet--carry score: "
        f"{tested_selections} selections and {tested_pairs} return types"
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


def two_digit_selection_table(
    support: tuple[int, ...],
    moment_order: int,
    prime: int,
    low_digit: int,
) -> dict[int, tuple[SelectionType, ...]]:
    """Group representatives for selections at order ``p*d+r``."""

    grouped: dict[
        int, dict[tuple[int, bool, int], SelectionType]
    ] = defaultdict(dict)
    for entries in compositions(moment_order, len(support)):
        remainders = tuple(entry % prime for entry in entries)
        remainder_sum = sum(remainders)
        assert remainder_sum % prime == low_digit
        carry = (remainder_sum - low_digit) // prime
        denominator = 1
        for entry in entries:
            denominator *= factorial(entry // prime)
            denominator *= factorial(entry % prime)
        full_support = all(entries)
        radial_x = sum(
            entry * level for entry, level in zip(entries, support)
        )
        key = (carry, full_support, denominator % prime)
        grouped[radial_x].setdefault(
            key,
            SelectionType(
                entries,
                carry,
                False,
                full_support,
            ),
        )
    return {
        radial_x: tuple(types.values())
        for radial_x, types in grouped.items()
    }


def verify_two_digit_window(
    radial_limit: int,
    order_limit: int,
    residue_limit: int,
) -> None:
    """Verify the exact affine ``p*d+r`` score and unit transform."""

    tested_types = 0
    relative_scores: dict[int, int] = defaultdict(int)
    one_high_digit_quotients: dict[tuple[int, int], int] = defaultdict(int)
    first_examples: dict[
        int,
        tuple[
            int,
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
        for high_digit in range(1, order_limit + 1):
            prime = next_prime(radial_degree * (high_digit + 1))
            assert prime > radial_degree * (high_digit + 1)
            for low_digit in range(min(prime - 1, residue_limit) + 1):
                moment_order = prime * high_digit + low_digit
                tables = {
                    support: two_digit_selection_table(
                        support,
                        moment_order,
                        prime,
                        low_digit,
                    )
                    for support in supports
                }

                for operator_support in supports:
                    operator_table = tables[operator_support]
                    for polynomial_support in supports:
                        polynomial_table = tables[polynomial_support]
                        common_radials = (
                            operator_table.keys()
                            & polynomial_table.keys()
                        )
                        for radial_x in common_radials:
                            radial_y = (
                                radial_degree * moment_order - radial_x
                            )
                            radial_remainders = (
                                radial_x % prime,
                                radial_y % prime,
                            )
                            radial_quotients = (
                                radial_x // prime,
                                radial_y // prime,
                            )
                            total_low_residue = (
                                radial_degree * low_digit
                            ) % prime
                            radial_carry = (
                                sum(radial_remainders)
                                - total_low_residue
                            ) // prime
                            assert radial_carry in (0, 1)

                            radial_value = (
                                factorial(radial_x)
                                * factorial(radial_y)
                            )
                            radial_valuation = valuation(
                                radial_value,
                                prime,
                            )
                            predicted_radial_valuation = (
                                radial_degree * high_digit
                                + (
                                    radial_degree * low_digit
                                ) // prime
                                - radial_carry
                            )
                            assert (
                                radial_valuation
                                == predicted_radial_valuation
                            )
                            radial_unit = unit_mod_prime(
                                radial_value,
                                prime,
                                radial_valuation,
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
                                assert sum(
                                    entry // prime for entry in operator.entries
                                ) == high_digit - operator.carry
                                operator_multinomial = multinomial(
                                    moment_order,
                                    operator.entries,
                                )
                                operator_valuation = valuation(
                                    operator_multinomial,
                                    prime,
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
                                    * factorial(high_digit)
                                    * factorial(low_digit)
                                    * pow(
                                        operator_denominator % prime,
                                        -1,
                                        prime,
                                    )
                                ) % prime
                                assert (
                                    operator_unit
                                    == predicted_operator_unit
                                )

                                for polynomial in polynomial_table[radial_x]:
                                    assert sum(
                                        entry // prime
                                        for entry in polynomial.entries
                                    ) == high_digit - polynomial.carry
                                    polynomial_multinomial = multinomial(
                                        moment_order,
                                        polynomial.entries,
                                    )
                                    polynomial_valuation = valuation(
                                        polynomial_multinomial,
                                        prime,
                                    )
                                    assert (
                                        polynomial_valuation
                                        == polynomial.carry
                                    )
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
                                        * factorial(high_digit)
                                        * factorial(low_digit)
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

                                    exact_valuation = (
                                        operator_valuation
                                        + polynomial_valuation
                                        + radial_valuation
                                    )
                                    predicted_valuation = (
                                        radial_degree * high_digit
                                        + (
                                            radial_degree * low_digit
                                        ) // prime
                                        + operator.carry
                                        + polynomial.carry
                                        - radial_carry
                                    )
                                    assert exact_valuation == predicted_valuation
                                    relative_score = (
                                        operator.carry
                                        + polynomial.carry
                                        - radial_carry
                                    )
                                    assert relative_score >= -1
                                    relative_scores[relative_score] += 1
                                    if high_digit == 1:
                                        assert operator.carry in (0, 1)
                                        assert polynomial.carry in (0, 1)
                                        one_high_digit_quotients[
                                            (
                                                1 - operator.carry,
                                                1 - polynomial.carry,
                                            )
                                        ] += 1
                                    first_examples.setdefault(
                                        relative_score,
                                        (
                                            radial_degree,
                                            high_digit,
                                            low_digit,
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
        "PASS two-digit affine jet--carry transform: "
        f"{tested_types} exact return types"
    )
    print(f"two-digit relative score distribution: {dict(relative_scores)}")
    print(
        "p+r high-quotient side counts: "
        f"{dict(one_high_digit_quotients)}"
    )
    for score, example in sorted(first_examples.items()):
        print(f"first relative-score {score} example: {example}")
    print(
        "STATUS: exact bounded regression for the general p*d+r formulas; "
        "factorial-compatible shell inheritance remains open"
    )


def signed_singleton_unit(
    low_correction: tuple[int, ...],
    mark: int,
    low_total: int,
    prime: int,
) -> tuple[int, int]:
    """Return the predicted valuation and unit for ``p*e_mark+b``.

    The off-mark entries of ``b`` are nonnegative, while its marked entry
    may be negative.  This is the signed recentering naturally produced by
    an exposed ``p+r`` singleton fibre.
    """

    marked_correction = low_correction[mark]
    off_factorial = 1
    for index, entry in enumerate(low_correction):
        if index != mark:
            assert entry >= 0
            off_factorial *= factorial(entry)

    if marked_correction >= 0:
        denominator = factorial(marked_correction) * off_factorial
        unit = (
            factorial(low_total)
            * pow(denominator % prime, -1, prime)
        ) % prime
        return 0, unit

    borrow = -marked_correction
    denominator = off_factorial % prime
    unit = (
        (-1) ** (borrow - 1)
        * factorial(low_total)
        * factorial(borrow - 1)
        * pow(denominator, -1, prime)
    ) % prime
    return 1, unit


def verify_exposed_singleton_signed_digits(
    radial_limit: int,
    residue_limit: int,
) -> None:
    """Check stable signed corrections at the two exposed level marks.

    For levels ``0,...,s`` and an endpoint mark ``h``, a state in the
    same mass/level fibre as ``p*e_h+a`` has the form ``p*e_h+b`` with
    ``b`` in a finite set independent of sufficiently large ``p``.  The
    written proof uses an exposing linear form; this routine enumerates
    the endpoint fibres and verifies the signed factorial-unit formula.
    """

    stable_fibres = 0
    tested_states = 0
    tested_units = 0

    for span in range(1, radial_limit + 1):
        levels = tuple(range(span + 1))
        for low_total in range(residue_limit + 1):
            first_prime = next_prime(2 * span + 2 * low_total + 3)
            second_prime = next_prime(first_prime)
            state_tables = {}
            for prime in (first_prime, second_prime):
                by_level: dict[int, list[tuple[int, ...]]] = defaultdict(list)
                for state in compositions(prime + low_total, len(levels)):
                    level = sum(
                        entry * channel
                        for entry, channel in zip(state, levels)
                    )
                    by_level[level].append(state)
                state_tables[prime] = by_level

            for mark in (0, span):
                for low_state in compositions(low_total, len(levels)):
                    low_level = sum(
                        entry * channel
                        for entry, channel in zip(low_state, levels)
                    )
                    corrections = []
                    for prime in (first_prime, second_prime):
                        target_level = prime * mark + low_level
                        fibre_corrections = set()
                        for state in state_tables[prime].get(target_level, ()):
                            correction = list(state)
                            correction[mark] -= prime
                            correction_tuple = tuple(correction)
                            fibre_corrections.add(correction_tuple)

                            assert sum(correction_tuple) == low_total
                            assert sum(
                                entry * channel
                                for entry, channel in zip(
                                    correction_tuple,
                                    levels,
                                )
                            ) == low_level

                            exact = multinomial(prime + low_total, state)
                            predicted_order, predicted_unit = (
                                signed_singleton_unit(
                                    correction_tuple,
                                    mark,
                                    low_total,
                                    prime,
                                )
                            )
                            assert valuation(exact, prime) == predicted_order
                            assert unit_mod_prime(
                                exact,
                                prime,
                                predicted_order,
                            ) == predicted_unit
                            tested_states += 1
                            tested_units += 1
                        corrections.append(fibre_corrections)

                    assert corrections[0] == corrections[1]
                    stable_fibres += 1

    print(
        "PASS exposed-singleton signed-digit stabilization: "
        f"{stable_fibres} endpoint fibres, {tested_states} states, "
        f"{tested_units} p-free units"
    )


def verify_interior_singleton_graver_family(prime_limit: int) -> None:
    """Replay the sharp positive-density centered-triple family.

    At the nonvertex level one, the empty-high state

        ((p-1)/2, 1, (p-1)/2)

    shares the mass and level of ``(0,p,0)``.  Their difference is an
    unbounded multiple of the fixed centered-triple Graver move.  Swapping
    the two states between the operator and polynomial sides preserves the
    complete factorial weight, so side-blind scalar tomography cannot turn
    this state-level fact into packet inheritance by itself.
    """

    tested_primes = 0
    hilbert_states = 0
    primitive = (1, -2, 1)
    pure_generator = (1, 0, 1, 0)
    centered_generator = (2, 1, 0, 1)
    for prime in range(5, prime_limit + 1):
        if not is_prime(prime) or prime == 2:
            continue
        for multiplicity in range(prime // 2 + 1):
            pure_count = prime - 2 * multiplicity
            state = (
                multiplicity,
                pure_count,
                multiplicity,
            )
            module_vector = (
                prime,
                *state,
            )
            reconstructed = tuple(
                pure_count * pure_entry
                + multiplicity * centered_entry
                for pure_entry, centered_entry in zip(
                    pure_generator,
                    centered_generator,
                    strict=True,
                )
            )
            assert module_vector == reconstructed
            assert sum(state) == prime
            assert state[1] + 2 * state[2] == prime
            hilbert_states += 1

        repeats = (prime - 1) // 2
        singleton = (0, prime, 0)
        distributed = (repeats, 1, repeats)
        assert sum(singleton) == sum(distributed) == prime
        assert singleton[1] + 2 * singleton[2] == prime
        assert distributed[1] + 2 * distributed[2] == prime
        assert tuple(
            right - left
            for left, right in zip(singleton, distributed)
        ) == tuple(repeats * entry for entry in primitive)

        singleton_multinomial = multinomial(prime, singleton)
        distributed_multinomial = multinomial(prime, distributed)
        assert valuation(singleton_multinomial, prime) == 0
        assert valuation(distributed_multinomial, prime) == 1

        radial_factorial = factorial(prime) ** 2
        first_packet_weight = (
            singleton_multinomial
            * distributed_multinomial
            * radial_factorial
        )
        second_packet_weight = (
            distributed_multinomial
            * singleton_multinomial
            * radial_factorial
        )
        assert first_packet_weight == second_packet_weight
        tested_primes += 1

    assert tested_primes
    print(
        "PASS interior singleton Graver family: "
        f"{tested_primes} odd primes through {prime_limit}; "
        f"{hilbert_states} states in the Hilbert basis "
        "(1;0,1,0),(2;1,0,1); "
        "difference=((p-1)/2)*(1,-2,1)"
    )


def falling_factorial(value: int, order: int) -> int:
    answer = 1
    for offset in range(order):
        answer *= value - offset
    return answer


def verify_radial_carry_hasse_compression(
    prime_limit: int,
    residue_limit: int,
) -> None:
    """Verify that every long binary radial-carry interval is one row.

    If the two low radial digits add to ``p+R``, their factorial unit is

        t! (p+R-t)! = (-1)^(t-R) t^(falling R+1)  (mod p).

    Consequently a weighted sum over all carry positions is minus the
    ``(R+1)``-st derivative of its generating polynomial at ``-1``.
    """

    tested_positions = 0
    tested_rows = 0
    for prime in range(3, prime_limit + 1):
        if not is_prime(prime):
            continue
        for residue in range(min(residue_limit, prime - 2) + 1):
            order = residue + 1
            coefficients = {
                position: (
                    position**3
                    + 2 * position**2
                    + 3 * position
                    + 5
                )
                for position in range(residue + 1, prime)
            }
            factorial_sum = 0
            derivative = 0
            for position, coefficient in coefficients.items():
                complement = prime + residue - position
                assert 0 <= complement < prime
                factorial_unit = (
                    factorial(position) * factorial(complement)
                ) % prime
                predicted = (
                    (-1) ** (position - residue)
                    * falling_factorial(position, order)
                ) % prime
                assert factorial_unit == predicted
                factorial_sum += coefficient * factorial_unit
                derivative += (
                    coefficient
                    * falling_factorial(position, order)
                    * (-1) ** (position - order)
                )
                tested_positions += 1
            assert factorial_sum % prime == (-derivative) % prime
            tested_rows += 1

    assert tested_rows
    print(
        "PASS binary radial-carry Hasse compression: "
        f"{tested_positions} positions in {tested_rows} rows; "
        "carry interval = -F^(R+1)(-1) mod p"
    )


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
    parser.add_argument("--residue-limit", type=int, default=3)
    arguments = parser.parse_args()
    if arguments.residue_limit < 0:
        parser.error("--residue-limit must be nonnegative")
    verify_primitive_bridge_support(arguments.bridge_limit)
    verify_nonhomogeneous_jet_score()
    verify_nonhomogeneous_two_digit_score(arguments.residue_limit)
    verify_window(arguments.radial_limit, arguments.order_limit)
    verify_two_digit_window(
        arguments.radial_limit,
        arguments.order_limit,
        arguments.residue_limit,
    )
    verify_exposed_singleton_signed_digits(
        arguments.radial_limit,
        arguments.residue_limit,
    )
    verify_interior_singleton_graver_family(arguments.bridge_limit)
    verify_radial_carry_hasse_compression(
        arguments.bridge_limit,
        arguments.residue_limit,
    )


if __name__ == "__main__":
    main()

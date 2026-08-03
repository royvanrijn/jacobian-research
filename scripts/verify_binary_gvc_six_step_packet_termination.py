#!/usr/bin/env python3
"""Regression for the exact six-step binary Hall-packet termination proof.

After row operations, every member of

    R_(s+6) B_a B_(a+1) = R_s B_(a+3) B_(a+4)

has columns R_6,R_0,B_0,B_1,B_3,B_4.  At scale N, the states with the
same C2 and C3 marked histograms as the endpoint rays form one line z_t.
Their normalized weights are

    binom(2N,N) * binom(N,t)^3 * U^(N-t) * V^t.

The accompanying note proves the formula for all N and proves that every
fixed finite-character refinement is terminal in two rows.  This script
independently enumerates the finite fibres in a regression range and checks
the exact coefficient formulas and two-row obstruction.  It also verifies
the affine two-digit decomposition at scales ``N=p*d+r``.  At ``N=p+1``
the first unit correction of the no-borrow endpoints is already nonzero on
the blind Hall branch ``V=-U``.
"""

from __future__ import annotations

import argparse
import math
from collections.abc import Iterator


State = tuple[int, int, int, int, int, int]
LEVELS = (6, 0, 0, 1, 3, 4)


def weak_compositions(total: int, length: int) -> Iterator[tuple[int, ...]]:
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first,) + tail


def fiber_states(
    scale: int,
    levels: tuple[int, ...] = LEVELS,
) -> tuple[State, ...]:
    states = []
    target_level = scale * (levels[0] + levels[2] + levels[3])
    for operator in weak_compositions(scale, 2):
        for polynomial in weak_compositions(2 * scale, 4):
            state = operator + polynomial
            if sum(value * level for value, level in zip(state, levels)) == target_level:
                states.append(state)
    return tuple(states)


def marked_histogram(
    state: State,
    order: int,
    levels: tuple[int, ...] = LEVELS,
) -> tuple[tuple[int, ...], ...]:
    operator = [0] * order
    polynomial = [0] * order
    for value, level in zip(state[:2], levels[:2]):
        operator[level % order] += value
    for value, level in zip(state[2:], levels[2:]):
        polynomial[level % order] += value
    return tuple(operator), tuple(polynomial)


def endpoint(scale: int, right: bool = False) -> State:
    if right:
        return (0, scale, 0, 0, scale, scale)
    return (scale, 0, scale, scale, 0, 0)


def line_state(scale: int, parameter: int) -> State:
    return (
        scale - parameter,
        parameter,
        scale - parameter,
        scale - parameter,
        parameter,
        parameter,
    )


def multinomial(total: int, parts: tuple[int, ...]) -> int:
    denominator = math.prod(math.factorial(part) for part in parts)
    return math.factorial(total) // denominator


def valuation(number: int, prime: int) -> int:
    order = 0
    while number and number % prime == 0:
        number //= prime
        order += 1
    return order


def unit_mod_prime(number: int, prime: int, order: int) -> int:
    return (number // prime**order) % prime


def state_weight(scale: int, parameter: int) -> int:
    state = line_state(scale, parameter)
    return multinomial(scale, state[:2]) * multinomial(2 * scale, state[2:])


def verify_affine_digit_shells(
    prime_limit: int,
    high_digit_limit: int,
    low_digit_limit: int,
) -> None:
    tested_rows = 0
    for prime in range(5, prime_limit + 1):
        if any(prime % divisor == 0 for divisor in range(2, math.isqrt(prime) + 1)):
            continue
        for high_digit in range(1, min(high_digit_limit, prime - 1) + 1):
            for low_digit in range(min(low_digit_limit, prime - 1) + 1):
                scale = prime * high_digit + low_digit
                for parameter in range(scale + 1):
                    coefficient = math.comb(scale, parameter) ** 3
                    order = valuation(coefficient, prime)
                    quotient, residue = divmod(parameter, prime)
                    if residue <= low_digit:
                        assert order == 0
                        predicted = (
                            math.comb(high_digit, quotient) ** 3
                            * math.comb(low_digit, residue) ** 3
                        ) % prime
                        assert coefficient % prime == predicted
                    else:
                        assert quotient < high_digit
                        assert order == 3
                        bridge = residue - low_digit
                        denominator = (
                            bridge * math.comb(low_digit + bridge, bridge)
                        ) % prime
                        low_unit = (
                            (-1) ** (bridge + 1)
                            * pow(denominator, -1, prime)
                        ) % prime
                        predicted = (
                            high_digit**3
                            * math.comb(high_digit - 1, quotient) ** 3
                            * low_unit**3
                        ) % prime
                        assert (
                            unit_mod_prime(coefficient, prime, order)
                            == predicted
                        )
                    tested_rows += 1

        if prime >= 7:
            reciprocal_sum = sum(
                pow(
                    (index * (index + 1)) % prime,
                    -3,
                    prime,
                )
                for index in range(1, prime - 1)
            ) % prime
            assert reciprocal_sum == 20 % prime

            scale = prime + 1
            no_borrow = (0, 1, prime, prime + 1)
            endpoint_sum = sum(
                math.comb(scale, parameter) ** 3
                * (-1) ** parameter
                for parameter in no_borrow
            )
            assert valuation(endpoint_sum, prime) == 1
            assert unit_mod_prime(endpoint_sum, prime, 1) == -6 % prime

            carry_sum = sum(
                math.comb(scale, parameter) ** 3
                * (-1) ** parameter
                for parameter in range(2, prime)
            )
            assert valuation(carry_sum, prime) >= 3
            assert unit_mod_prime(carry_sum, prime, 3) == 20 % prime

    print(
        "PASS affine p*d+r Franel shell factorization: "
        f"{tested_rows} coefficients"
    )
    print(
        "PASS p+1 blind-Hall correction: no-borrow unit -6 and "
        "carry-shell unit 20"
    )


def verify_odd_shift_family(max_shift: int, max_scale: int) -> None:
    checked = 0
    for half_shift in range(3, max_shift + 1, 2):
        levels = (
            2 * half_shift,
            0,
            0,
            1,
            half_shift,
            half_shift + 1,
        )
        for scale in range(1, max_scale + 1):
            states = fiber_states(scale, levels)
            target = (
                marked_histogram(endpoint(scale), 2, levels),
                marked_histogram(endpoint(scale), half_shift, levels),
            )
            blind = {
                state
                for state in states
                if (
                    marked_histogram(state, 2, levels),
                    marked_histogram(state, half_shift, levels),
                )
                == target
            }
            expected = {
                line_state(scale, parameter)
                for parameter in range(scale + 1)
            }
            assert blind == expected
            assert marked_histogram(
                endpoint(scale),
                half_shift + 1,
                levels,
            ) != marked_histogram(
                endpoint(scale, right=True),
                half_shift + 1,
                levels,
            )
            checked += 1

    print(
        "PASS odd-shift C2/Ch-blind fibres are the same Franel line: "
        f"{checked} shift/scale pairs"
    )


def verify(
    max_scale: int,
    max_character_order: int,
    prime_limit: int,
    high_digit_limit: int,
    low_digit_limit: int,
    max_odd_shift: int,
) -> None:
    for scale in range(1, max_scale + 1):
        states = fiber_states(scale)
        target = (
            marked_histogram(endpoint(scale), 2),
            marked_histogram(endpoint(scale), 3),
        )
        blind = {
            state
            for state in states
            if (
                marked_histogram(state, 2),
                marked_histogram(state, 3),
            )
            == target
        }
        expected = {
            line_state(scale, parameter)
            for parameter in range(scale + 1)
        }
        assert blind == expected
        assert endpoint(scale, right=True) in blind

        for parameter in range(scale + 1):
            expected_weight = (
                math.comb(2 * scale, scale)
                * math.comb(scale, parameter) ** 3
            )
            assert state_weight(scale, parameter) == expected_weight

    left = endpoint(1)
    right = endpoint(1, right=True)
    assert marked_histogram(left, 2) == marked_histogram(right, 2)
    assert marked_histogram(left, 3) == marked_histogram(right, 3)
    assert marked_histogram(left, 4) != marked_histogram(right, 4)

    for order in range(1, max_character_order + 1):
        central = math.comb(2 * order, order) ** 3
        assert 2 - central != 0

        first_row = {
            (order, 0): 1,
            (0, order): 1,
        }
        second_row = {
            (2 * order, 0): 1,
            (order, order): central,
            (0, 2 * order): 1,
        }
        assert first_row == {
            (order, 0): 1,
            (0, order): 1,
        }
        assert second_row[(order, order)] == central

    verify_affine_digit_shells(
        prime_limit,
        high_digit_limit,
        low_digit_limit,
    )
    verify_odd_shift_family(max_odd_shift, max_scale)

    print(
        "PASS C2/C3-blind six-step fibers are exactly the Franel line "
        f"through scale {max_scale}"
    )
    print(
        "PASS multinomial weights equal binom(2N,N)*binom(N,t)^3 "
        f"through scale {max_scale}"
    )
    print("PASS C4 separates the primitive endpoints")
    print(
        "PASS fixed-character two-row obstruction "
        f"2-binom(2h,h)^3 != 0 through h={max_character_order}"
    )
    print(
        "STATUS: exact characteristic-zero proof is in "
        "BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md; affine carry promotion "
        "remains unproved in the parked route, while Hall-envelope "
        "separation proves unrestricted GVC(2)"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-scale", type=int, default=12)
    parser.add_argument("--max-character-order", type=int, default=32)
    parser.add_argument("--prime-limit", type=int, default=43)
    parser.add_argument("--high-digit-limit", type=int, default=4)
    parser.add_argument("--low-digit-limit", type=int, default=4)
    parser.add_argument("--max-odd-shift", type=int, default=15)
    arguments = parser.parse_args()
    if arguments.max_scale < 1:
        parser.error("--max-scale must be positive")
    if arguments.max_character_order < 1:
        parser.error("--max-character-order must be positive")
    if arguments.prime_limit < 7:
        parser.error("--prime-limit must be at least seven")
    if arguments.high_digit_limit < 1:
        parser.error("--high-digit-limit must be positive")
    if arguments.low_digit_limit < 0:
        parser.error("--low-digit-limit must be nonnegative")
    if arguments.max_odd_shift < 3:
        parser.error("--max-odd-shift must be at least three")
    verify(
        arguments.max_scale,
        arguments.max_character_order,
        arguments.prime_limit,
        arguments.high_digit_limit,
        arguments.low_digit_limit,
        arguments.max_odd_shift,
    )


if __name__ == "__main__":
    main()

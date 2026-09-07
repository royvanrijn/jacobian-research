#!/usr/bin/env python3
"""Exact replay of the Cobham obstruction for Hall carry states.

For an odd prime ``p`` let ``carry_p(n)`` be the indicator that adding
``n+n`` in base ``p`` produces at least one carry.  By Kummer this is the
indicator that ``p`` divides ``binomial(2*n, n)``.  It is a two-state
``p``-automatic sequence, but it is not ultimately periodic.

Consequently no finite sequence automatic in two multiplicatively
independent bases can refine the complete carry states for either base:
Cobham would make the refinement, and hence ``carry_p``, ultimately
periodic.  The calculation below checks the digit/Kummer formulas, explicit
witnesses against every proposed bounded period, and the loss of information
on sparse prime-power rays.

This is an exact regression for the elementary proof in
``extended-geometry/BINARY_GVC_PRIMITIVE_TRANSLATION_OBSERVABILITY.md``.
It is an obstruction to a proposed proof route, not a GVC counterexample.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable


DEFAULT_PRIMES = (3, 5, 7, 11, 13)


def is_prime(number: int) -> bool:
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def base_digits(number: int, base: int) -> tuple[int, ...]:
    """Return the least-significant-first base-``base`` digits."""

    assert number >= 0
    assert base >= 2
    if number == 0:
        return (0,)
    digits: list[int] = []
    while number:
        number, digit = divmod(number, base)
        digits.append(digit)
    return tuple(digits)


def equal_addition_carry_positions(number: int, base: int) -> tuple[int, ...]:
    """Positions carrying out while adding ``number + number`` in a base."""

    incoming = 0
    positions: list[int] = []
    for position, digit in enumerate(base_digits(number, base)):
        outgoing = (2 * digit + incoming) // base
        assert outgoing in (0, 1)
        if outgoing:
            positions.append(position)
        incoming = outgoing
    return tuple(positions)


def factorial_valuation(number: int, prime: int) -> int:
    """Legendre valuation of ``number!``."""

    value = 0
    while number:
        number //= prime
        value += number
    return value


def central_binomial_valuation(number: int, prime: int) -> int:
    """The ``prime``-valuation of ``binomial(2*number, number)``."""

    return factorial_valuation(2 * number, prime) - 2 * factorial_valuation(
        number, prime
    )


def carry_indicator(number: int, prime: int) -> int:
    """Indicator of at least one carry in ``number + number`` base ``prime``."""

    return int(bool(equal_addition_carry_positions(number, prime)))


def digit_automaton_indicator(number: int, prime: int) -> int:
    """Output of the two-state DFAO recognizing a central-binomial carry."""

    threshold = (prime - 1) // 2
    state = 0
    for digit in base_digits(number, prime):
        state = int(state or digit > threshold)
    return state


def prime_valuation(number: int, prime: int) -> int:
    """Return ``v_prime(number)`` for a positive integer."""

    assert number > 0
    value = 0
    while number % prime == 0:
        number //= prime
        value += 1
    return value


def nonperiodicity_witness(
    prime: int,
    period: int,
    minimum_index: int,
) -> dict[str, int]:
    """Construct equal-residue, unequal-output terms beyond an index.

    Write ``period = prime**a * unit`` and choose ``m`` with
    ``m*unit == -1 (mod prime)``.  The base-prime digit of ``m*period`` at
    position ``a`` is then ``prime-1``.  A sufficiently remote ``prime**k``
    has no carry, while ``prime**k + m*period`` does.
    """

    assert prime % 2 == 1
    assert period > 0
    exponent = prime_valuation(period, prime)
    unit = period // prime**exponent
    multiplier = ((prime - 1) * pow(unit, -1, prime)) % prime
    assert 1 <= multiplier < prime
    displacement = multiplier * period

    power = 1
    power_exponent = 0
    lower_bound = max(minimum_index, displacement)
    while power <= lower_bound:
        power *= prime
        power_exponent += 1

    first = power
    second = power + displacement
    assert first >= minimum_index
    assert second >= minimum_index
    assert (second - first) % period == 0
    assert carry_indicator(first, prime) == 0
    assert carry_indicator(second, prime) == 1
    return {
        "period": period,
        "valuation_of_period": exponent,
        "multiplier": multiplier,
        "first_index": first,
        "second_index": second,
        "first_output": 0,
        "second_output": 1,
        "power_exponent": power_exponent,
    }


def digest_rows(rows: Iterable[dict[str, int]]) -> str:
    payload = "\n".join(
        json.dumps(row, sort_keys=True, separators=(",", ":")) for row in rows
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def verify_prime(
    prime: int,
    maximum_index: int,
    maximum_period: int,
    minimum_periodic_index: int,
    maximum_ray_coefficient: int,
    maximum_ray_exponent: int,
) -> dict[str, object]:
    assert prime % 2 == 1

    carry_histogram = {0: 0, 1: 0}
    maximum_carry_count = 0
    for number in range(maximum_index + 1):
        positions = equal_addition_carry_positions(number, prime)
        valuation = central_binomial_valuation(number, prime)
        assert valuation == len(positions), (prime, number, valuation, positions)
        indicator = int(bool(positions))
        assert indicator == digit_automaton_indicator(number, prime)
        carry_histogram[indicator] += 1
        maximum_carry_count = max(maximum_carry_count, len(positions))

    witnesses = [
        nonperiodicity_witness(prime, period, minimum_periodic_index)
        for period in range(1, maximum_period + 1)
    ]

    ray_checks = 0
    for coefficient in range(maximum_ray_coefficient + 1):
        expected = carry_indicator(coefficient, prime)
        for exponent in range(1, maximum_ray_exponent + 1):
            assert carry_indicator(coefficient * prime**exponent, prime) == expected
            ray_checks += 1

    sample_periods = sorted(
        {1, min(2, maximum_period), min(prime, maximum_period), maximum_period}
    )
    return {
        "dfa": {
            "states": 2,
            "transition": "seen_high_digit OR digit>(p-1)/2",
            "accepting_output": 1,
        },
        "index_census": {
            "maximum_index": maximum_index,
            "output_histogram": {str(key): value for key, value in carry_histogram.items()},
            "maximum_carry_count": maximum_carry_count,
            "kummer_legendre_checks": maximum_index + 1,
        },
        "nonperiodicity": {
            "periods_checked": maximum_period,
            "minimum_index": minimum_periodic_index,
            "witness_digest": digest_rows(witnesses),
            "sample_witnesses": [witnesses[period - 1] for period in sample_periods],
            "largest_second_index": max(row["second_index"] for row in witnesses),
        },
        "prime_power_rays": {
            "maximum_coefficient": maximum_ray_coefficient,
            "maximum_exponent": maximum_ray_exponent,
            "checks": ray_checks,
            "identity": "carry_p(q*p^e)=carry_p(q)",
        },
    }


def first_disagreement(first_prime: int, second_prime: int) -> dict[str, int]:
    assert first_prime < second_prime
    number = (first_prime + 1) // 2
    assert carry_indicator(number, first_prime) == 1
    assert carry_indicator(number, second_prime) == 0
    for earlier in range(number):
        assert carry_indicator(earlier, first_prime) == carry_indicator(
            earlier, second_prime
        )
    return {
        "first_prime": first_prime,
        "second_prime": second_prime,
        "first_disagreement": number,
        "first_output": 1,
        "second_output": 0,
    }


def build_certificate(args: argparse.Namespace) -> dict[str, object]:
    primes = tuple(sorted(set(args.primes)))
    assert primes and all(prime >= 3 and is_prime(prime) for prime in primes)

    prime_certificates = {
        str(prime): verify_prime(
            prime,
            args.maximum_index,
            args.maximum_period,
            args.minimum_periodic_index,
            args.maximum_ray_coefficient,
            args.maximum_ray_exponent,
        )
        for prime in primes
    }
    disagreements = [
        first_disagreement(first, second)
        for index, first in enumerate(primes)
        for second in primes[index + 1 :]
    ]

    return {
        "theorem": "prime-specific central-binomial carry obstruction to Cobham promotion",
        "status": "proved theorem and exact bounded regression; not a GVC(2) counterexample",
        "parameters": {
            "primes": list(primes),
            "maximum_index": args.maximum_index,
            "maximum_period": args.maximum_period,
            "minimum_periodic_index": args.minimum_periodic_index,
            "maximum_ray_coefficient": args.maximum_ray_coefficient,
            "maximum_ray_exponent": args.maximum_ray_exponent,
        },
        "prime_certificates": prime_certificates,
        "cross_prime_first_disagreements": disagreements,
        "proved_implication": (
            "a finite sequence automatic in multiplicatively independent bases "
            "cannot refine either complete carry state, since a letter image "
            "would be this non-ultimately-periodic indicator"
        ),
        "sparse_probe_obstruction": (
            "every p-power ray is stationary after digit shift although the "
            "full-scale carry sequence is not ultimately periodic"
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--primes",
        type=lambda value: tuple(int(item) for item in value.split(",")),
        default=DEFAULT_PRIMES,
    )
    parser.add_argument("--maximum-index", type=int, default=10000)
    parser.add_argument("--maximum-period", type=int, default=512)
    parser.add_argument("--minimum-periodic-index", type=int, default=10000)
    parser.add_argument("--maximum-ray-coefficient", type=int, default=32)
    parser.add_argument("--maximum-ray-exponent", type=int, default=8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    assert args.maximum_index >= 0
    assert args.maximum_period >= 1
    assert args.minimum_periodic_index >= 0
    assert args.maximum_ray_coefficient >= 0
    assert args.maximum_ray_exponent >= 1
    certificate = build_certificate(args)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")
    print(
        "PASS Cobham carry obstruction: exact Kummer/DFA formulas, "
        f"periods 1..{args.maximum_period}, and stationary p-power rays"
    )


if __name__ == "__main__":
    main()

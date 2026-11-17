#!/usr/bin/env python3
"""Exact regressions for first-ghost source collapse and ray rigidity."""

from __future__ import annotations

from collections import Counter
from math import comb, factorial, gcd


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def multinomial(total: int, entries: tuple[int, ...]) -> int:
    value = factorial(total)
    for entry in entries:
        value //= factorial(entry)
    return value


def verify_source_collapse() -> tuple[int, int]:
    checks = 0
    fibres = 0
    for prime in (5, 7, 11):
        levels = (0, 1, 2, 3)
        by_total: dict[int, list[tuple[int, ...]]] = {}
        for state in compositions(prime, len(levels)):
            if max(state) == prime:
                continue
            total = sum(n * level for n, level in zip(state, levels))
            by_total.setdefault(total, []).append(state)
            unit = factorial(prime - 1)
            for n in state:
                unit //= factorial(n)
            assert multinomial(prime, state) == prime * unit
            for order in range(prime):
                direct = (multinomial(prime, state) // prime) * comb(total, order)
                predicted = unit * comb(total, order)
                assert direct == predicted
                checks += 1
        for total, states in by_total.items():
            if len(states) <= 1:
                continue
            fibres += 1
            signatures = set()
            for state in states:
                unit = factorial(prime - 1)
                for n in state:
                    unit //= factorial(n)
                signature = tuple(
                    ((multinomial(prime, state) // prime) * comb(total, order)) // unit
                    for order in range(prime)
                )
                signatures.add(signature)
            assert len(signatures) == 1
    return checks, fibres


def verify_centered_witness() -> int:
    checked = 0
    for prime in (5, 7, 11, 13, 17, 19):
        left = (0, 0, 2, prime - 2)
        right = (0, 1, 0, prime - 1)
        levels = (0, 1, 2, 3)
        assert sum(left) == sum(right) == prime
        total_left = sum(n * level for n, level in zip(left, levels))
        total_right = sum(n * level for n, level in zip(right, levels))
        assert total_left == total_right == 3 * prime - 2
        assert tuple(a - b for a, b in zip(left, right)) == (0, -1, 2, -1)
        for order in range(prime):
            assert comb(total_left, order) == comb(total_right, order)
        checked += 1
    return checked


def verify_vandermonde() -> int:
    checked = 0
    for left in range(8):
        for right in range(8):
            for order in range(1, left + right + 1):
                two_block = sum(
                    comb(left, first) * comb(right, order - first)
                    for first in range(1, order)
                    if first <= left and order - first <= right
                )
                endpoints = comb(left, order) + comb(right, order)
                assert endpoints + two_block == comb(left + right, order)
                checked += 1
    return checked


def partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    maximum = total if maximum is None else min(maximum, total)
    for first in range(maximum, 0, -1):
        for tail in partitions(total - first, first):
            yield (first,) + tail


def root_signature(parts: tuple[int, ...]) -> Counter[tuple[int, int]]:
    signature: Counter[tuple[int, int]] = Counter()
    for part in parts:
        for numerator in range(1, part + 1):
            divisor = gcd(numerator, part)
            signature[(numerator // divisor, part // divisor)] += 1
    return signature


def recover(signature: Counter[tuple[int, int]]) -> tuple[int, ...]:
    remaining = signature.copy()
    result: list[int] = []
    while remaining:
        largest = max(
            denominator
            for (numerator, denominator), multiplicity in remaining.items()
            if numerator == 1 and multiplicity
        )
        multiplicity = remaining[(1, largest)]
        result.extend([largest] * multiplicity)
        for numerator in range(1, largest + 1):
            divisor = gcd(numerator, largest)
            key = (numerator // divisor, largest // divisor)
            remaining[key] -= multiplicity
            assert remaining[key] >= 0
            if remaining[key] == 0:
                del remaining[key]
    return tuple(sorted(result, reverse=True))


def scaled_factorial_product(parts: tuple[int, ...], scale: int) -> int:
    value = 1
    for part in parts:
        value *= factorial(scale * part)
    return value


def verify_ray_rigidity() -> int:
    checked = 0
    for total in range(1, 23):
        seen: set[tuple[tuple[tuple[int, int], int], ...]] = set()
        for part in partitions(total):
            signature = root_signature(part)
            assert recover(signature) == part
            key = tuple(sorted(signature.items()))
            assert key not in seen
            seen.add(key)
            checked += 1

    for left, right in (
        ((4, 1, 1, 1), (3, 2, 2)),
        ((6, 1, 1), (5, 3)),
        ((10, 1, 1, 1), (7, 6)),
    ):
        assert scaled_factorial_product(left, 1) == scaled_factorial_product(right, 1)
        assert root_signature(left) != root_signature(right)
        assert any(
            scaled_factorial_product(left, scale)
            != scaled_factorial_product(right, scale)
            for scale in range(2, 7)
        )
    return checked


def main() -> None:
    rows, fibres = verify_source_collapse()
    witnesses = verify_centered_witness()
    vandermonde = verify_vandermonde()
    partitions_checked = verify_ray_rigidity()
    print(
        "PASS first-ghost source-total collapse: "
        f"{rows} normalized Hasse rows across {fibres} nontrivial fibres"
    )
    print(f"PASS uniform invisible centered-triple family: {witnesses} primes")
    print(f"PASS Bell/Vandermonde recombination: {vandermonde} identities")
    print(
        "PASS scaled-factorial root reconstruction: "
        f"{partitions_checked} partitions through total 22"
    )
    print("STATUS: local first-ghost exposure is false; pure rays are rigid")


if __name__ == "__main__":
    main()

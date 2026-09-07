#!/usr/bin/env python3
"""Exact small-support searches at the binary torsion--torus frontier.

This script has two deliberately separate modes.

``pair``
    Search for Laurent polynomials ``f,g`` such that

        CT(f**n) + CT(g**n) = 0

    through a prescribed depth, while the origin belongs to both Newton
    intervals.  This tests the unrestricted ``Z x C_2`` trace statement.

``shared``
    Search for one Laurent polynomial ``H`` and two distinct target slopes
    ``r,s`` such that

        [z**(n*r)] H**n + [z**(n*s)] H**n = 0.

    This is the one-free-variable shared-factor specialization of the
    structured Cartesian trace packet.

A bounded search is not a proof.  Its purpose is to find a counterexample
candidate or identify the first exact obstruction row for the next
symbolic saturation.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import defaultdict
from collections.abc import Iterable


Laurent = dict[int, int]


def multiply(left: Laurent, right: Laurent) -> Laurent:
    answer: Laurent = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = left_exponent + right_exponent
            answer[exponent] = (
                answer.get(exponent, 0) + left_value * right_value
            )
    return {
        exponent: value
        for exponent, value in answer.items()
        if value
    }


def coefficient_rows(
    polynomial: Laurent,
    targets: Iterable[int],
    depth: int,
) -> dict[int, tuple[int, ...]]:
    targets = tuple(targets)
    rows = {target: [] for target in targets}
    current: Laurent = {0: 1}
    for moment in range(1, depth + 1):
        current = multiply(current, polynomial)
        for target in targets:
            rows[target].append(current.get(moment * target, 0))
    return {
        target: tuple(values)
        for target, values in rows.items()
    }


def coefficient_tuples(width: int, height: int) -> Iterable[tuple[int, ...]]:
    values = range(-height, height + 1)
    for coefficients in itertools.product(values, repeat=2 * width + 1):
        if any(coefficients):
            yield coefficients


def tuple_to_polynomial(coefficients: tuple[int, ...], width: int) -> Laurent:
    return {
        exponent: coefficient
        for exponent, coefficient in zip(
            range(-width, width + 1),
            coefficients,
            strict=True,
        )
        if coefficient
    }


def contains_origin(polynomial: Laurent) -> bool:
    exponents = polynomial.keys()
    return min(exponents) <= 0 <= max(exponents)


def primitive_tuple(coefficients: tuple[int, ...]) -> bool:
    common = 0
    for coefficient in coefficients:
        common = math.gcd(common, abs(coefficient))
    return common == 1


def first_nonzero(values: tuple[int, ...]) -> int | None:
    for index, value in enumerate(values, 1):
        if value:
            return index
    return None


def search_pair(width: int, height: int, depth: int) -> None:
    vectors: dict[tuple[int, ...], tuple[int, ...]] = {}
    best_depth = 0
    best_data = None
    count = 0

    for coefficients in coefficient_tuples(width, height):
        if not primitive_tuple(coefficients):
            continue
        polynomial = tuple_to_polynomial(coefficients, width)
        if not contains_origin(polynomial):
            continue
        rows = coefficient_rows(polynomial, (0,), depth)[0]
        count += 1

        opposite = tuple(-value for value in rows)
        if opposite in vectors:
            print(
                "SURVIVOR unrestricted C2 trace through requested depth"
            )
            print(f"f={vectors[opposite]}")
            print(f"g={coefficients}")
            print(f"rows={rows}")
            return

        for other_rows, other_coefficients in vectors.items():
            summed = tuple(
                left + right
                for left, right in zip(rows, other_rows, strict=True)
            )
            obstruction = first_nonzero(summed)
            survived = depth if obstruction is None else obstruction - 1
            if survived > best_depth:
                best_depth = survived
                best_data = (
                    other_coefficients,
                    coefficients,
                    summed,
                    obstruction,
                )
        vectors.setdefault(rows, coefficients)

    print(f"NO survivor in this torsion--torus search through depth {depth}")
    print(f"primitive origin-containing polynomials={count}")
    print(f"best_zero_prefix={best_depth}")
    if best_data is not None:
        left, right, summed, obstruction = best_data
        print(f"best_f={left}")
        print(f"best_g={right}")
        print(f"summed_rows={summed}")
        print(f"first_obstruction={obstruction}")


def search_shared(width: int, height: int, depth: int) -> None:
    best_depth = 0
    best_data = None
    count = 0

    for coefficients in coefficient_tuples(width, height):
        if not primitive_tuple(coefficients):
            continue
        polynomial = tuple_to_polynomial(coefficients, width)
        active_min = min(polynomial)
        active_max = max(polynomial)
        slopes = tuple(range(active_min, active_max + 1))
        if len(slopes) < 2:
            continue
        rows = coefficient_rows(polynomial, slopes, depth)
        count += 1

        for left_index, left in enumerate(slopes):
            for right in slopes[left_index + 1 :]:
                summed = tuple(
                    left_value + right_value
                    for left_value, right_value in zip(
                        rows[left],
                        rows[right],
                        strict=True,
                    )
                )
                obstruction = first_nonzero(summed)
                survived = (
                    depth if obstruction is None else obstruction - 1
                )
                if survived > best_depth:
                    best_depth = survived
                    best_data = (
                        coefficients,
                        left,
                        right,
                        summed,
                        obstruction,
                    )
                if obstruction is None:
                    print(
                        "SURVIVOR shared-factor trace through requested depth"
                    )
                    print(f"H={coefficients}")
                    print(f"slopes=({left},{right})")
                    print(f"rows={summed}")
                    return

    print(f"NO shared-factor survivor through depth {depth}")
    print(f"primitive polynomials with two slopes={count}")
    print(f"best_zero_prefix={best_depth}")
    if best_data is not None:
        coefficients, left, right, summed, obstruction = best_data
        print(f"best_H={coefficients}")
        print(f"best_slopes=({left},{right})")
        print(f"summed_rows={summed}")
        print(f"first_obstruction={obstruction}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("pair", "shared"))
    parser.add_argument("--width", type=int, default=2)
    parser.add_argument("--height", type=int, default=1)
    parser.add_argument("--depth", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.width < 1 or args.height < 1 or args.depth < 1:
        raise SystemExit("width, height, and depth must be positive")
    if args.mode == "pair":
        search_pair(args.width, args.height, args.depth)
    else:
        search_shared(args.width, args.height, args.depth)


if __name__ == "__main__":
    main()

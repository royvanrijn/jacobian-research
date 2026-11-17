#!/usr/bin/env python3
"""Degree-pattern census for higher Hessian--Ritt intersections.

The canonical Hessian atlas identifies C_(a,b) with the normalized
degree-(a,b) polynomial-composition locus.  This script enumerates all ordered
multiplicative degree words.  A word carrying both cuts is a common
refinement, while swapping adjacent coprime blocks costs one
(one power/Chebyshev Ritt core).

For each pair of factor loci it reports the least number of cores required.
It also prints the prime-refinement shadow, which shows where compatibility
between several cores can first enter higher intersections.  This is a
combinatorial candidate census, not a primary-decomposition certificate.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import math
from functools import reduce


DEFAULT_DEGREES = (18, 24, 30)


def prime_factors(n: int) -> tuple[int, ...]:
    factors: list[int] = []
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            factors.append(divisor)
            n //= divisor
        divisor += 1
    if n > 1:
        factors.append(n)
    return tuple(factors)


def distinct_permutations(values: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(set(itertools.permutations(values))))


def ordered_factorizations(n: int) -> tuple[tuple[int, ...], ...]:
    """All ordered multiplicative words of product ``n`` and length >= 2."""

    result: set[tuple[int, ...]] = set()

    def extend(prefix: tuple[int, ...], remainder: int) -> None:
        if prefix and remainder >= 2:
            result.add(prefix + (remainder,))
        for divisor in range(2, remainder):
            if remainder % divisor == 0:
                extend(prefix + (divisor,), remainder // divisor)

    extend((), n)
    return tuple(sorted(result, key=lambda word: (len(word), word)))


def proper_factor_pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((a, n // a) for a in range(2, n) if n % a == 0)


def product(values: tuple[int, ...]) -> int:
    return reduce(lambda left, right: left * right, values, 1)


def cuts(word: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    n = product(word)
    return frozenset(
        (product(word[:index]), n // product(word[:index]))
        for index in range(1, len(word))
    )


def ritt_neighbors(word: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    result = []
    for index in range(len(word) - 1):
        if math.gcd(word[index], word[index + 1]) != 1:
            continue
        swapped = list(word)
        swapped[index], swapped[index + 1] = swapped[index + 1], swapped[index]
        result.append(tuple(swapped))
    return tuple(result)


def shortest_path(
    universe: frozenset[tuple[int, ...]],
    starts: frozenset[tuple[int, ...]],
    targets: frozenset[tuple[int, ...]],
) -> tuple[int, tuple[tuple[int, ...], ...]] | None:
    queue: list[tuple[int, int, tuple[int, ...]]] = []
    predecessor: dict[tuple[int, ...], tuple[int, ...] | None] = {}
    distance: dict[tuple[int, ...], int] = {}
    for word in starts:
        distance[word] = 0
        predecessor[word] = None
        heapq.heappush(queue, (0, len(word), word))

    endpoint = None
    while queue:
        cost, _, word = heapq.heappop(queue)
        if cost != distance[word]:
            continue
        if word in targets:
            endpoint = word
            break
        for candidate in ritt_neighbors(word):
            if candidate not in universe:
                continue
            new_cost = cost + 1
            if new_cost >= distance.get(candidate, 10**9):
                continue
            distance[candidate] = new_cost
            predecessor[candidate] = word
            heapq.heappush(queue, (new_cost, len(candidate), candidate))
    if endpoint is None:
        return None

    path = [endpoint]
    while predecessor[path[-1]] is not None:
        path.append(predecessor[path[-1]])
    return distance[endpoint], tuple(reversed(path))


def swap_label(left: tuple[int, ...], right: tuple[int, ...]) -> str | None:
    if len(left) != len(right) or sorted(left) != sorted(right):
        return None
    changed = [index for index in range(len(left)) if left[index] != right[index]]
    if len(changed) != 2 or changed[1] != changed[0] + 1:
        return None
    index = changed[0]
    prefix = product(left[:index])
    suffix = product(left[index + 2 :])
    core = f"{left[index]}x{left[index + 1]}"
    decorations = []
    if prefix != 1:
        decorations.append(f"left {prefix}")
    if suffix != 1:
        decorations.append(f"right {suffix}")
    return core if not decorations else f"{core} ({', '.join(decorations)})"


def pair_census(n: int) -> list[dict[str, object]]:
    words = ordered_factorizations(n)
    universe = frozenset(words)
    pairs = proper_factor_pairs(n)
    carriers = {
        pair: frozenset(word for word in words if pair in cuts(word)) for pair in pairs
    }
    rows: list[dict[str, object]] = []
    for first, second in itertools.combinations(pairs, 2):
        common = tuple(sorted(carriers[first] & carriers[second]))
        result = shortest_path(universe, carriers[first], carriers[second])
        assert result is not None
        distance, path = result
        cores = tuple(
            label
            for left, right in zip(path, path[1:])
            if (label := swap_label(left, right)) is not None
        )
        rows.append(
            {
                "loci": (first, second),
                "common": common,
                "distance": distance,
                "path": path,
                "cores": cores,
            }
        )
    return rows


def format_pair(pair: tuple[int, int]) -> str:
    return f"{pair[0]}o{pair[1]}"


def report(n: int) -> None:
    factors = prime_factors(n)
    prime_words = distinct_permutations(factors)
    words = ordered_factorizations(n)
    pairs = proper_factor_pairs(n)
    rows = pair_census(n)
    print(
        f"DEGREE {n}: primes={factors}, ordered words={len(words)}, "
        f"prime words={len(prime_words)}, loci={len(pairs)}"
    )
    print("  PRIME-REFINEMENT SHADOW")
    for word in prime_words:
        memberships = ",".join(format_pair(pair) for pair in sorted(cuts(word)))
        print(f"    {'o'.join(map(str, word))}: {memberships}")
    print("  PAIRS")
    for row in rows:
        first, second = row["loci"]
        label = f"{format_pair(first)} & {format_pair(second)}"
        if row["distance"] == 0:
            coarsest_length = min(map(len, row["common"]))
            coarsest = [word for word in row["common"] if len(word) == coarsest_length]
            refinements = ",".join("o".join(map(str, word)) for word in coarsest)
            print(f"    {label}: common refinement [{refinements}]")
        else:
            path = " -> ".join("o".join(map(str, word)) for word in row["path"])
            cores = "; ".join(row["cores"])
            print(f"    {label}: {row['distance']} core(s), {path}; {cores}")

    multistep = [row for row in rows if row["distance"] > 1]
    print(f"  PAIRS REQUIRING MORE THAN ONE CORE: {len(multistep)}")
    maximum_cuts = max(len(cuts(word)) for word in prime_words)
    print(f"  MOST LOCI ON ONE PRIME WORD: {maximum_cuts}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("degrees", nargs="*", type=int, default=DEFAULT_DEGREES)
    args = parser.parse_args()
    for degree in args.degrees:
        report(degree)


if __name__ == "__main__":
    main()

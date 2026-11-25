#!/usr/bin/env python3
"""Verify sharp split CM numerical p-ranks from endpoint bounds.

The default replays all prior minimum p-ranks, independently checks a
small endpoint census, and emits aggregate counts and a digest through
dimension 100. Use --all-minima for the complete exact rows.
"""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import comb
from pathlib import Path


PRIOR_CENSUS = Path(__file__).with_name("cm_newton_polygon_census.json")


def sharp_p_rank(height, signature):
    """Return the least numerical p-rank compatible with all inequalities."""
    mixed = tuple(m for m in signature if m > 0)
    pairs = len(mixed)
    definite_rank = height * (len(signature) - pairs)
    if not pairs:
        return definite_rank

    half = height // 2
    if max(mixed) == 1:
        if height == 2:
            return definite_rank + 2 * ((pairs + 2) // 2)
        if height % 2:
            return definite_rank + half * (pairs + 1) + 1
        return definite_rank + (half - 1) * (pairs + 1) + (pairs + 2) // 2 + 1

    if height == 4 and min(mixed) == 2:
        return definite_rank + pairs + 3

    slack = [max(0, 2 * pairs + 1 - sum(min(j, m) for m in mixed))
             for j in range(1, half + 1)]
    if height % 2 == 0:
        slack[-1] = (slack[-1] + 1) // 2
    return definite_rank + sum(slack)


def endpoint_options(height, multiplicity):
    """Generate the cost and rounded prefix profile for each endpoint choice."""
    zero_cap = height - multiplicity
    yield height, (0,) * (height - 1)
    for zeros in range(zero_cap):
        for ones in range(multiplicity):
            profile = tuple(
                int(zeros < j < height - ones) + int(ones < j < height - zeros)
                for j in range(1, height)
            )
            yield zeros + ones, profile


def brute_minimum(height, signature):
    if not any(signature):
        return height * len(signature)
    capacities = tuple(
        sum(min(j, height - j, m) for m in signature) - 1
        for j in range(1, height)
    )
    best = height * len(signature)
    for choices in product(*(endpoint_options(height, m) for m in signature)):
        rank = sum(cost for cost, _ in choices)
        if rank < best and all(
            sum(profile[j] for _, profile in choices) <= cap
            for j, cap in enumerate(capacities)
        ):
            best = rank
    return best


def first_forcing_class(height, signature):
    """A disjoint partition of the ordinary-forcing signatures."""
    mixed = tuple(m for m in signature if m > 0)
    if not mixed:
        return "all definite"
    if mixed == (1,):
        return "one minimum-one mixed pair"
    if height == 2 and mixed == (1, 1):
        return "rank two with two balanced pairs"
    if height == 4 and mixed == (2,):
        return "rank four with one balanced pair"
    return None


def rows(bound):
    for height in range(1, bound + 1):
        for pairs in range(1, bound // height + 1):
            for signature in combinations_with_replacement(range(height // 2 + 1), pairs):
                yield height, pairs, signature, sharp_p_rank(height, signature)


def verification(bound):
    previous = json.loads(PRIOR_CENSUS.read_text(encoding="utf-8"))
    expected = {(h, e, tuple(sig)): (minimum, count)
                for h, e, sig, _, count, minimum in previous["rows"]}
    minimum_digest = sha256()
    count = 0
    old_rows = 0
    brute_rows = 0
    forcing_counts = Counter()
    forcing_counts_30 = Counter()
    number_per_dimension = Counter()
    forcing_per_dimension = Counter()

    for height, pairs, signature, minimum in rows(bound):
        key = height, pairs, signature
        dimension = height * pairs
        forcing = first_forcing_class(height, signature)
        assert (minimum == dimension) == (forcing is not None)
        count += 1
        number_per_dimension[dimension] += 1
        if forcing is not None:
            forcing_counts[forcing] += 1
            forcing_per_dimension[dimension] += 1
        if dimension <= previous["max_dimension"]:
            wanted, surviving = expected[key]
            assert minimum == wanted and (minimum == dimension) == (surviving == 1), key
            old_rows += 1
            if forcing is not None:
                forcing_counts_30[forcing] += 1
        if dimension <= 20:
            assert minimum == brute_minimum(height, signature), key
            brute_rows += 1
        minimum_digest.update((json.dumps((height, pairs, signature, minimum),
                                          separators=(",", ":")) + "\n").encode())

    assert count == sum(comb(pairs + height // 2, height // 2)
                        for height in range(1, bound + 1)
                        for pairs in range(1, bound // height + 1))
    assert old_rows == len(expected) == 1188
    assert forcing_counts_30 == {
        "all definite": 111,
        "one minimum-one mixed pair": 81,
        "rank two with two balanced pairs": 14,
        "rank four with one balanced pair": 7,
    }
    return {
        "max_dimension": bound,
        "signature_count": count,
        "forcing_counts": dict(forcing_counts),
        "forcing_counts_through_30": dict(forcing_counts_30),
        "prior_minimum_rows_checked": old_rows,
        "brute_endpoint_rows_through_20": brute_rows,
        "minima_json_lines_sha256": minimum_digest.hexdigest(),
        "dimensions": [[g, number_per_dimension[g], forcing_per_dimension[g]]
                       for g in range(1, bound + 1)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=100)
    parser.add_argument("--all-minima", action="store_true")
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without -O so the verification assertions execute")
    if args.bound < 30:
        parser.error("choose a bound of at least 30 for the historical comparison")
    if args.all_minima:
        for row in rows(args.bound):
            print(json.dumps(row, separators=(",", ":")))
    else:
        print(json.dumps(verification(args.bound), indent=2))


if __name__ == "__main__":
    main()

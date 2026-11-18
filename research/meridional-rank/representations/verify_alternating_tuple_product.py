#!/usr/bin/env python3
"""Exact finite controls for the alternating tuple's ordered-product theorem.

The quantified proof is the insertion argument: right multiplication by
(1,2b)(6,2b+1) inserts 2b immediately after 1 and 2b+1 immediately after 6
in the existing long cycle. Finite replay here checks exact permutations,
the retained source tuple, its base generation words and braid stabilizers.
It does not replace that argument or independently review cover realization.
Uses only Python's standard library; default range is 4 <= b <= 64.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def cycle(n, *blocks):
    result = list(range(n))
    for block in blocks:
        for left, right in zip(block, block[1:] + block[:1]):
            result[left] = right
    return tuple(result)


def compose(left, right):
    return tuple(left[right[j]] for j in range(len(left)))


def product(permutations):
    result = tuple(range(len(permutations[0])))
    for permutation in permutations:
        result = compose(result, permutation)
    return result


def inverse(permutation):
    result = [0] * len(permutation)
    for left, right in enumerate(permutation):
        result[right] = left
    return tuple(result)


def tuple_for(b):
    n = 2 * b
    return [
        cycle(n, [0, 1], [2, 3]),
        cycle(n, [1, 4], [2, 3]),
        cycle(n, [0, 5], [1, 6]),
        cycle(n, [1, 2], [6, 7]),
    ] + [cycle(n, [1, 2*j], [6, 2*j+1]) for j in range(4, b)]


def expected_cycle(b):
    return (
        [0, 5, 1] + list(range(2*b-2, 7, -2)) + [2, 6]
        + list(range(2*b-1, 8, -2)) + [7, 4]
    )


def hurwitz(word, initial):
    current = list(initial)
    for letter in word:
        j = abs(letter) - 1
        left, right = current[j:j+2]
        if letter > 0:
            current[j:j+2] = [right, compose(compose(right, left), inverse(right))]
        else:
            current[j:j+2] = [compose(compose(left, right), inverse(left)), left]
    return current


def rejects(callback):
    try:
        callback()
    except AssertionError:
        return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-b", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.max_b >= 4, "maximum b must be at least four")
    data_path = HERE / "alternating_target_prime_families_inputs.json"
    data = json.loads(data_path.read_text())
    seeds = tuple_for(4)
    require(seeds == [tuple(p) for p in data["seed_quartet"]], "retained tuple binding")
    require(product(seeds) == cycle(8, [0, 5, 1, 2, 6, 7, 4]), "base product")

    word_a = [3, 0, 1, 3, 0, 2, 0, 1, 2]
    word_b = [0, 3, 1, 2, 3, 1, 3, 2, 1, 3]
    require(product([seeds[j] for j in word_a]) == cycle(8, list(range(7))), "A8 base word a")
    require(product([seeds[j] for j in word_b]) == cycle(8, [0, 1, 7]), "A8 base word b")

    beta = data["beta"]
    require(isinstance(beta, list) and all(isinstance(k, int) for k in beta), "literal beta")
    ramification_checks = []
    for b in range(4, args.max_b + 1):
        n = 2*b
        generators = tuple_for(b)
        for generator in generators:
            require(compose(generator, generator) == tuple(range(n)), "involution")
            require(sum(generator[j] != j for j in range(n)) == 4, "double transposition")
        ordered = product(generators)
        long_cycle = expected_cycle(b)
        require(len(long_cycle) == n-1 and set(long_cycle) == set(range(n)) - {3}, "cycle support")
        require(ordered == cycle(n, long_cycle), "ordered cycle formula")
        branch_tuple = generators + [inverse(ordered)]
        require(product(branch_tuple) == tuple(range(n)), "product-one branch tuple")
        require(2*b + (2*b-2) == 2*n - 2, "genus-zero ramification")

        braid = beta + [j for j in range(4, b) for _ in range(3)]
        require(hurwitz(braid, generators) == generators, "braid tuple stabilizer")
        if b < args.max_b:
            extended = tuple(p + (n, n+1) for p in generators)
            previous = product(extended)
            new_generator = cycle(n+2, [1, n], [6, n+1])
            inserted = compose(previous, new_generator)
            require(inserted[1] == n and inserted[n] == previous[1], "first insertion")
            require(inserted[6] == n+1 and inserted[n+1] == previous[6], "second insertion")
            require(all(inserted[j] == previous[j] for j in range(n) if j not in (1, 6)), "other arrows")
            require(inserted == cycle(n+2, expected_cycle(b+1)), "inductive cycle expression")
        if b in (4, args.max_b):
            ramification_checks.append({"b": b, "degree": n, "ramification": 2*n-2})

    altered = list(seeds)
    altered[3] = cycle(8, [0, 2], [6, 7])
    require(rejects(lambda: require(product(altered) == cycle(8, expected_cycle(4)), "mutated product")),
            "corrupt tuple rejected")
    require(rejects(lambda: require(hurwitz(beta[:-1], seeds) == seeds, "mutated braid")),
            "corrupt braid rejected")
    report = {
        "schema": "alternating-tuple-product-controls-v1",
        "status": "pass",
        "exact_arithmetic": "integer permutations; Python standard library",
        "checked_range": {"minimum_b": 4, "maximum_b": args.max_b},
        "checks": [
            "retained tuple binding", "base product", "two A8 generation words",
            "double-transposition cycle types", "ordered-product formula",
            "product-one completion", "genus-zero ramification count",
            "braid tuple stabilizer", "two insertion arrows", "corrupt tuple and braid rejection",
        ],
        "boundary_cases": ramification_checks,
        "sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [Path(__file__).resolve(), data_path, HERE / "ALTERNATING_TARGET_PRIME_FAMILIES.tex"]
        },
        "scope": "Exact bounded replay supplements the written all-b insertion proof; no cover coefficients, algebraic descent, or independent review are certified.",
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()


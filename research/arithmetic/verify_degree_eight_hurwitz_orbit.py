#!/usr/bin/env python3
"""Complete finite Hurwitz-orbit certificate in degree eight.

Run: python3 research/arithmetic/verify_degree_eight_hurwitz_orbit.py
Optional --write-receipt PATH. Uses only the Python standard library.
All 210 double transpositions in S8 are included. Meet-in-the-middle
enumeration includes every quadruple with the fixed literal knot product.
No sampling, floating-point continuation, or group-library oracle is used.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "research/meridional-rank/representations/alternating_target_prime_families_inputs.json"
IDENTITY = tuple(range(8))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def compose(left, right):
    return tuple(left[right[j]] for j in range(8))


def inverse(value):
    result = [0] * 8
    for j, image in enumerate(value):
        result[image] = j
    return tuple(result)


def permutation(*cycles):
    result = list(IDENTITY)
    for cycle in cycles:
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            result[a] = b
    return tuple(result)


def ordered_product(values):
    result = IDENTITY
    for value in values:
        result = compose(result, value)
    return result


def generated_order(values):
    found = {IDENTITY}
    queue = [IDENTITY]
    for value in queue:
        for generator in values:
            candidate = compose(value, generator)
            if candidate not in found:
                found.add(candidate)
                queue.append(candidate)
    return len(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    conjugacy_class = []
    for a, b, c, d in combinations(range(8), 4):
        for pairing in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
            conjugacy_class.append(permutation(*pairing))
    conjugacy_class.sort()
    index = {value: j for j, value in enumerate(conjugacy_class)}
    require(len(index) == len(conjugacy_class) == 210, "complete DT class")
    for value in conjugacy_class:
        require(compose(value, value) == IDENTITY, "involutions")
        require(sum(j != value[j] for j in range(8)) == 4, "DT support")

    inputs = json.loads(INPUT.read_text())
    seed_values = tuple(tuple(value) for value in inputs["seed_quartet"])
    seed = tuple(index[value] for value in seed_values)
    target = permutation((0, 5, 1, 2, 6, 7, 4))
    require(ordered_product(seed_values) == target, "literal knot product")

    pairs = defaultdict(list)
    for a, left in enumerate(conjugacy_class):
        for b, right in enumerate(conjugacy_class):
            pairs[compose(left, right)].append((a, b))
    # abcd = target iff cd = (ab)^(-1) target. Every ordered pair is
    # included exactly once in pairs; this is an exhaustive finite reduction.
    all_tuples = set()
    for ab, left_pairs in pairs.items():
        right_pairs = pairs.get(compose(inverse(ab), target), ())
        for a, b in left_pairs:
            for c, d in right_pairs:
                all_tuples.add((a, b, c, d))
    require(len(all_tuples) == 93296, "complete quadruple count")
    require(seed in all_tuples, "seed retained")

    # Independent direct completion on all tuples with this first factor.
    first = seed[0]
    direct = set()
    for b, second in enumerate(conjugacy_class):
        ab = compose(conjugacy_class[first], second)
        for c, third in enumerate(conjugacy_class):
            fourth = compose(inverse(compose(ab, third)), target)
            if fourth in index:
                direct.add((first, b, c, index[fourth]))
    require(direct == {value for value in all_tuples if value[0] == first},
            "independent direct-completion control")

    conjugation = [
        [index[compose(compose(left, right), left)] for right in conjugacy_class]
        for left in conjugacy_class
    ]

    def hurwitz(value, j, sign):
        a, b = value[j:j+2]
        pair = (b, conjugation[b][a]) if sign == 1 else (conjugation[a][b], a)
        return value[:j] + pair + value[j+2:]

    for j in range(3):
        require(hurwitz(hurwitz(seed, j, 1), j, -1) == seed,
                "Hurwitz inverse control")
        moved = hurwitz(seed, j, 1)
        require(ordered_product([conjugacy_class[i] for i in moved]) == target,
                "Hurwitz product control")

    remaining = all_tuples.copy()
    components = []
    while remaining:
        initial = seed if not components else min(remaining)
        remaining.remove(initial)
        queue = [initial]
        for value in queue:
            for j in range(3):
                for sign in (-1, 1):
                    candidate = hurwitz(value, j, sign)
                    require(candidate in all_tuples, "Hurwitz closure")
                    if candidate in remaining:
                        remaining.remove(candidate)
                        queue.append(candidate)
        values = tuple(conjugacy_class[j] for j in initial)
        order = generated_order(values)
        digest = hashlib.sha256(b"".join(bytes(v) for v in sorted(queue))).hexdigest()
        components.append({
            "size": len(queue),
            "generated_group_order": order,
            "representative_indices": list(initial),
            "representative_permutations": [list(v) for v in values],
            "sorted_tuple_sha256": digest,
        })
    require(len(components) == 8, "component count")
    require(sum(v["size"] for v in components) == len(all_tuples), "partition")
    full = [value for value in components if value["generated_group_order"] == 20160]
    require(len(full) == 1 and full[0]["size"] == 47334, "unique A8 component")
    require(components[0] == full[0], "literal knot lies in A8 component")
    require(sorted((v["size"], v["generated_group_order"]) for v in components) ==
            [(1029, 168), (1029, 168), (1029, 1344), (1029, 1344),
             (10976, 2520), (15435, 2520), (15435, 2520), (47334, 20160)],
            "complete component census")

    # Deliberately invalid controls: wrong product, unconjugated swap,
    # and a proper generated subgroup cannot certify full A8.
    wrong = (seed[1],) + seed[1:]
    require(wrong not in all_tuples, "wrong-product rejection")
    swapped = (seed[1], seed[0]) + seed[2:]
    require(swapped not in all_tuples, "unconjugated-swap rejection")
    require(any(v["generated_group_order"] == 1344 for v in components),
            "proper transitive subgroup distinguished")

    report = {
        "schema_version": 1,
        "kind": "complete_finite_hurwitz_certificate",
        "python": platform.python_version(),
        "arithmetic": "exact permutations of eight symbols; standard library",
        "command": "python3 research/arithmetic/verify_degree_eight_hurwitz_orbit.py",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "class_size": len(conjugacy_class),
        "pair_product_count": len(pairs),
        "quadruple_count": len(all_tuples),
        "target_product": list(target),
        "components": components,
        "unique_A8_orbit": True,
        "A8_orbit_size": 47334,
        "A8_absolute_Nielsen_class_size": 6762,
        "controls": ["direct completion for fixed first factor", "Hurwitz inverses",
                     "product preservation", "wrong-product rejection",
                     "unconjugated-swap rejection", "proper subgroup distinction"],
        "scope": (
            "Complete B4 orbit computation on DT quadruples with fixed 7-cycle product. "
            "Together with the written reduction it proves the unique absolute A8 "
            "Nielsen orbit and identifies every characteristic-zero A8 cover with "
            "this passport with the knot tuple up to braid and sheet relabeling. "
            "It does not prove equality of individual covers or an all-degree orbit theorem."
        ),
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
        "result": "PASS",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.write_receipt:
        args.write_receipt.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()

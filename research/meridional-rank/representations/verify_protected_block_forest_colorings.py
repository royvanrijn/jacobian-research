#!/usr/bin/env python3
"""Complete exact tree-component census for the protected U block.

The main census fixes the product cycle and its cyclic first-color orbits.
A separate reference census through degree six fixes the first color under
the full symmetric group and tests graph connectivity directly instead.
The note proves the finite reduction and its global cycle-space consequence.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("This exact census requires Python assertions enabled.")

import argparse
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from math import comb
from pathlib import Path
import platform


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
U = (1, 2, -3, -3, 2, -3, -3, 2, 1, 1, 2, 2, 2, -3)
EXPECTED_MAIN = {
    2: (2, 4, 0), 3: (2, 12, 3), 4: (5, 49, 0),
    5: (6, 163, 0), 6: (13, 647, 0), 7: (19, 1736, 0),
    8: (33, 4395, 0), 9: (42, 6561, 0),
}


def require(value, message):
    if not value:
        raise AssertionError(message)


def multiply(a, b):
    return tuple(a[j] for j in b)


def transposition_length(p):
    seen = set()
    cycles = 0
    for i in range(len(p)):
        if i in seen:
            continue
        cycles += 1
        while i not in seen:
            seen.add(i)
            i = p[i]
    return len(p) - cycles


def involutions(n):
    answer = [(tuple(range(n)), 0)]
    for a, b in combinations(range(n), 2):
        p = list(range(n))
        p[a], p[b] = b, a
        answer.append((tuple(p), 1))
    for a, b, c, d in combinations(range(n), 4):
        for i, j, k, ell in ((a, b, c, d), (a, c, b, d), (a, d, b, c)):
            p = list(range(n))
            p[i], p[j] = j, i
            p[k], p[ell] = ell, k
            answer.append((tuple(p), 2))
    require(len(answer) == 1 + comb(n, 2) + 3*comb(n, 4),
            "complete zero/single/double-transposition domain")
    require(len(dict(answer)) == len(answer), "no duplicate involutions")
    for p, edges in answer:
        require(multiply(p, p) == tuple(range(n)), "actual involution")
        require(transposition_length(p) == edges, "exact edge count")
    return answer


def rotate(p, shift):
    n = len(p)
    return tuple((p[(j - shift) % n] + shift) % n for j in range(n))


def braid_action(values, braid=U):
    values = list(values)
    for letter in braid:
        j = abs(letter) - 1
        a, b = values[j:j + 2]
        # All propagated entries remain involutions, so inverse equals itself.
        values[j:j + 2] = (
            [b, multiply(multiply(b, a), b)] if letter > 0
            else [multiply(multiply(a, b), a), a]
        )
    return tuple(values)


def protected(values):
    output = braid_action(values)
    return output[0] == values[0] and output[3] == values[3]


def graph_cycle_rank(values):
    """Occurrence edges remain distinct even when their endpoint pairs repeat."""
    n = len(values[0])
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    edges = 0
    for p in values:
        for a, b in enumerate(p):
            if a < b:
                edges += 1
                parent[find(a)] = find(b)
    components = len({find(i) for i in range(n)})
    return edges - n + components, components


def cycle_normalized_census():
    records = []
    witnesses = {}
    for n in range(2, 10):
        inv = involutions(n)
        lookup = dict(inv)
        cycle = tuple((j + 1) % n for j in range(n))

        @lru_cache(None)
        def successors(residual):
            length = transposition_length(residual)
            answer = []
            for g, edges in inv:
                if edges > length:
                    continue
                new = multiply(g, residual)
                if transposition_length(new) == length - edges:
                    answer.append((g, edges, new))
            return tuple(answer)

        firsts = factorizations = returning = 0
        for a, a_edges in inv:
            if a != min(rotate(a, shift) for shift in range(n)):
                continue
            if a_edges + 6 < n - 1:
                continue
            firsts += 1
            residual = multiply(a, cycle)
            if transposition_length(residual) != n - 1 - a_edges:
                continue
            for b, b_edges, residual2 in successors(residual):
                if n - 1 - a_edges - b_edges > 4:
                    continue
                for c, c_edges, d in successors(residual2):
                    if d not in lookup:
                        continue
                    d_edges = lookup[d]
                    vector = (a_edges, b_edges, c_edges, d_edges)
                    if sum(vector) != n - 1:
                        continue
                    values = (a, b, c, d)
                    require(multiply(multiply(multiply(a, b), c), d) == cycle,
                            "forced fourth factor and product convention")
                    require(graph_cycle_rank(values) == (0, 1),
                            "every accepted minimal product is an actual tree")
                    factorizations += 1
                    if protected(values):
                        returning += 1
                        witnesses.setdefault(vector, values)
        record = (firsts, factorizations, returning)
        require(record == EXPECTED_MAIN[n], "complete cycle-normalized census")
        records.append({"vertices": n, "first_color_representatives": firsts,
                        "minimal_factorizations": factorizations,
                        "outer_returning_factorizations": returning,
                        "cached_residuals": successors.cache_info().currsize})
    require(set(witnesses) == {(0, 0, 1, 1)}, "only permitted tree-component vector")
    return records, witnesses


def first_color_normalized_reference():
    """Independent finite reduction: no cycle normalization or length pruning."""
    records = []
    all_vectors = set()
    for n in range(2, 7):
        inv = involutions(n)
        by_edges = {i: [p for p, edges in inv if edges == i] for i in range(3)}
        identity = tuple(range(n))
        first = [(identity, 0)]
        p = list(identity)
        p[0], p[1] = 1, 0
        first.append((tuple(p), 1))
        if n >= 4:
            p = list(identity)
            p[0], p[1], p[2], p[3] = 1, 0, 3, 2
            first.append((tuple(p), 2))
        trees = returning = 0
        for a, a_edges in first:
            for b, b_edges in inv:
                for c, c_edges in inv:
                    remaining = n - 1 - a_edges - b_edges - c_edges
                    if remaining not in by_edges:
                        continue
                    for d in by_edges[remaining]:
                        values = (a, b, c, d)
                        if graph_cycle_rank(values) != (0, 1):
                            continue
                        trees += 1
                        if protected(values):
                            returning += 1
                            all_vectors.add((a_edges, b_edges, c_edges, remaining))
        require((returning > 0) == (n == 3),
                "reference normalization has exactly the same positive degree")
        records.append({"vertices": n, "tree_factorizations": trees,
                        "outer_returning_factorizations": returning})
    require(all_vectors == {(0, 0, 1, 1)}, "independent first-color census vector")
    return records


def permutation(n, *edges):
    p = list(range(n))
    for a, b in edges:
        p[a], p[b] = b, a
    return tuple(p)


def controls_and_packing(witnesses):
    old = (permutation(8, (0, 1), (2, 3)),
           permutation(8, (1, 4), (2, 3)),
           permutation(8, (0, 5), (1, 6)),
           permutation(8, (1, 2), (6, 7)))
    require(protected(old), "old alternating eight-symbol seed still permitted")
    require(graph_cycle_rank(old) == (1, 1), "old seed has exactly one occurrence cycle")
    for values in (old, witnesses[(0, 0, 1, 1)]):
        inverse_braid = tuple(-s for s in U[::-1])
        require(braid_action(braid_action(values), inverse_braid) == values,
                "signed Hurwitz moves invert exactly")
    identity = tuple(range(2))
    first_only = (identity, identity, (1, 0), identity)
    image = braid_action(first_only)
    require(image[0] == first_only[0] and image[3] != first_only[3],
            "one protected outer relation alone is insufficient")
    path9 = (permutation(9, (0, 1), (2, 3)),
             permutation(9, (1, 2), (4, 5)),
             permutation(9, (3, 4), (6, 7)),
             permutation(9, (5, 6), (7, 8)))
    require(graph_cycle_rank(path9) == (0, 1), "four-double-transposition path tree")
    require(not protected(path9), "explicit nine-symbol tree violates outer closure")
    # Complete finite disjoint-component packing. Only occurrence counts matter.
    reachable = {(0, 0, 0, 0)}
    for _ in range(8):
        reachable |= {tuple(x + y for x, y in zip(v, w))
                      for v in tuple(reachable) for w in witnesses
                      if all(x + y <= 2 for x, y in zip(v, w))}
    require((2, 2, 2, 2) not in reachable, "no forest made from admitted components")
    require(len(reachable) == 3, "complete packing state count")
    return {"old_seed_occurrence_cycle_rank": 1,
            "old_seed_outer_return": True,
            "single_outer_relation_negative_control": True,
            "explicit_nine_symbol_tree_rejected": True,
            "packing_states": [list(v) for v in sorted(reachable)],
            "four_double_transposition_forest_exists": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    records, witnesses = cycle_normalized_census()
    reference = first_color_normalized_reference()
    controls = controls_and_packing(witnesses)
    pins = [Path(__file__).resolve(), HERE / "MULTIPLE_PROTECTED_OBSTRUCTIONS.md"]
    result = {
        "status": "passed", "python": platform.python_version(),
        "braid": list(U), "arithmetic": "exact finite permutations and integer graph ranks",
        "source_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in pins},
        "complete_cycle_normalized_census": records,
        "independent_first_color_normalized_reference_through": 6,
        "reference_census": reference,
        "permitted_component_vectors": [list(v) for v in sorted(witnesses)],
        "component_witnesses": [{"vector": list(v), "tuple": [list(p) for p in t]}
                                for v, t in sorted(witnesses.items())],
        "controls": controls,
        "scope": {
            "complete_domain": "all connected tree components of a four-involution tuple, each color a product of at most two disjoint transpositions; total at most eight occurrence edges",
            "two_outer_relations_only": True,
            "written_consequence": "Every local four-double-transposition U coloring has a nonzero occurrence-cycle; k disjoint blocks force transitive support v<=2b-k+1",
            "no_transitive_nine_symbol_coloring": True,
            "all_local_support_at_most_eight": "not asserted for disconnected actions",
            "prime_gap_two_realization": False,
            "formal_verification": False,
            "external_review": False,
            "historical_priority": False,
        },
    }
    if args.write_receipt:
        require(args.write_receipt.suffix == ".json", "receipt must have .json suffix")
        args.write_receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "census": records, "reference": reference,
                      "permitted_vectors": result["permitted_component_vectors"],
                      "forest_tuple_exists": False}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact controls for the multiple-protected-obstruction proof.

The note supplies unrestricted proofs. This checker does not replay the
imported local signed obstruction or prove prime/gap-two realization.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
U = [1, 2, -3, -3, 2, -3, -3, 2, 1, 1, 2, 2, 2, -3]


def require(value, message):
    if not value:
        raise AssertionError(message)


def reduced(word):
    answer = []
    for letter in word:
        if answer and answer[-1] == -letter:
            answer.pop()
        else:
            answer.append(letter)
    return tuple(answer)


def inverse(word):
    return tuple(-letter for letter in word[::-1])


def multiply(*words):
    return reduced(letter for word in words for letter in word)


def hurwitz(values, braid):
    values = list(values)
    for letter in braid:
        j = abs(letter) - 1
        a, c = values[j:j + 2]
        values[j:j + 2] = (
            [c, multiply(inverse(c), a, c)] if letter > 0
            else [multiply(a, c, inverse(a)), a]
        )
    return values


def band(i, j):
    require(1 <= i < j, "ordered band endpoints")
    return list(range(i, j)) + list(range(-(j - 2), -i + 1))


def prefix(k):
    return [letter + 4*i if letter > 0 else letter - 4*i
            for i in range(k) for letter in U]


def mobile_positions(k, b):
    protected = {4*i + a for i in range(k) for a in (1, 4)}
    return [i for i in range(1, b + 1) if i not in protected]


def tail(k, b):
    mobile = mobile_positions(k, b)
    return [letter for i, j in zip(mobile, mobile[1:])
            for letter in band(i, j)]


def strand_permutation(b, braid):
    answer = list(range(b))
    for letter in braid:
        j = abs(letter) - 1
        answer[j], answer[j + 1] = answer[j + 1], answer[j]
    return answer


def cycles(p):
    seen = set()
    answer = []
    for i in range(len(p)):
        if i in seen:
            continue
        orbit = []
        while i not in seen:
            seen.add(i)
            orbit.append(i)
            i = p[i]
        answer.append(orbit)
    return answer


def rank(rows):
    if not rows:
        return 0
    matrix = [[Fraction(x) for x in row] for row in rows]
    pivot = 0
    for j in range(len(matrix[0])):
        selected = next((i for i in range(pivot, len(matrix))
                         if matrix[i][j]), None)
        if selected is None:
            continue
        matrix[pivot], matrix[selected] = matrix[selected], matrix[pivot]
        scale = matrix[pivot][j]
        matrix[pivot] = [x / scale for x in matrix[pivot]]
        for i in range(len(matrix)):
            if i != pivot and matrix[i][j]:
                scale = matrix[i][j]
                matrix[i] = [x - scale*y for x, y in
                             zip(matrix[i], matrix[pivot])]
        pivot += 1
        if pivot == len(matrix):
            break
    return pivot


def dot(a, c):
    return sum(x*y for x, y in zip(a, c))


def band_controls():
    count = 0
    for b in range(2, 13):
        values = [(i,) for i in range(1, b + 1)]
        for i in range(1, b):
            for j in range(i + 1, b + 1):
                W = tuple(range(i + 1, j))
                expected = values.copy()
                expected[i - 1] = multiply(W, (j,), inverse(W))
                expected[j - 1] = multiply((-j,), inverse(W), (i,), W, (j,))
                require(hurwitz(values, band(i, j)) == expected,
                        "exact general-band free-group formula")
                count += 1
    # A swapped conjugation convention does not fix the protected middle word.
    require(hurwitz([(1,), (2,), (3,)], [-1, 2, 1])[1] != (2,),
            "reject opposite band conjugation")
    return count


def knot_controls():
    require(strand_permutation(4, U) == [1, 0, 3, 2],
            "literal U matching permutation")
    count = 0
    explicit = None
    for k in range(1, 9):
        for extra in range(4):
            b = 4*k + extra
            P, R = prefix(k), tail(k, b)
            mobile = mobile_positions(k, b)
            s = len(mobile)
            protected = set(range(1, b + 1)) - set(mobile)
            values = [(i,) for i in range(1, b + 1)]
            after = hurwitz(values, R)
            require(all(after[j - 1] == (j,) for j in protected),
                    "all protected meridian words fixed")
            require(len([c for c in cycles(strand_permutation(b, R))
                         if len(c) > 1]) == 1,
                    "one mobile permutation cycle")
            for ell in range(4):
                braid = P + R*(1 + s*ell)
                require(len(cycles(strand_permutation(b, braid))) == 1,
                        "literal full braid has one-component closure")
                count += 1
            if k == 4 and extra == 0:
                explicit = {"k": k, "strands": b, "ell": 0,
                            "braid": P + R,
                            "augmentation_upper_bound": b - (k + 1)//2,
                            "meridional_lower_bound": "not supplied",
                            "primality": "not supplied"}
    return count, explicit


def defect_controls():
    # Exact checks of every assignment of k local alternative sides, k<=12.
    count = 0
    for k in range(1, 13):
        for mask in range(1 << k):
            root = mask.bit_count()
            covector = k - root
            require(root + covector == k, "one chosen side per block")
            require(max(root, covector) >= (k + 1)//2,
                    "irreducible common deficit bound")
            count += 1
    # Two blocks need not yield TWO root defects and TWO covector defects.
    roots = [[int(i == j) for j in range(7)] for i in range(3)]
    roots += [[1, 1, 1, 0, 0, 0, 0]]
    roots += [[int(i == j) for j in range(7)] for i in range(3, 7)]
    covectors = [
        [2, 1, -1, 2, 1, -1, 0],
        [-2, 2, -1, 2, 2, -1, 3],
        [1, 2, 2, -2, 1, 0, -2],
        [1, 0, 1, -2, 3, 1, 0],
        [-2, -1, 2, 2, -2, 2, -2],
        [1, 3, 2, -2, 2, 3, -1],
        [3, 2, 1, 2, 3, 2, 5],
        [2, 4, 5, 2, 3, 7, 2],
    ]
    require([dot(u, f) for u, f in zip(roots, covectors)] == [2]*8,
            "all rank-one matrices have determinant mu=-1")
    require([rank(roots[:4]), rank(covectors[:4]),
             rank(roots[4:]), rank(covectors[4:])] == [3, 4, 4, 3],
            "opposite one-sided local defects")
    require(rank(roots) == rank(covectors) == 7,
            "both full spans have dimension b-1")
    # If a nonzero invariant W contains one root, nonzero evaluations
    # propagate it along this graph; full covector rank supplies that
    # first root. Strong connectivity therefore proves irreducibility.
    for start in range(8):
        reached = {start}
        while True:
            added = reached | {j for i in reached for j in range(8)
                               if dot(covectors[j], roots[i])}
            if added == reached:
                break
            reached = added
        require(len(reached) == 8, "strongly connected pairing graph")
    return count, {"roots": roots, "covectors": covectors,
                   "global_r_q": [7, 7],
                   "block_r_q": [[3, 4], [4, 3]],
                   "mu": -1,
                   "irreducibility": "full spans and strongly connected pairing graph",
                   "scope": "linear dependence alternatives only; no U relations claimed"}


def edge_controls():
    seed = [((0, 1), (2, 3)), ((1, 4), (2, 3)),
            ((0, 5), (1, 6)), ((1, 2), (6, 7))]
    edges = {tuple(sorted(edge)) for matching in seed for edge in matching}
    require(len(edges) == 7, "seven distinct edges in eight occurrences")
    reached = {0}
    while True:
        added = reached | {a for a, b in edges if b in reached}
        added |= {b for a, b in edges if a in reached}
        if added == reached:
            break
        reached = added
    require(len(reached) == 8, "seed support graph connected")
    require(7*4 < 2*16 - 1, "four literal seeds cannot color transitive A32")
    for k in range(1, 65):
        for extra in range(5):
            b = 4*k + extra
            require(7*k + 2*extra == 2*b - k, "literal edge budget")
            optimistic_lower = (2*b - k + 1)//2
            augmentation_upper = b - (k + 1)//2
            require(optimistic_lower - augmentation_upper == k % 2,
                    "these bounds certify at most gap one")
    return {"seed": seed, "distinct_edges": sorted(edges),
            "per_block_occurrences": 8, "per_block_distinct_edges": 7,
            "universal_budget": "v <= 2*b-k+1",
            "limitation": "support-graph lower bound and stacked upper bound only"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    bands = band_controls()
    knots, explicit = knot_controls()
    assignments, sharpness = defect_controls()
    edges = edge_controls()
    pins = [HERE / "MULTIPLE_PROTECTED_OBSTRUCTIONS.md", Path(__file__),
            ROOT / "results/prime-knot-augmentation-gaps/representations.tex"]
    result = {
        "status": "passed",
        "python": platform.python_version(),
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in pins},
        "exact_free_group_band_controls": bands,
        "bounded_knot_permutation_controls": knots,
        "bounded_defect_assignment_controls": assignments,
        "explicit_sixteen_strand_knot": explicit,
        "two_block_linear_algebra_control": sharpness,
        "literal_seed_coloring_obstruction": edges,
        "scope": {
            "written_unrestricted_results": [
                "r+q <= 2*b-k for disjoint protected blocks",
                "irreducible dimension <= b-ceil(k/2)",
                "explicit knot closure and protection formula",
                "literal old-seed transitive support v <= 2*b-k+1"],
            "local_signed_obstruction_replayed": False,
            "matching_meridional_lower_bound": False,
            "prime_or_hyperbolic_realization": False,
            "gap_two_knot": False,
            "formal_verification": False,
            "external_review": False,
        },
    }
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "band_controls": bands,
                      "knot_controls": knots, "defect_assignments": assignments,
                      "explicit_ar_upper": explicit["augmentation_upper_bound"],
                      "gap_two_proved": False}, indent=2))


if __name__ == "__main__":
    main()

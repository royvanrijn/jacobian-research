#!/usr/bin/env python3
"""Exact signed-transposition lifts of a specified symmetric quotient.

Input is a finite presentation, with *every* generator a meridian, and its
transposition image (zero-based endpoints). Relator letters are signed,
one-based generator indices. This module verifies the quotient relations,
not that a user-supplied presentation is actually a given knot group.
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def row_basis(rows: list[int]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for value in rows:
        while value:
            p = value.bit_length() - 1
            if p in basis:
                value ^= basis[p]
            else:
                basis[p] = value
                break
    return basis


def contains(basis: dict[int, int], value: int) -> bool:
    while value:
        p = value.bit_length() - 1
        if p not in basis:
            return False
        value ^= basis[p]
    return True


def nullspace(rows: list[int], n: int) -> list[int]:
    basis = row_basis(rows)
    result = []
    for free in range(n):
        if free in basis:
            continue
        value = 1 << free
        for p in sorted(basis):
            if (basis[p] & value).bit_count() % 2:
                value ^= 1 << p
        result.append(value)
    assert all((row & v).bit_count() % 2 == 0 for row in rows for v in result)
    return result


def signed_swap(x: int, p: int, q: int, a: int) -> int:
    if ((x >> p) ^ (x >> q) ^ a) & 1:
        x ^= (1 << p) | (1 << q)
    return x


def sign_potentials(n: int, edges: list[tuple[int, int]], signs: int):
    """Solve b_p+b_q=a_e in ambient F2^n, including parallel edges."""
    adjacency = [[] for _ in range(n)]
    for i, (p, q) in enumerate(edges):
        a = (signs >> i) & 1
        adjacency[p].append((q, a))
        adjacency[q].append((p, a))
    values = [None] * n
    components = 0
    consistent = True
    for root in range(n):
        if values[root] is not None:
            continue
        components += 1
        values[root] = 0
        stack = [root]
        while stack:
            p = stack.pop()
            for q, a in adjacency[p]:
                expected = values[p] ^ a
                if values[q] is None:
                    values[q] = expected
                    stack.append(q)
                elif values[q] != expected:
                    consistent = False
    return components, values if consistent else None


def orbit(n: int, edges: list[tuple[int, int]], signs: int, start: int = 0) -> set[int]:
    seen = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for i, (p, q) in enumerate(edges):
            y = signed_swap(x, p, q, (signs >> i) & 1)
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return seen


def relator_rows(n: int, edges: list[tuple[int, int]], relators: list[list[int]]) -> list[int]:
    """Evaluate each affine relator at zero with symbolic binary signs."""
    rows = []
    for ri, word in enumerate(relators):
        perm = list(range(n))
        shifts = [0] * n  # each entry is a linear form in the signs
        for letter in word:
            if not isinstance(letter, int) or letter == 0 or abs(letter) > len(edges):
                raise ValueError(f"invalid generator letter in relator {ri}")
            i = abs(letter) - 1
            p, q = edges[i]
            # Signed swaps are involutions: positive and negative letters agree.
            shifts[perm[p]] ^= 1 << i
            shifts[perm[q]] ^= 1 << i
            perm[p], perm[q] = perm[q], perm[p]
        if perm != list(range(n)):
            raise ValueError(f"the proposed symmetric quotient fails relator {ri}")
        rows.extend(row for row in shifts if row)
    return rows


def solve(data: dict) -> dict:
    n = data["degree"]
    if not isinstance(n, int) or n < 3:
        raise ValueError("symmetric degree must be at least three")
    edges = [tuple(e) for e in data["transpositions"]]
    if any(len(e) != 2 or not all(isinstance(p, int) and 0 <= p < n for p in e)
           or e[0] == e[1] for e in edges):
        raise ValueError("each meridian image must be a transposition")
    if sign_potentials(n, edges, 0)[0] != 1:
        raise ValueError("transposition graph is disconnected: image is not S_n")
    rows = relator_rows(n, edges, data["relators"])
    lifts = nullspace(rows, len(edges))
    cuts = [sum((int(p == v) ^ int(q == v)) << i
                for i, (p, q) in enumerate(edges)) for v in range(n)]
    cuts_basis = row_basis(cuts)
    assert len(cuts_basis) == n - 1
    assert all((row & c).bit_count() % 2 == 0 for row in rows for c in cuts)
    witness = next((a for a in lifts if not contains(cuts_basis, a)), None)
    result = {
        "label": data.get("label"),
        "degree": n,
        "meridional_generators": len(edges),
        "relators_checked": len(data["relators"]),
        "symmetric_quotient_verified": True,
        "lift_space_dimension": len(lifts),
        "meridional_rank_lower_bound": len(lifts),
        "bound_mechanism": "counting all lifts of the fixed symmetric quotient; two possible images per meridian",
        "ambient_cut_space_dimension": n - 1,
        "unbalanced_lift_dimension": len(lifts) - (n - 1),
        "transitive_signed_lift_exists": witness is not None,
        "conditional_on_presentation_and_meridian_provenance": True,
    }
    if witness is not None:
        assert sign_potentials(n, edges, witness)[1] is None
        result.update({
            "witness_signs": [(witness >> i) & 1 for i in range(len(edges))],
            "sheets": 1 << (n - 1),
            "single_transitive_signed_quotient_lower_bound": n,
            "image": f"W(D_{n}) = F2^{n-1} semidirect S_{n}",
        })
        if n <= 14:
            result["independent_orbit_size"] = len(orbit(n, edges, witness))
            assert result["independent_orbit_size"] == 1 << (n - 1)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("presentation", type=Path)
    args = parser.parse_args()
    print(json.dumps(solve(json.loads(args.presentation.read_text())), indent=2))


if __name__ == "__main__":
    main()

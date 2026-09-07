#!/usr/bin/env python3
"""Finite support census for degree-four diagonal G_m-equivariant maps of A^3.

The script enumerates the hyperplane arrangement

    (alpha-e_i).w = 0,  2 <= |alpha| <= 4,

up to coordinate permutation and simultaneous reversal of all weights.  It
records both rank-two intersections (isolated primitive weight rays) and
rank-one strata (generic points of a single relation plane).  A directed
edge j -> i means that an allowed nonlinear monomial in the i-th component
uses x_j.  Self edges are omitted; acyclicity of this cross-dependency graph
is the combinatorial triangularity criterion used in the accompanying note.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


Exponent = tuple[int, int, int]
Weight = tuple[int, int, int]
Relation = tuple[int, int, int]


@dataclass(frozen=True)
class Monomial:
    component: int
    exponent: Exponent
    relation: Relation


def primitive(vector: Iterable[int], identify_sign: bool = False) -> Weight:
    result = tuple(vector)
    divisor = math.gcd(*(abs(entry) for entry in result))
    if divisor:
        result = tuple(entry // divisor for entry in result)
    if identify_sign:
        first = next((entry for entry in result if entry), 0)
        if first < 0:
            result = tuple(-entry for entry in result)
    return result  # type: ignore[return-value]


def nonlinear_monomials(degree_bound: int = 4) -> tuple[Monomial, ...]:
    result: list[Monomial] = []
    for component in range(3):
        for degree in range(2, degree_bound + 1):
            for alpha_0 in range(degree + 1):
                for alpha_1 in range(degree - alpha_0 + 1):
                    exponent = (
                        alpha_0,
                        alpha_1,
                        degree - alpha_0 - alpha_1,
                    )
                    relation = list(exponent)
                    relation[component] -= 1
                    result.append(
                        Monomial(component, exponent, tuple(relation))
                    )
    return tuple(result)


MONOMIALS = nonlinear_monomials()
RELATIONS = tuple(sorted({monomial.relation for monomial in MONOMIALS}))


def dot(left: Iterable[int], right: Iterable[int]) -> int:
    return sum(a * b for a, b in zip(left, right))


def cross(left: Relation, right: Relation) -> Weight:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def support(weight: Weight) -> tuple[Monomial, ...]:
    return tuple(
        monomial for monomial in MONOMIALS
        if dot(monomial.relation, weight) == 0
    )


def relation_rank(entries: Iterable[Relation]) -> int:
    entries = tuple(entries)
    if not entries:
        return 0
    pivot = entries[0]
    if all(not any(cross(pivot, entry)) for entry in entries[1:]):
        return 1
    return 2


def dependency_edges(entries: Iterable[Monomial]) -> frozenset[tuple[int, int]]:
    return frozenset(
        (variable, monomial.component)
        for monomial in entries
        for variable, exponent in enumerate(monomial.exponent)
        if variable != monomial.component and exponent
    )


def has_directed_cycle(edges: frozenset[tuple[int, int]]) -> bool:
    for first in range(3):
        for second in range(first + 1, 3):
            if (first, second) in edges and (second, first) in edges:
                return True
    return (
        {(0, 1), (1, 2), (2, 0)} <= edges
        or {(0, 2), (2, 1), (1, 0)} <= edges
    )


def permute_weight(weight: Weight, permutation: tuple[int, int, int]) -> Weight:
    return tuple(weight[index] for index in permutation)  # type: ignore[return-value]


def canonical_weight(weight: Weight) -> Weight:
    return min(
        primitive(permute_weight(weight, permutation), identify_sign=True)
        for permutation in itertools.permutations(range(3))
    )


def canonical_support_key(entries: Iterable[Monomial]) -> tuple[tuple[int, Exponent], ...]:
    raw = tuple((entry.component, entry.exponent) for entry in entries)
    candidates = []
    for permutation in itertools.permutations(range(3)):
        inverse = {old: new for new, old in enumerate(permutation)}
        transformed = []
        for component, exponent in raw:
            transformed.append(
                (
                    inverse[component],
                    tuple(exponent[old] for old in permutation),
                )
            )
        candidates.append(tuple(sorted(transformed)))
    return min(candidates)


def isolated_weight_rays() -> set[Weight]:
    result = set()
    for left, right in itertools.combinations(RELATIONS, 2):
        candidate = cross(left, right)
        if any(candidate):
            result.add(primitive(candidate, identify_sign=True))
    return result


def rank_one_supports() -> dict[tuple[tuple[int, Exponent], ...], Relation]:
    """Return support types on generic points of each relation plane."""
    result = {}
    for relation in RELATIONS:
        line = primitive(relation, identify_sign=True)
        entries = tuple(
            monomial for monomial in MONOMIALS
            if not any(cross(line, monomial.relation))
        )
        key = canonical_support_key(entries)
        result.setdefault(key, line)
    return result


def monomial_label(entry: Monomial) -> str:
    names = "xyz"
    factors = []
    for name, exponent in zip(names, entry.exponent):
        if exponent == 1:
            factors.append(name)
        elif exponent:
            factors.append(f"{name}^{exponent}")
    return f"F{entry.component + 1}:{''.join(factors)}"


def census() -> dict:
    isolated_by_support: dict[tuple[tuple[int, Exponent], ...], set[Weight]] = (
        defaultdict(set)
    )
    for weight in isolated_weight_rays():
        entries = support(weight)
        if relation_rank(entry.relation for entry in entries) == 2:
            isolated_by_support[canonical_support_key(entries)].add(
                canonical_weight(weight)
            )

    records = []
    for key, weights in isolated_by_support.items():
        representative = min(weights)
        entries = support(representative)
        edges = dependency_edges(entries)
        records.append(
            {
                "kind": "rank2",
                "weight": representative,
                "support_size": len(entries),
                "cyclic": has_directed_cycle(edges),
                "edges": sorted(edges),
                "support": [monomial_label(entry) for entry in entries],
            }
        )

    for key, relation in rank_one_supports().items():
        entries = tuple(
            Monomial(component, exponent, tuple(
                exponent[index] - (1 if index == component else 0)
                for index in range(3)
            ))
            for component, exponent in key
        )
        edges = dependency_edges(entries)
        records.append(
            {
                "kind": "rank1",
                "relation_line": relation,
                "support_size": len(entries),
                "cyclic": has_directed_cycle(edges),
                "edges": sorted(edges),
                "support": [monomial_label(entry) for entry in entries],
            }
        )

    records.sort(
        key=lambda record: (
            record["kind"],
            record.get("weight", record.get("relation_line")),
        )
    )
    return {
        "degree_bound": 4,
        "relation_count_with_labels": len(MONOMIALS),
        "distinct_relation_vectors": len(RELATIONS),
        "rank2_support_types": sum(r["kind"] == "rank2" for r in records),
        "rank2_cyclic_types": sum(
            r["kind"] == "rank2" and r["cyclic"] for r in records
        ),
        "rank1_support_types": sum(r["kind"] == "rank1" for r in records),
        "rank1_cyclic_types": sum(
            r["kind"] == "rank1" and r["cyclic"] for r in records
        ),
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    parser.add_argument("--cyclic-only", action="store_true")
    arguments = parser.parse_args()
    result = census()
    if arguments.cyclic_only:
        result["records"] = [
            record for record in result["records"] if record["cyclic"]
        ]
    rendered = json.dumps(result, indent=2)
    if arguments.json:
        arguments.json.write_text(rendered + "\n")
    else:
        print(rendered)


if __name__ == "__main__":
    main()

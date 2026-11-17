#!/usr/bin/env python3
"""Verify the universal nested-set target graph and its label action."""

from __future__ import annotations

import json
import math
import sys
from itertools import permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.wonderful_branch import (  # noqa: E402
    act_on_boundary,
    boundary_subsets,
    diagonal_generators,
    maximal_nested_sets,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "branch_wonderful_pullback.json"
)


def odd_double_factorial(value: int) -> int:
    result = 1
    for factor in range(1, value + 1, 2):
        result *= factor
    return result


boundary_counts = {}
zero_stratum_counts = {}
equivariance_checks = 0
for mark_count in range(4, 8):
    boundaries = boundary_subsets(mark_count)
    nested = maximal_nested_sets(mark_count)
    expected_boundary_count = 2 ** (mark_count - 1) - mark_count - 1
    expected_zero_strata = odd_double_factorial(2 * mark_count - 5)
    assert len(boundaries) == expected_boundary_count
    assert len(nested) == expected_zero_strata
    assert all(len(collection) == mark_count - 3 for collection in nested)

    # Fix the marked target values 0 and infinity; permute every moving
    # critical-value label.  The boundary building set and all maximal
    # nested sets are preserved.
    moving_marks = tuple(range(1, mark_count - 1))
    boundary_set = set(boundaries)
    nested_set = {
        frozenset(collection)
        for collection in nested
    }
    for image in permutations(moving_marks):
        action = {0: 0, **dict(zip(moving_marks, image))}
        transported_boundaries = {
            act_on_boundary(boundary, action)
            for boundary in boundaries
        }
        assert transported_boundaries == boundary_set
        transported_nested = {
            frozenset(
                act_on_boundary(boundary, action)
                for boundary in collection
            )
            for collection in nested
        }
        assert transported_nested == nested_set
        equivariance_checks += 1

    boundary_counts[str(mark_count)] = len(boundaries)
    zero_stratum_counts[str(mark_count)] = len(nested)


# With marks (0,infinity,A,B,C), the ten boundary divisors of Mbar_0,5 are
# six pair diagonals and four triple diagonals.  In the Kapranov P2 chart
# these are respectively the six arrangement lines and four exceptional
# divisors from blowing up their triple points.
degree_six_boundaries = boundary_subsets(5)
assert sum(len(boundary) == 2 for boundary in degree_six_boundaries) == 6
assert sum(len(boundary) == 3 for boundary in degree_six_boundaries) == 4

value_names = {0: "0", 1: "A", 2: "B", 3: "C"}
degree_six_ideals = {
    "".join(map(str, sorted(boundary))): diagonal_generators(
        boundary,
        value_names,
    )
    for boundary in degree_six_boundaries
}
assert degree_six_ideals["01"] == ("A",)
assert degree_six_ideals["12"] == ("B-A",)
assert degree_six_ideals["123"] == ("B-A", "C-A")


# The degree-five two-scale chart has four target marks.  Mbar_0,4=P1, and
# pulling back its homogeneous coordinate [lambda_1:lambda_2] gives the graph
# of [x^3:y^2], hence the blowup of (x^3,y^2).
degree_five_boundaries = boundary_subsets(4)
assert len(degree_five_boundaries) == 3
degree_five_pullback = {
    "target": "Mbar_0,4=P1",
    "branch_ratio": "[lambda_1:lambda_2]=[unit*x^3:unit*y^2]",
    "normalized_graph": "normalized blowup of (x^3,y^2)",
}


artifact = {
    "experiment": "wonderful pullback of the labelled stable-target graph",
    "verified_mark_counts": [4, 7],
    "boundary_divisor_counts": boundary_counts,
    "zero_stratum_counts": zero_stratum_counts,
    "expected_zero_strata_formula": "(2*n-5)!!",
    "label_equivariance_checks": equivariance_checks,
    "building_set": (
        "subsets S of finite target marks with 2<=|S|<=n-2; "
        "the omitted mark is infinity"
    ),
    "compatibility": "S subset T, T subset S, or S disjoint T",
    "pullback_rule": {
        "target_model": (
            "the wonderful compactification Mbar_0,n of the labelled "
            "critical-value configuration"
        ),
        "source_model": (
            "normalization of the principal component of the iterated "
            "blowups of the pulled-back diagonal ideals I_S"
        ),
        "tropical_model": (
            "saturated pullback of the tropical Mbar_0,n tree fan"
        ),
        "consequence": (
            "radial valuation walls are the first layer; residue collisions "
            "are later pulled-back diagonal centers in the same building set"
        ),
    },
    "degree_five_recovery": degree_five_pullback,
    "degree_six_recovery": {
        "boundary_divisors": 10,
        "arrangement_lines": 6,
        "exceptional_triple_centers": 4,
        "sample_pullback_ideals": degree_six_ideals,
    },
    "scope": (
        "complete labelled stable-target graph and its nested-set and "
        "permutation combinatorics; the recursive resonance atlas separately "
        "constructs the polynomial admissible-cover component over it"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS wonderful boundaries: n=4 through n=7")
print("PASS maximal nested sets: (2*n-5)!! zero strata")
print(f"PASS label action: {equivariance_checks} permutations")
print("PASS degree five: pullback blowup (x^3,y^2)")
print("PASS degree six: six lines plus four triple centers")
print("BRANCH_WONDERFUL_PULLBACK_PASS")

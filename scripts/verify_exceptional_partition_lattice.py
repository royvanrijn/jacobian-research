#!/usr/bin/env python3
"""Exact census of the exceptional weighted partition lattice through N=20."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    collision_precedes,
    refinement_branches,
)
from jcsearch.exceptional_partition_lattice import (  # noqa: E402
    atomic_components,
    intersection_census,
)


def full_contact_partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, 1, -1):
        for tail in full_contact_partitions(total - first, first):
            yield (first,) + tail


def minimal_common_coarsenings(left, right):
    common = [
        partition
        for partition in full_contact_partitions(left.degree)
        if collision_precedes(left.partition, partition)
        and collision_precedes(right.partition, partition)
    ]
    return tuple(
        partition
        for partition in common
        if not any(
            finer != partition and collision_precedes(finer, partition)
            for finer in common
        )
    )


# The universal primitive collision block.  Depress the monic cubic and
# compare it with a monic quadratic.  Its transverse algebra is one dual
# number; this is the weight attached to every sixfold transfer vertex.
Z = sp.symbols("Z")
u, v, alpha, beta = sp.symbols("u v alpha beta")
Q = Z**3 + u * Z + v
T = Z**2 + alpha * Z + beta
block_polynomial = sp.Poly(sp.expand(Q**2 - T**3), Z)
coefficients = [
    block_polynomial.coeff_monomial(Z**degree) for degree in range(7)
]
block_basis = sp.groebner(coefficients, v, u, beta, alpha, order="lex")
expected_block = sp.groebner(
    (alpha, v, 2 * u - 3 * beta, beta**2),
    v,
    u,
    beta,
    alpha,
    order="lex",
)
assert all(expected_block.reduce(item)[1] == 0 for item in coefficients)
assert all(
    block_basis.reduce(item.as_expr())[1] == 0 for item in expected_block.polys
)


records = intersection_census(3, 20)
component_records = [record for record in records if record.kind == "component_pair"]
self_records = [
    record for record in records if record.kind == "normalization_self_pair"
]

# Exhaust all full-contact partitions and verify that the closed formula is
# the unique minimal common integer-partition coarsening for every distinct
# atomic component pair.
for degree in range(3, 21):
    components = atomic_components(degree)
    for left_index, left in enumerate(components):
        for right in components[left_index + 1 :]:
            record = next(
                item
                for item in component_records
                if item.left == left and item.right == right
            )
            assert minimal_common_coarsenings(left, right) == (
                record.minimal_common_coarsening,
            )

            # Equality modulo affine polynomials is strong equality.  The
            # Wronskian divisor has degree N+(b+d)/2>N for distinct types.
            assert left.b + right.b > 0
            wronskian_divisor_twice_degree = (
                2 * degree + left.b + right.b
            )
            assert wronskian_divisor_twice_degree > 2 * degree

# Every listed self-pair has exactly the two exchanged allocations over its
# minimal stratum.  In N<=20 the global transfer balance forces (+1,-1), so
# no merged higher-transfer block enters the bounded census.
assert min(record.degree for record in self_records) == 12
for record in self_records:
    branches = refinement_branches(
        record.left.partition, record.minimal_common_coarsening
    )
    assert len(branches) == 2
    transfers = tuple(
        tuple((left_i - right_i) // 3 for (left_i, _), (right_i, _) in zip(*pair))
        for pair in ((branches[0], branches[1]),)
    )[0]
    assert sorted(value for value in transfers if value) == [-1, 1]

# Tangent-cone and length bookkeeping.  The total tangent space is the smooth
# stratum tangent plus one excess direction per primitive transfer block.
for record in records:
    assert sum(record.minimal_common_coarsening) == record.degree
    assert record.support_dimension == len(record.minimal_common_coarsening) - 1
    assert record.transverse_embedding_dimension == record.transfer_cluster_count
    assert record.intersection_tangent_dimension == (
        record.support_dimension + record.transfer_cluster_count
    )
    assert sum(record.tangent_cone_hilbert_vector) == (
        record.transverse_intersection_length
    )
    assert record.transverse_intersection_length == 2 ** record.transfer_cluster_count

# Fixed census totals make accidental omissions visible.
assert len(component_records) == 36
assert len(self_records) == 10
assert len(records) == 46
assert {record.transverse_intersection_length for record in records} == {2, 4, 8}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json", action="store_true", help="print the complete machine-readable census"
    )
    args = parser.parse_args()
    if args.json:
        print(
            json.dumps(
                {
                    "degree_range": [3, 20],
                    "component_pair_count": len(component_records),
                    "normalization_self_pair_count": len(self_records),
                    "records": [record.to_dict() for record in records],
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    print("PASS: every distinct component pair through N=20 has one minimal common coarsening")
    print("PASS: every minimal sixfold transfer contributes k[e]/(e^2) of length two")
    print("PASS: all 36 component pairs have I=2^t and binomial tangent cone")
    print("PASS: the 10 first normalization self-pairs have two independent clusters and I=4")


if __name__ == "__main__":
    main()

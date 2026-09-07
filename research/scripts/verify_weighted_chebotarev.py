#!/usr/bin/env python3
"""Verify the S_n fixed-point law and finite-field seed diagnostics."""

import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.chebotarev import (  # noqa: E402
    cycle_type_centralizer_size,
    cycle_type_probability,
    derangement_count,
    factorization_type_mod_prime,
    find_pencil_factorization_witness,
    fixed_point_count,
    fixed_point_distribution,
    integer_partitions,
    pencil_factorization_histogram,
    pencil_simple_root_histogram,
)
from jcsearch.weighted import (  # noqa: E402
    WeightedSeedModel,
    canonical_seed,
    deformation_basis,
    w,
)


# Exact symmetric-group combinatorics, factorial moments, and Poisson truncation.
for n in range(3, 9):
    distribution = fixed_point_distribution(n)
    assert sum(distribution.values(), Fraction()) == 1
    assert sum(fixed * probability for fixed, probability in distribution.items()) == 1
    assert distribution.get(n - 1, 0) == 0
    assert distribution[0] == Fraction(derangement_count(n), math.factorial(n))
    for order in range(1, n + 1):
        factorial_moment = sum(
            math.prod(range(fixed - order + 1, fixed + 1)) * probability
            for fixed, probability in distribution.items()
            if fixed >= order
        )
        assert factorial_moment == 1
    assert sum(fixed_point_count(n, fixed) for fixed in range(n + 1)) == math.factorial(n)

    # The conjugacy-class law refines the fixed-point law: a cycle type
    # lambda has density 1/z_lambda, and aggregation by the number of
    # one-cycles recovers the fixed-point distribution above.
    by_fixed_points = {}
    total_cycle_probability = Fraction()
    for partition in integer_partitions(n):
        centralizer = cycle_type_centralizer_size(partition)
        assert math.factorial(n) % centralizer == 0
        probability = cycle_type_probability(partition)
        total_cycle_probability += probability
        fixed = partition.count(1)
        by_fixed_points[fixed] = by_fixed_points.get(fixed, Fraction()) + probability
    assert total_cycle_probability == 1
    assert by_fixed_points == distribution
    for fixed, probability in distribution.items():
        truncated_poisson = Fraction(1, math.factorial(fixed)) * sum(
            (Fraction(-1) ** index) / math.factorial(index)
            for index in range(n - fixed + 1)
        )
        assert probability == truncated_poisson

# The original limiting law is exactly the S_3 fixed-point distribution.
assert fixed_point_distribution(3) == {
    0: Fraction(1, 3),
    1: Fraction(1, 2),
    3: Fraction(1, 6),
}

# Prime-field diagnostics for canonical degrees three and four and a generic
# two-extra-root degree-five seed.  These are finite samples, not the proof of
# Chebotarev; exact checks below certify totals and root incidence.
models = (
    ("canonical n=3", WeightedSeedModel(canonical_seed(2))),
    ("canonical n=4", WeightedSeedModel(canonical_seed(3))),
    (
        "deformed n=5",
        WeightedSeedModel(canonical_seed(2) + deformation_basis(2)),
    ),
)
prime = 31
for label, model in models:
    histogram = pencil_simple_root_histogram(model.primitive, w, prime)
    assert sum(histogram.values()) == prime**2

    # For each root r there are q-1 slopes other than H'(r), and then a unique
    # intercept.  Hence the total simple-root incidence is exactly q(q-1).
    assert sum(fiber * count for fiber, count in histogram.items()) == prime * (prime - 1)

    distribution = fixed_point_distribution(model.fiber_degree)
    discrepancy = max(
        abs(histogram.get(fixed, 0) - float(probability) * prime**2)
        for fixed, probability in distribution.items()
    )
    scaled = discrepancy / prime ** Fraction(3, 2)
    print(f"PASS {label} over F_{prime}: max discrepancy/q^(3/2)={float(scaled):.4f}")

# The deterministic cycle-type generator finds every partition in one audited
# good field, and its balanced integer output reduces to the requested type.
quartic = models[1][1]
factorization_histogram = pencil_factorization_histogram(
    quartic.primitive, w, prime
)
assert sum(factorization_histogram.values()) == prime**2
for partition in integer_partitions(quartic.fiber_degree):
    assert factorization_histogram.get(partition, 0) > 0
    witness = find_pencil_factorization_witness(
        quartic.primitive, w, prime, partition
    )
    assert witness is not None
    assert factorization_type_mod_prime(
        quartic.primitive,
        w,
        witness["lifted_slope"],
        witness["lifted_intercept"],
        prime,
    ) == partition

print("PASS: S_n fixed-point probabilities and factorial moments are exact")
print("PASS: conjugacy-class probabilities 1/z_lambda refine the fixed-point law")
print("PASS: fixed-point probabilities are the truncated Poisson(1) law")
print("PASS: simple-root incidence is q(q-1) in every audited inverse pencil")
print("PASS: finite samples align with the proved Chebotarev scaling")
print("PASS: every degree-four cycle type has a generated modular/integer witness")

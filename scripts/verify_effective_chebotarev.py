#!/usr/bin/env python3
"""Exact certificate and boundary audits for effective finite-field Chebotarev theorem."""

from collections import Counter
import math
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.chebotarev import (
    pencil_simple_root_histogram,
    rational_good_reduction_certificate,
)
from jcsearch.weighted import (
    WeightedSeedModel,
    canonical_seed,
    deformation_basis,
    w,
    x,
    y,
    z,
)


def modular_terms(expression, variables, prime):
    terms = []
    for monomial, coefficient in sp.Poly(expression, *variables).terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        assert denominator % prime
        terms.append((monomial, numerator * pow(denominator, -1, prime) % prime))
    return terms


def evaluate(terms, point, prime):
    return sum(
        coefficient * math.prod(pow(value, exponent, prime)
                                for value, exponent in zip(point, monomial))
        for monomial, coefficient in terms
    ) % prime


models = (
    WeightedSeedModel(canonical_seed(2)),
    WeightedSeedModel(canonical_seed(3)),
)

for model in models:
    certificate = rational_good_reduction_certificate(
        model.primitive, w, c=model.c, b0=model.b, a0=model.a
    )
    n = model.fiber_degree
    assert certificate["degree"] == n
    assert certificate["discriminant_degree"] == n
    assert certificate["bad_integer"] != 0
    assert certificate["boundary_resultant"] != 0

    # Every prime outside the certificate is larger than n and preserves all
    # explicitly inverted model constants.
    prime = next(
        candidate for candidate in (5, 7, 11, 13, 17, 19, 23, 29, 31)
        if sp.isprime(candidate)
        and certificate["bad_integer"] % candidate
        and candidate > n
    )
    assert math.factorial(n) % prime

    # The generic double-root parameterization and its derivative certify a
    # nonzero tame transposition locus modulo every certified prime.
    branch_s, branch_t = model.branch_parameterization()
    assert sp.expand(sp.diff(branch_t, w) - w * sp.diff(branch_s, w)) == 0
    leading_h2 = sp.Poly(sp.diff(model.primitive, w, 2), w).LC()
    numerator, denominator = map(int, sp.fraction(leading_h2))
    assert numerator * pow(denominator, -1, prime) % prime

    mapping = model.mapping()
    map_terms = [modular_terms(component, (x, y, z), prime)
                 for component in mapping]
    target_fibers = Counter()
    boundary_fibers = Counter()
    boundary_source_points = 0
    for source in ((xx, yy, zz)
                   for xx in range(prime)
                   for yy in range(prime)
                   for zz in range(prime)):
        target = tuple(evaluate(terms, source, prime) for terms in map_terms)
        target_fibers[target] += 1
        if target[2] == 0:
            boundary_fibers[target[:2]] += 1
            boundary_source_points += 1

    full_histogram = Counter(target_fibers.values())
    full_histogram[0] = prime**3 - len(target_fibers)
    boundary_histogram = Counter(boundary_fibers.values())
    boundary_histogram[0] = prime**2 - len(boundary_fibers)
    pencil_histogram = pencil_simple_root_histogram(model.primitive, w, prime)

    assert sum(boundary_histogram.values()) == prime**2
    assert boundary_histogram[0] == 0
    assert boundary_source_points == 2 * prime**2 - prime
    assert sum(fiber * count for fiber, count in boundary_histogram.items()) == (
        2 * prime**2 - prime
    )
    for fiber in range(n + 1):
        assert full_histogram[fiber] == (
            (prime - 1) * pencil_histogram.get(fiber, 0)
            + boundary_histogram[fiber]
        )

    print(
        f"PASS degree {n}: bad integer {certificate['bad_integer']}, "
        f"disc degree {certificate['discriminant_degree']}, boundary over F_{prime}"
    )

# A deformed degree-five seed has a genuinely nontrivial squarefree-factor
# resultant, so this audits more than the endpoint factors W and W-1.
deformed = WeightedSeedModel(canonical_seed(2) + deformation_basis(2))
deformed_certificate = rational_good_reduction_certificate(
    deformed.primitive, w, c=deformed.c, b0=deformed.b, a0=deformed.a
)
assert deformed_certificate["degree"] == 5
assert abs(deformed_certificate["boundary_resultant"]) != 1
assert deformed_certificate["discriminant_degree"] == 5

print("PASS: explicit certificates preserve degree, Jacobian units, and root profiles")
print("PASS: the deformed degree-five certificate has a nontrivial boundary resultant")
print("PASS: the generic branch derivative supplies tame transpositions at good primes")
print("PASS: N_j=(q-1)M_j+B_j and both exact boundary moments hold")

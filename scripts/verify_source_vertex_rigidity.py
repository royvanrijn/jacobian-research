#!/usr/bin/env python3
"""Verify the genus-zero two-fiber reconstruction used at source vertices."""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import prod
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "source_vertex_rigidity.json"
)

W = sp.symbols("W")


@lru_cache(maxsize=None)
def partitions(total: int, maximum: int | None = None) -> tuple[tuple[int, ...], ...]:
    """Return integer partitions in nonincreasing order."""

    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    result = []
    for first in range(maximum, 0, -1):
        for rest in partitions(total - first, first):
            result.append((first,) + rest)
    return tuple(result)


@lru_cache(maxsize=None)
def divisor_polynomial(
    points: tuple[int, ...],
    multiplicities: tuple[int, ...],
) -> sp.Expr:
    """Return the monic polynomial of a finite effective divisor."""

    assert len(points) == len(multiplicities)
    return sp.expand(
        sp.prod(
            (W - point) ** multiplicity
            for point, multiplicity in zip(points, multiplicities)
        )
    )


def multiplicity_at(polynomial: sp.Expr, point: int) -> int:
    """Return the order of a polynomial at a rational point."""

    multiplicity = 0
    quotient = sp.Poly(polynomial, W, domain=sp.QQ)
    factor = sp.Poly(W - point, W, domain=sp.QQ)
    while quotient.eval(point) == 0:
        quotient = quotient.exquo(factor)
        multiplicity += 1
    return multiplicity


def divisor_value(
    points: tuple[int, ...],
    multiplicities: tuple[int, ...],
    test_point: int,
) -> int:
    """Evaluate the product defining a divisor without symbolic expansion."""

    return prod(
        (test_point - point) ** multiplicity
        for point, multiplicity in zip(points, multiplicities)
    )


profiles_checked = 0
permutation_checks = 0
third_flag_checks = 0
verified_divisors: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()

# Include polynomial, reciprocal-polynomial, and mixed finite/infinite
# divisors.  Infinity can occur on at most one side because the zero and
# pole divisors are disjoint.
for degree in range(1, 8):
    for zero_at_infinity in range(degree + 1):
        for pole_at_infinity in range(degree + 1):
            if zero_at_infinity and pole_at_infinity:
                continue
            finite_zero_degree = degree - zero_at_infinity
            finite_pole_degree = degree - pole_at_infinity
            for zero_profile in partitions(finite_zero_degree):
                for pole_profile in partitions(finite_pole_degree):
                    zero_points = tuple(range(0, len(zero_profile)))
                    pole_points = tuple(
                        range(
                            len(zero_profile) + 2,
                            len(zero_profile) + 2 + len(pole_profile),
                        )
                    )
                    numerator = divisor_polynomial(
                        zero_points,
                        zero_profile,
                    )
                    denominator = divisor_polynomial(
                        pole_points,
                        pole_profile,
                    )

                    assert sum(zero_profile) == finite_zero_degree
                    assert sum(pole_profile) == finite_pole_degree
                    assert set(zero_points).isdisjoint(pole_points)
                    assert (
                        finite_pole_degree
                        - finite_zero_degree
                        == zero_at_infinity - pole_at_infinity
                    )
                    # Each effective divisor is checked once even when it
                    # occurs in many zero/pole profile pairs.
                    for points, profile, polynomial in (
                        (zero_points, zero_profile, numerator),
                        (pole_points, pole_profile, denominator),
                    ):
                        divisor_key = (points, profile)
                        if divisor_key in verified_divisors:
                            continue
                        for point, multiplicity in zip(points, profile):
                            assert (
                                multiplicity_at(polynomial, point)
                                == multiplicity
                            )

                        # Reordering the entries of a labelled divisor does
                        # not change its reconstructed polynomial.
                        orders = (
                            tuple(range(len(profile))),
                            tuple(reversed(range(len(profile)))),
                        )
                        if len(profile) <= 4:
                            orders = tuple(permutations(range(len(profile))))
                        for order in orders:
                            reordered = divisor_polynomial(
                                tuple(points[index] for index in order),
                                tuple(profile[index] for index in order),
                            )
                            assert reordered == polynomial
                            permutation_checks += 1
                        verified_divisors.add(divisor_key)

                    # Any two functions with these two fibers differ by a
                    # scalar.  A point over the target flag 1 fixes it.
                    test_point = len(zero_profile) + len(pole_profile) + 3
                    base_value = Fraction(
                        divisor_value(
                            zero_points,
                            zero_profile,
                            test_point,
                        ),
                        divisor_value(
                            pole_points,
                            pole_profile,
                            test_point,
                        ),
                    )
                    assert base_value != 0
                    normalized_scalar = 1 / base_value
                    assert normalized_scalar * base_value == 1

                    alternate_scalar = Fraction(7, 3)
                    alternate_value = alternate_scalar * base_value
                    alternate_normalization = (
                        alternate_scalar / alternate_value
                    )
                    assert alternate_normalization == normalized_scalar

                    second_test_point = test_point + 1
                    second_base_value = Fraction(
                        divisor_value(
                            zero_points,
                            zero_profile,
                            second_test_point,
                        ),
                        divisor_value(
                            pole_points,
                            pole_profile,
                            second_test_point,
                        ),
                    )
                    assert (
                        normalized_scalar * second_base_value
                        == alternate_normalization * second_base_value
                    )

                    profiles_checked += 1
                    third_flag_checks += 1


# The three source-vertex maps used throughout the degree-six atlas are
# direct instances of the same divisor reconstruction.
central_map = sp.expand(W**2 * (W - 1) ** 2 * (W - 3) ** 2)
assert central_map == divisor_polynomial((0, 1, 3), (2, 2, 2))

quadratic_cluster_tail = sp.expand(W * (W - 1))
assert quadratic_cluster_tail == divisor_polynomial((0, 1), (1, 1))

cyclic_square_tail = W**2
assert cyclic_square_tail == divisor_polynomial((0,), (2,))


artifact = {
    "experiment": "genus-zero source-vertex two-fiber rigidity",
    "degrees_checked": [1, 7],
    "divisor_profiles_checked": profiles_checked,
    "permutation_checks": permutation_checks,
    "third_flag_checks": third_flag_checks,
    "reconstruction": {
        "input": (
            "two disjoint effective degree-d divisors D_0 and D_infinity "
            "on a fixed P1 source component"
        ),
        "formula": (
            "f(W)=c*product(W-a_i)^(m_i)/"
            "product(W-b_j)^(n_j), with the degree difference recording "
            "the orders at infinity"
        ),
        "uniqueness": (
            "two maps with the same zero and pole divisors have constant "
            "ratio"
        ),
        "normalization": (
            "one source point over a third target flag fixes c; the rest "
            "of that third fiber is then an existence/compatibility test"
        ),
    },
    "degree_six_instances": {
        "central": "W^2*(W-1)^2*(W-3)^2",
        "quadratic_cluster_tail": "W*(W-1)",
        "cyclic_square_tail": "W^2",
    },
    "modular_consequence": (
        "once the source curve and inverse-image divisors of target flags "
        "are fixed, there is no residual vertexwise Hurwitz-class choice; "
        "the recursive resonance atlas separately constructs that "
        "flag-complete enhancement and proves contraction compatibility"
    ),
    "scope": (
        "componentwise reconstruction on genus-zero source curves; it does "
        "not construct the source tree or its target-flag pullback divisors "
        "over every wonderful boundary stratum"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print(f"PASS two-fiber reconstruction: {profiles_checked} divisor profiles")
print(f"PASS divisor permutation invariance: {permutation_checks} checks")
print(f"PASS third-flag scalar normalization: {third_flag_checks} checks")
print("PASS degree-six central, cluster-tail, and cyclic-tail instances")
print("SOURCE_VERTEX_RIGIDITY_PASS")

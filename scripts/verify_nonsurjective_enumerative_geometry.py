#!/usr/bin/env python3
"""Exact enumerative consequences of the 2/3 component decomposition."""

import math
import sys
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    collision_precedes,
    component_decomposition_count,
    maximal_two_three_partitions,
)


def full_contact_partitions(total, maximum=None):
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, 1, -1):
        for tail in full_contact_partitions(total - first, first):
            yield (first,) + tail


def atomic_pairs(degree):
    return tuple(
        (double_count, triple_count)
        for triple_count in range(degree // 3 + 1)
        if (degree - 3 * triple_count) % 2 == 0
        for double_count in ((degree - 3 * triple_count) // 2,)
    )


@lru_cache(maxsize=None)
def precedes(finer, coarser):
    return collision_precedes(finer, coarser)


# The component count is a period-six quasipolynomial.
for degree in range(121):
    count = len(atomic_pairs(degree))
    quotient, residue = divmod(degree, 6)
    expected = quotient if residue == 1 else quotient + 1
    assert count == expected
    if degree >= 2:
        assert count == len(maximal_two_three_partitions(degree))

x, u, v, t, y = sp.symbols("x u v t y")
count_series = sp.series(1 / ((1 - x**2) * (1 - x**3)), x, 0, 61).removeO()
for degree in range(61):
    assert sp.expand(count_series).coeff(x, degree) == len(atomic_pairs(degree))

# Dimensions and codimensions are consecutive staircases, and their refined
# series is the product of the two atomic geometric series.
refined = (
    1 / ((1 - u * v * x**2) * (1 - u * v**2 * x**3))
    - 1
    - u * v * x**2
) / (u * v**2)
refined = sp.series(refined, x, 0, 41).removeO().expand()
for degree in range(3, 41):
    pairs = atomic_pairs(degree)
    for index, (double_count, triple_count) in enumerate(pairs):
        dimension = double_count + triple_count - 1
        codimension = double_count + 2 * triple_count - 2
        assert dimension == degree // 2 - 1 - index
        assert codimension == (degree - 3) // 2 + index
        assert refined.coeff(x, degree).coeff(u, dimension).coeff(v, codimension) == 1
    assert sum(
        refined.coeff(x, degree).coeff(u, double_count + triple_count - 1).coeff(
            v, double_count + 2 * triple_count - 2
        )
        for double_count, triple_count in pairs
    ) == len(pairs)

# A common coarsening of two extreme components is automatically a
# coarsening of every intermediate component.  The maximal-support common
# type replaces each 3*2 <-> 2*3 exchange by one sixfold root.
INTERSECTION_CHECK_LIMIT = 22
interval_counts = defaultdict(lambda: defaultdict(int))
for degree in range(3, INTERSECTION_CHECK_LIMIT + 1):
    pairs = atomic_pairs(degree)
    partitions = tuple((3,) * b + (2,) * a for a, b in pairs)
    all_contacts = tuple(full_contact_partitions(degree))
    for left_index, (left_a, left_b) in enumerate(pairs):
        for right_index in range(left_index + 1, len(pairs)):
            right_a, right_b = pairs[right_index]
            gap = right_index - left_index
            assert right_a == left_a - 3 * gap
            assert right_b == left_b + 2 * gap
            top = tuple(sorted((2,) * right_a + (3,) * left_b + (6,) * gap, reverse=True))
            assert precedes(partitions[left_index], top)
            assert precedes(partitions[right_index], top)
            assert len(top) - 1 == left_a + left_b - 2 * gap - 1
            interval_counts[gap][degree] += 1

            for collision in all_contacts:
                common_to_extremes = precedes(
                    partitions[left_index], collision
                ) and precedes(partitions[right_index], collision)
                if not common_to_extremes:
                    continue
                for middle in partitions[left_index : right_index + 1]:
                    assert precedes(middle, collision)

            # The middle normalization fibers have the binomial profile.
            for offset in range(gap + 1):
                middle_a, middle_b = pairs[left_index + offset]
                branches = component_decomposition_count(top, middle_a, middle_b)
                assert branches == math.comb(gap, offset)

for gap in range(1, 6):
    gap_series = sp.series(
        x ** (6 * gap) / ((1 - x**2) * (1 - x**3)),
        x,
        0,
        INTERSECTION_CHECK_LIMIT + 1,
    ).removeO().expand()
    for degree in range(INTERSECTION_CHECK_LIMIT + 1):
        assert gap_series.coeff(x, degree) == interval_counts[gap][degree]

# Refined interval series: t records top-stratum dimension plus one and y
# records the component-index gap.
interval_series = (
    1
    / ((1 - t * x**2) * (1 - t * x**3))
    * (y * t * x**6 / (1 - y * t * x**6))
)
interval_series = sp.series(
    interval_series, x, 0, INTERSECTION_CHECK_LIMIT + 1
).removeO().expand()
for degree in range(3, INTERSECTION_CHECK_LIMIT + 1):
    pairs = atomic_pairs(degree)
    for left_index, (left_a, left_b) in enumerate(pairs):
        for right_index in range(left_index + 1, len(pairs)):
            gap = right_index - left_index
            top_dimension = left_a + left_b - 2 * gap - 1
            coefficient = (
                interval_series.coeff(x, degree)
                .coeff(y, gap)
                .coeff(t, top_dimension + 1)
            )
            # Several intervals can share the same degree, gap, and dimension;
            # compare against a direct count below rather than to one.
            direct = sum(
                1
                for start, (start_a, start_b) in enumerate(pairs)
                if start + gap < len(pairs)
                and start_a + start_b - 2 * gap - 1 == top_dimension
            )
            assert coefficient == direct

# The projective coincident-root locus has degree binomial(a+b,a)*2^a*3^b.
# Its degree/dimension series is 1/(1-2*t*x^2-3*t*x^3).
degree_series = sp.series(1 / (1 - 2 * t * x**2 - 3 * t * x**3), x, 0, 31)
degree_series = degree_series.removeO().expand()
for degree in range(31):
    for double_count, triple_count in atomic_pairs(degree):
        support = double_count + triple_count
        projective_degree = (
            math.comb(support, double_count)
            * 2**double_count
            * 3**triple_count
        )
        assert degree_series.coeff(x, degree).coeff(t, support) == projective_degree

# Chow codimension series for one hyperplane incidence.
chow_series = sp.series(y / (1 - 2 * y * x**2 - 3 * y**2 * x**3), x, 0, 31)
chow_series = chow_series.removeO().expand()
for degree in range(31):
    for double_count, triple_count in atomic_pairs(degree):
        support = double_count + triple_count
        projective_degree = (
            math.comb(support, double_count)
            * 2**double_count
            * 3**triple_count
        )
        codimension = double_count + 2 * triple_count + 1
        assert chow_series.coeff(x, degree).coeff(y, codimension) == projective_degree

# Source multidegrees of the pulled-back incidence divisor.
for double_count in range(5):
    for triple_count in range(5):
        support = double_count + triple_count
        if support == 0:
            continue
        for p in range(double_count + 1):
            for q in range(triple_count + 1):
                r = support - 1 - p - q
                if r < 0:
                    continue
                expression = sp.expand(
                    u**p * v**q * (2 * u + 3 * v) ** (r + 1)
                )
                actual = expression.coeff(u, double_count).coeff(v, triple_count)
                expected = (
                    math.comb(r + 1, double_count - p)
                    * 2 ** (double_count - p)
                    * 3 ** (triple_count - q)
                )
                assert actual == expected

print("PASS: component counts have the period-six quasipolynomial and rational series")
print("PASS: component dimensions and codimensions form the exact staircases")
print("PASS: component intersections depend only on their extreme indices")
print("PASS: interval intersections have the x^6-shifted generating functions")
print("PASS: top collision fibers have the binomial branch profile")
print("PASS: coincident-root degrees, Chow classes, and multidegrees agree")

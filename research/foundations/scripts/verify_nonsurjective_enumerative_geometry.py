#!/usr/bin/env python3
"""Exact enumerative consequences of the 2/3 component decomposition."""

import math
import sys
from collections import defaultdict
from functools import lru_cache
from itertools import product
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    collision_precedes,
    component_decomposition_count,
    contact_partition_incidence,
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

# The retained-coefficient projection has one boundary base point.  Its local
# Segre correction is the weighted complete-intersection length
# (support+1)!/(a!b!), giving the exact normalized seed-component degree.
corrected_degree_series = sp.series(
    1 / (1 - 2 * t * x**2 - 3 * t * x**3)
    - 1 / (1 - t * x**2 - t * x**3) ** 2,
    x,
    0,
    31,
).removeO().expand()
for degree in range(31):
    for double_count, triple_count in atomic_pairs(degree):
        support = double_count + triple_count
        correction = (support + 1) * math.comb(support, double_count)
        component_degree = math.comb(support, double_count) * (
            2**double_count * 3**triple_count - support - 1
        )
        assert correction == math.factorial(support + 1) // (
            math.factorial(double_count) * math.factorial(triple_count)
        )
        assert (
            corrected_degree_series.coeff(x, degree).coeff(t, support)
            == component_degree
        )

assert {
    pair: math.comb(sum(pair), pair[0])
    * (2 ** pair[0] * 3 ** pair[1] - sum(pair) - 1)
    for pair in ((0, 1), (2, 0), (1, 1), (3, 0), (0, 2))
} == {(0, 1): 1, (2, 0): 1, (1, 1): 6, (3, 0): 4, (0, 2): 6}

# All graph multidegrees involving a source hyperplane avoid the unique base
# point and retain their raw values.  Only the top target degree is corrected.
A, B, S, T, Z, epsilon = sp.symbols("A B S T Z epsilon")
graph_multidegree_series = (
    (2 * A + 3 * B)
    / ((1 - A * S) * (1 - B * T) * (1 - Z * (2 * A + 3 * B)))
    - ((1 - Z * (A + B)) ** -2 - 1) / Z
)
GRAPH_SUPPORT_CHECK = 5
truncated_graph_series = sp.series(
    graph_multidegree_series.subs({A: epsilon * A, B: epsilon * B}),
    epsilon,
    0,
    GRAPH_SUPPORT_CHECK + 1,
).removeO().expand()
for support in range(1, GRAPH_SUPPORT_CHECK + 1):
    homogeneous_piece = truncated_graph_series.coeff(epsilon, support)
    for double_count in range(support + 1):
        triple_count = support - double_count
        actual = homogeneous_piece.coeff(A, double_count).coeff(B, triple_count)
        expected = 0
        for p in range(double_count + 1):
            for q in range(triple_count + 1):
                r = support - 1 - p - q
                if r < 0:
                    continue
                raw = (
                    math.comb(r + 1, double_count - p)
                    * 2 ** (double_count - p)
                    * 3 ** (triple_count - q)
                )
                correction = (
                    (support + 1) * math.comb(support, double_count)
                    if p == q == 0
                    else 0
                )
                expected += (raw - correction) * S**p * T**q * Z**r
        assert sp.expand(actual - expected) == 0


def standard_monomial_count(groebner_basis):
    """Length of a zero-dimensional quotient from its leading monomials."""
    leading = tuple(
        polynomial.LM(order=groebner_basis.order).exponents
        for polynomial in groebner_basis.polys
    )
    variable_count = len(groebner_basis.gens)
    bounds = []
    for index in range(variable_count):
        pure_powers = [
            monomial[index]
            for monomial in leading
            if monomial[index] > 0
            and all(
                exponent == 0
                for other_index, exponent in enumerate(monomial)
                if other_index != index
            )
        ]
        assert pure_powers
        bounds.append(min(pure_powers))
    return sum(
        1
        for exponents in product(*(range(bound) for bound in bounds))
        if not any(
            all(
                exponents[index] >= monomial[index]
                for index in range(variable_count)
            )
            for monomial in leading
        )
    )


# Direct Groebner regressions for every local weighted system with at most
# four variables.  These are finite checks of the Hilbert-series calculation,
# not a replacement for the uniform regular-sequence proof.
W = sp.symbols("W_local_degree_check")
for support in range(1, 5):
    for double_count in range(support + 1):
        triple_count = support - double_count
        degree = 2 * double_count + 3 * triple_count
        if degree < 3:
            continue
        q_variables = sp.symbols(f"q1:{double_count + 1}") if double_count else ()
        r_variables = sp.symbols(f"r1:{triple_count + 1}") if triple_count else ()
        local_variables = q_variables + r_variables
        Q = 1 + sum(
            variable * W**index
            for index, variable in enumerate(q_variables, start=1)
        )
        R = 1 + sum(
            variable * W**index
            for index, variable in enumerate(r_variables, start=1)
        )
        M = sp.expand(Q**2 * R**3)
        local_system = [
            M.coeff(W, coefficient_degree)
            for coefficient_degree in range(2, support + 2)
        ]
        basis = sp.groebner(local_system, *local_variables, order="grevlex")
        assert basis.is_zero_dimensional
        assert standard_monomial_count(basis) == math.factorial(support + 1) // (
            math.factorial(double_count) * math.factorial(triple_count)
        )

# The local correction is independent of the atom multiplicities.  Test the
# universal product M=prod Q_nu**e_nu for several exponent/degree profiles.
for exponents, factor_degrees in (
    ((4,), (2,)),
    ((2, 4), (1, 1)),
    ((3, 5), (2, 1)),
    ((2, 3, 5), (1, 1, 1)),
):
    factors = []
    generalized_variables = ()
    for factor_index, factor_degree in enumerate(factor_degrees):
        variables = sp.symbols(
            f"general_factor_{factor_index}_1:{factor_degree + 1}"
        )
        generalized_variables += variables
        factors.append(
            1
            + sum(
                variable * W**index
                for index, variable in enumerate(variables, start=1)
            )
        )
    generalized_M = sp.expand(
        sp.prod(
            factor**exponent for factor, exponent in zip(factors, exponents)
        )
    )
    generalized_support = sum(factor_degrees)
    generalized_system = [
        generalized_M.coeff(W, coefficient_degree)
        for coefficient_degree in range(2, generalized_support + 2)
    ]
    generalized_basis = sp.groebner(
        generalized_system, *generalized_variables, order="grevlex"
    )
    predicted_length = math.factorial(generalized_support + 1)
    for factor_degree in factor_degrees:
        predicted_length //= math.factorial(factor_degree)
    assert generalized_basis.is_zero_dimensional
    assert standard_monomial_count(generalized_basis) == predicted_length

# Independent projection regression: the degree-six (3,3) normalization is a
# rational curve.  Clearing denominators in [1:h_3:h_4:h_5] gives four
# coprime forms of common homogenized degree six, exactly as (22) predicts.
three_three = contact_partition_incidence(6, (3, 3))
q1, q2 = three_three.quotient_coordinates
q = sp.symbols("q_projective_degree_check")
parameter = {
    q1: 1 - 3 * q**2 / (3 * q - 1),
    q2: -3 * q**3 / (3 * q - 1),
}
primitive = sp.Poly(three_three.primitive, three_three.variable)
affine_seed_coordinates = [
    sp.cancel(primitive.coeff_monomial(three_three.variable**power).subs(parameter))
    for power in range(3, 6)
]
numerators_and_denominators = [
    sp.fraction(coordinate) for coordinate in affine_seed_coordinates
]
common_denominator = sp.lcm(
    [denominator for _, denominator in numerators_and_denominators]
)
projective_coordinates = [common_denominator] + [
    sp.cancel(numerator * common_denominator / denominator)
    for numerator, denominator in numerators_and_denominators
]
common_factor = sp.gcd_list(projective_coordinates)
projective_coordinates = [
    sp.cancel(coordinate / common_factor) for coordinate in projective_coordinates
]
assert sp.gcd_list(projective_coordinates) == 1
assert max(sp.degree(coordinate, q) for coordinate in projective_coordinates) == 6

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
print("PASS: weighted boundary corrections give the exact projective component degrees")
print("PASS: local complete-intersection lengths agree through support four")
print("PASS: the projected degree-six (3,3) normalization is a sextic curve")
print("PASS: the full graph-multidegree series agrees through support five")
print("PASS: the local correction is universal across tested multiplicity atoms")

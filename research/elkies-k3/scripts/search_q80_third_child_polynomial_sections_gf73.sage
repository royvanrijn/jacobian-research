#!/usr/bin/env sage
"""Search the lattice-predicted polynomial sections on the CM24 third child.

The exact model has fibers ``2I7+3I2+4I1``.  The saturated optimal MW basis
from ``analyze_q80_fourth_q12_cm24_marking.sage`` has component corrections

* P1: A1 support (0,1,1), A6 corrections (10/7,10/7);
* P2: A1 support (1,0,1), A6 corrections (12/7,6/7);
* P3: A1 support (0,0,0), A6 corrections (10/7,10/7).

Thus P1 and P2 are polynomial sections whose x-coordinate passes through
both I7 nodes and two I2 nodes.  Four node values leave only the leading
coefficient of a quartic free, so all candidates require just ``3*73`` exact
square tests.  ``--two-node`` additionally scans the three-parameter family
for the P3 support; it is bounded and optional.

Node incidence does not distinguish the two orientations of an In component.
Candidates emitted here must still be matched by the exact component labels
and MW pairings before they define the fourth-neighbor marking.
"""

import argparse
from itertools import combinations, product

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


parser = argparse.ArgumentParser()
parser.add_argument("--two-node", action="store_true")
parser.add_argument("--match-mw", action="store_true")
arguments = parser.parse_args()
if arguments.match_mw:
    arguments.two_node = True

finite = GF(73)
polynomial_ring = PolynomialRing(finite, "V")
V = polynomial_ring.gen()
fiber_ring = PolynomialRing(finite, "x")
x = fiber_ring.gen()

A = (
    6*V**8 + 16*V**7 + 47*V**6 + 33*V**5 + 58*V**4
    + 2*V**3 + 63*V**2 + 17*V + 23
)
B = (
    33*V**12 + 64*V**10 + 61*V**9 + 45*V**8 + 14*V**7
    + 20*V**6 + 54*V**5 + 8*V**4 + 50*V**3 + 57*V**2
    + 47*V + 43
)
Delta = -finite(16)*(4*A**3+27*B**2)

i7_roots = (finite(-20), finite(-67))
i2_roots = (finite(-17), finite(-30), finite(-68))
for root in i7_roots:
    assert Delta.valuation(V-root) == 7
for root in i2_roots:
    assert Delta.valuation(V-root) == 2


def node_x(root):
    cubic = x**3 + A(root)*x + B(root)
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    return -common[0]/common[1]


nodes = {root: node_x(root) for root in i7_roots + i2_roots}


def polynomial_square_roots(polynomial):
    """Return the polynomial square roots of degree at most six."""
    assert polynomial.degree() <= 12
    if polynomial == 0:
        return (polynomial_ring.zero(),)
    shift = next(value for value in finite if polynomial(value) != 0)
    shifted = polynomial(V+shift)
    constant = shifted[0]
    if not constant.is_square():
        return ()
    roots = []
    for first in constant.sqrt(all=True):
        coefficients = [first]
        for degree in range(1, 7):
            known = sum(
                coefficients[left]*coefficients[degree-left]
                for left in range(1, degree)
            )
            coefficients.append((shifted[degree]-known)/(2*first))
        candidate_shifted = polynomial_ring(coefficients)
        if candidate_shifted**2 == shifted:
            roots.append(candidate_shifted(V-shift))
    return tuple(roots)


def constrained_x_family(selected_roots):
    interpolation = polynomial_ring.lagrange_polynomial(
        [(root, nodes[root]) for root in selected_roots]
    )
    vanishing = polynomial_ring.one()
    for root in selected_roots:
        vanishing *= V-root
    assert interpolation.degree() < len(selected_roots)
    return interpolation, vanishing


four_node_candidates = []
for i2_subset in combinations(i2_roots, 2):
    selected_roots = i7_roots + i2_subset
    interpolation, vanishing = constrained_x_family(selected_roots)
    assert vanishing.degree() == 4
    for leading in finite:
        X = interpolation + leading*vanishing
        square_roots = polynomial_square_roots(X**3+A*X+B)
        for Y in square_roots:
            assert Y**2 == X**3+A*X+B
            four_node_candidates.append(
                (tuple(int(root) for root in i2_subset), X, Y)
            )

print(
    "Q80THIRDPOLYSEARCH|prime=73|"
    f"I7_nodes={tuple((int(root), int(nodes[root])) for root in i7_roots)}|"
    f"I2_nodes={tuple((int(root), int(nodes[root])) for root in i2_roots)}|"
    f"four_node_tests={3*73}|four_node_sections={len(four_node_candidates)}",
    flush=True,
)
for i2_subset, X, Y in four_node_candidates:
    print(
        f"Q80THIRDPOLYSEARCH|I2_support={i2_subset}|"
        f"X={tuple(map(int, X.list()))}|Y={tuple(map(int, Y.list()))}",
        flush=True,
    )

two_node_candidates = []
if arguments.two_node:
    interpolation, vanishing = constrained_x_family(i7_roots)
    assert vanishing.degree() == 2
    for coefficients in product(finite, repeat=3):
        X = interpolation + polynomial_ring(coefficients)*vanishing
        square_roots = polynomial_square_roots(X**3+A*X+B)
        for Y in square_roots:
            assert Y**2 == X**3+A*X+B
            two_node_candidates.append((X, Y))
    print(
        "Q80THIRDPOLYSEARCH|"
        f"two_node_tests={73**3}|two_node_sections={len(two_node_candidates)}",
        flush=True,
    )
    for X, Y in two_node_candidates:
        actual_i2_support = tuple(
            int(root) for root in i2_roots if X(root) == nodes[root]
        )
        print(
            f"Q80THIRDPOLYSEARCH|I2_support={actual_i2_support}|"
            f"X={tuple(map(int, X.list()))}|Y={tuple(map(int, Y.list()))}",
            flush=True,
        )


if arguments.match_mw:
    function_field = polynomial_ring.fraction_field()
    curve = EllipticCurve(
        function_field,
        [0, 0, 0, function_field(A), function_field(B)],
    )
    points = tuple(
        curve(function_field(X), function_field(Y))
        for X, Y in two_node_candidates
    )
    assert len(set(points)) == len(points) == 30

    def hits_node(point, root):
        if point.is_zero():
            return False
        x_coordinate = point[0]
        y_coordinate = point[1]
        if x_coordinate.denominator()(root) == 0:
            return False
        return (
            x_coordinate(root) == nodes[root]
            and y_coordinate(root) == 0
        )

    def section_pole(point):
        assert not point.is_zero()
        x_coordinate = point[0]
        numerator_degree = x_coordinate.numerator().degree()
        denominator_degree = x_coordinate.denominator().degree()
        assert denominator_degree % 2 == 0
        infinity_excess = max(
            0, numerator_degree-denominator_degree-4
        )
        assert infinity_excess % 2 == 0
        return denominator_degree//2 + infinity_excess//2

    # Choose one nonidentity point as label 1 independently at each I7.
    # The component map is a homomorphism, so Q-k*reference is on the
    # identity component for exactly one k modulo seven.
    reference = points[0]

    def i7_label(point, root):
        answers = tuple(
            multiplier
            for multiplier in range(7)
            if not hits_node(point-multiplier*reference, root)
        )
        assert len(answers) == 1
        return answers[0]

    def labels(point):
        return (
            tuple(i7_label(point, root) for root in i7_roots),
            tuple(1 if hits_node(point, root) else 0 for root in i2_roots),
        )

    raw_point_labels = tuple(labels(point) for point in points)

    def self_correction(point_labels_row):
        i7, i2 = point_labels_row
        return sum(QQ(label*(7-label))/7 for label in i7) + QQ(sum(i2))/2

    def local_pairing(left_labels, right_labels):
        left_i7, left_i2 = left_labels
        right_i7, right_i2 = right_labels
        answer = sum(
            QQ(min(left, right)*(7-max(left, right)))/7
            for left, right in zip(left_i7, right_i7)
        )
        answer += QQ(sum(left*right for left, right in zip(left_i2, right_i2)))/2
        return answer

    target_heights = (QQ(1)/7, QQ(3)/7, QQ(8)/7)
    scaled_matches = []
    for scales in product(range(1, 7), repeat=2):
        point_labels = tuple(
            (
                tuple(
                    (scale*label) % 7
                    for scale, label in zip(scales, i7_labels)
                ),
                i2_labels,
            )
            for i7_labels, i2_labels in raw_point_labels
        )
        heights = tuple(
            QQ(4) + 2*section_pole(point) - self_correction(label_row)
            for point, label_row in zip(points, point_labels)
        )
        by_height = tuple(
            tuple(
                index for index, height in enumerate(heights)
                if height == target
            )
            for target in target_heights
        )

        def height_pairing(left_index, right_index):
            assert left_index != right_index
            intersection_value = section_pole(
                points[left_index]-points[right_index]
            )
            return (
                QQ(2)-intersection_value
                - local_pairing(
                    point_labels[left_index], point_labels[right_index]
                )
            )

        for first in by_height[0]:
            for second in by_height[1]:
                if height_pairing(first, second) != QQ(1)/14:
                    continue
                for third in by_height[2]:
                    if height_pairing(first, third) != 0:
                        continue
                    if height_pairing(second, third) != -QQ(1)/7:
                        continue
                    first_i2 = point_labels[first][1]
                    second_i2 = point_labels[second][1]
                    third_i2 = point_labels[third][1]
                    if (
                        sum(first_i2), sum(second_i2), sum(third_i2)
                    ) != (2, 2, 0):
                        continue
                    if sum(
                        a*b for a, b in zip(first_i2, second_i2)
                    ) != 1:
                        continue
                    scaled_matches.append(
                        (scales, point_labels, by_height, first, second, third)
                    )

    assert scaled_matches
    unique_horizontals = {}
    for _, match_labels, _, match_first, match_second, match_third in scaled_matches:
        candidate_horizontal = (
            -points[match_first]-3*points[match_second]-points[match_third]
        )
        candidate_i7_labels = tuple(
            (
                -match_labels[match_first][0][index]
                -3*match_labels[match_second][0][index]
                -match_labels[match_third][0][index]
            ) % 7
            for index in range(2)
        )
        candidate_i2_labels = tuple(
            (
                -match_labels[match_first][1][index]
                -3*match_labels[match_second][1][index]
                -match_labels[match_third][1][index]
            ) % 2
            for index in range(3)
        )
        key = (
            tuple(candidate_horizontal[0].numerator().list()),
            tuple(candidate_horizontal[0].denominator().list()),
            tuple(candidate_horizontal[1].numerator().list()),
            tuple(candidate_horizontal[1].denominator().list()),
        )
        unique_horizontals[key] = (
            candidate_horizontal,
            (candidate_i7_labels, candidate_i2_labels),
        )
    assert len(unique_horizontals) == 4
    scales, point_labels, by_height, first, second, third = scaled_matches[0]
    assert scales == (1, 3)
    assert (first, second, third) == (20, 0, 2)
    assert (
        point_labels[first], point_labels[second], point_labels[third]
    ) == (
        ((5, 2), (1, 0, 1)),
        ((1, 3), (1, 1, 0)),
        ((5, 5), (0, 0, 0)),
    )
    horizontal = -points[first]-3*points[second]-points[third]
    assert section_pole(horizontal) == 2

    def rational_coefficients(value):
        return (
            tuple(map(int, value.numerator().list())),
            tuple(map(int, value.denominator().list())),
        )

    assert rational_coefficients(horizontal[0]) == (
        (38, 69, 69, 7, 45, 4, 4, 3, 70),
        (35, 9, 56, 47, 1),
    )
    assert rational_coefficients(horizontal[1]) == (
        (50, 62, 5, 24, 56, 68, 27, 62, 46, 38, 56, 40, 34),
        (30, 22, 39, 20, 9, 34, 1),
    )

    print(
        "Q80THIRDPOLYSEARCH|"
        f"height_counts={tuple(len(indices) for indices in by_height)}|"
        f"label_scales={scales}|marked_basis_matches={len(scaled_matches)}|"
        f"chosen_indices={(first, second, third)}|"
        f"chosen_labels={(point_labels[first], point_labels[second], point_labels[third])}",
        flush=True,
    )
    for candidate_index, (candidate_horizontal, candidate_labels) in enumerate(
        unique_horizontals.values(), 1
    ):
        print(
            f"Q80THIRDPOLYSEARCH|horizontal_candidate={candidate_index}|"
            f"component_labels={candidate_labels}|"
            f"X_num_den={rational_coefficients(candidate_horizontal[0])}|"
            f"Y_num_den={rational_coefficients(candidate_horizontal[1])}",
            flush=True,
        )
    print(
        "Q80THIRDPOLYSEARCH|horizontal=-P1-3P2-P3|height=33/7|P.O=2|"
        f"X_num_den={rational_coefficients(horizontal[0])}|"
        f"Y_num_den={rational_coefficients(horizontal[1])}|"
        "status=PASS_EXACT_MW_MARKING",
        flush=True,
    )

assert four_node_candidates
print(
    "Q80THIRDPOLYSEARCH|status=PASS_EXACT_NODE_CONSTRAINED_SEARCH",
    flush=True,
)

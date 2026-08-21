#!/usr/bin/env sage
"""Test the fourth q12 horizontal section at the smooth fixed-j residue.

The two split-prime embeddings of the CM24 third child have equal j only at
``V=-27`` away from their discriminants.  The original embedding already has
the marked fourth-neighbor horizontal section ``Q=-P1-3P2-P3``.  This script
independently enumerates the thirty polynomial sections on the conjugate
embedding, recovers every Mordell--Weil marking compatible with the exact
height Gram, and asks whether the x-coordinate of the marked fourth section
matches after the unique short-Weierstrass quadratic twist at the fixed-j
fiber.

This is an exact characteristic-73 descent diagnostic.  A match would be a
specialization seed, not a characteristic-zero descent or a rational curve.
"""

from itertools import product

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


finite = GF(73)
base = PolynomialRing(finite, "V")
V = base.gen()
fiber_ring = PolynomialRing(finite, "x")
x = fiber_ring.gen()

original_A = (
    6*V**8 + 16*V**7 + 47*V**6 + 33*V**5 + 58*V**4
    + 2*V**3 + 63*V**2 + 17*V + 23
)
original_B = (
    33*V**12 + 64*V**10 + 61*V**9 + 45*V**8 + 14*V**7
    + 20*V**6 + 54*V**5 + 8*V**4 + 50*V**3 + 57*V**2
    + 47*V + 43
)
conjugate_A = (
    54*V**8 + 69*V**7 + V**6 + 70*V**5 + 2*V**4
    + 35*V**3 + 15*V**2 + 14*V + 56
)
conjugate_B = (
    26*V**12 + 5*V**11 + 10*V**10 + 49*V**9 + 64*V**8
    + 54*V**7 + 4*V**6 + 39*V**5 + 10*V**4 + 69*V**3
    + 23*V**2 + 49*V + 56
)

# The convenient short normalization reconstructed from the conjugate j-map
# is a constant nonsquare twist of the model carrying the CM sections.  Undo
# that harmless normalization twist before enumerating the section scheme.
normalization_twist = next(value for value in finite if value and not value.is_square())
marked_conjugate_A = normalization_twist**2*conjugate_A
marked_conjugate_B = normalization_twist**3*conjugate_B


def polynomial_square_roots(polynomial):
    """Return all polynomial square roots of degree at most six."""
    assert polynomial.degree() <= 12
    if polynomial == 0:
        return (base.zero(),)
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
        candidate_shifted = base(coefficients)
        if candidate_shifted**2 == shifted:
            roots.append(candidate_shifted(V-shift))
    return tuple(roots)


def fiber_nodes(A, B):
    delta = -finite(16)*(4*A**3+27*B**2)
    i7_roots = []
    i2_roots = []
    for factor, multiplicity in delta.factor():
        if factor.degree() != 1:
            continue
        root = -factor[0]/factor[1]
        if multiplicity == 7:
            i7_roots.append(root)
        elif multiplicity == 2:
            i2_roots.append(root)
    assert (len(i7_roots), len(i2_roots)) == (2, 3)

    def node_x(root):
        cubic = x**3+A(root)*x+B(root)
        common = cubic.gcd(cubic.derivative())
        assert common.degree() == 1
        return -common[0]/common[1]

    roots = tuple(i7_roots), tuple(i2_roots)
    nodes = {root: node_x(root) for root in sum(roots, ())}
    return roots[0], roots[1], nodes


def polynomial_sections(A, B):
    """Enumerate all quartic-x polynomial sections through both I7 nodes."""
    i7_roots, i2_roots, nodes = fiber_nodes(A, B)
    interpolation = base.lagrange_polynomial(
        [(root, nodes[root]) for root in i7_roots]
    )
    vanishing = base.one()
    for root in i7_roots:
        vanishing *= V-root
    sections = []
    for coefficients in product(finite, repeat=3):
        X = interpolation+base(coefficients)*vanishing
        for Y in polynomial_square_roots(X**3+A*X+B):
            sections.append((X, Y))
    if len(sections) != 30:
        raise AssertionError(
            f"expected 30 polynomial sections, found {len(sections)}"
        )
    return i7_roots, i2_roots, nodes, tuple(sections)


def marked_fourth_sections(A, B):
    """Recover all marking-compatible copies of -P1-3P2-P3."""
    i7_roots, i2_roots, nodes, sections = polynomial_sections(A, B)
    function_field = base.fraction_field()
    curve = EllipticCurve(
        function_field, [0, 0, 0, function_field(A), function_field(B)]
    )
    points = tuple(
        curve(function_field(X), function_field(Y)) for X, Y in sections
    )
    assert len(set(points)) == 30

    def hits_node(point, root):
        if point.is_zero():
            return False
        x_coordinate, y_coordinate = point[0], point[1]
        return (
            x_coordinate.denominator()(root) != 0
            and x_coordinate(root) == nodes[root]
            and y_coordinate(root) == 0
        )

    def section_pole(point):
        assert not point.is_zero()
        x_coordinate = point[0]
        numerator_degree = x_coordinate.numerator().degree()
        denominator_degree = x_coordinate.denominator().degree()
        infinity_excess = max(
            0, numerator_degree-denominator_degree-4
        )
        assert denominator_degree % 2 == infinity_excess % 2 == 0
        return denominator_degree//2+infinity_excess//2

    reference = points[0]

    def raw_labels(point):
        labels = []
        for root in i7_roots:
            answers = tuple(
                multiplier for multiplier in range(7)
                if not hits_node(point-multiplier*reference, root)
            )
            assert len(answers) == 1
            labels.append(answers[0])
        return tuple(labels), tuple(
            1 if hits_node(point, root) else 0 for root in i2_roots
        )

    raw_point_labels = tuple(raw_labels(point) for point in points)

    def self_correction(labels):
        i7, i2 = labels
        return sum(QQ(label*(7-label))/7 for label in i7)+QQ(sum(i2))/2

    def local_pairing(left, right):
        left_i7, left_i2 = left
        right_i7, right_i2 = right
        answer = sum(
            QQ(min(a, b)*(7-max(a, b)))/7
            for a, b in zip(left_i7, right_i7)
        )
        return answer+QQ(sum(a*b for a, b in zip(left_i2, right_i2)))/2

    target_heights = (QQ(1)/7, QQ(3)/7, QQ(8)/7)
    marked = []
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
            QQ(4)+2*section_pole(point)-self_correction(labels)
            for point, labels in zip(points, point_labels)
        )
        by_height = tuple(
            tuple(index for index, height in enumerate(heights) if height == target)
            for target in target_heights
        )

        def height_pairing(left_index, right_index):
            return (
                QQ(2)-section_pole(points[left_index]-points[right_index])
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
                    supports = tuple(
                        point_labels[index][1]
                        for index in (first, second, third)
                    )
                    if tuple(map(sum, supports)) != (2, 2, 0):
                        continue
                    if sum(a*b for a, b in zip(supports[0], supports[1])) != 1:
                        continue
                    horizontal = (
                        -points[first]-3*points[second]-points[third]
                    )
                    assert section_pole(horizontal) == 2
                    marked.append((scales, (first, second, third), horizontal))
    assert marked
    unique = {}
    for scales, indices, point in marked:
        unique.setdefault(point, (scales, indices))
    return len(marked), tuple(
        (point, unique[point][0], unique[point][1]) for point in unique
    )


original_matches, original_fourths = marked_fourth_sections(
    original_A, original_B
)
conjugate_matches, conjugate_fourths = marked_fourth_sections(
    marked_conjugate_A, marked_conjugate_B
)

fixed_base = finite(-27)
assert -finite(16)*(
    4*original_A(fixed_base)**3+27*original_B(fixed_base)**2
) != 0
assert -finite(16)*(
    4*marked_conjugate_A(fixed_base)**3
    + 27*marked_conjugate_B(fixed_base)**2
) != 0

# If E': y^2=x^3+A'x+B' is the quadratic twist of E by d in this
# normalization, then A'=d^2*A and B'=d^3*B.  The x-coordinate map over a
# square root of d is x' = d*x, independent of the sign of that root.
a_ratio = marked_conjugate_A(fixed_base)/original_A(fixed_base)
b_ratio = marked_conjugate_B(fixed_base)/original_B(fixed_base)
twist = b_ratio/a_ratio
assert twist**2 == a_ratio and twist**3 == b_ratio

original_values = tuple(
    (point[0](fixed_base), point[1](fixed_base), scales, indices)
    for point, scales, indices in original_fourths
)
conjugate_values = tuple(
    (point[0](fixed_base), point[1](fixed_base), scales, indices)
    for point, scales, indices in conjugate_fourths
)
matched_x = []
for original_index, original in enumerate(original_values):
    for conjugate_index, conjugate in enumerate(conjugate_values):
        if conjugate[0] == twist*original[0]:
            matched_x.append((original_index, conjugate_index))

print(
    "Q80FOURTHFIXEDJGF73|prime=73|fixed_base=-27|"
    f"conjugate_normalization_untwist={int(normalization_twist)}|"
    f"twist={int(twist)}|twist_square={int(twist.is_square())}|"
    f"original_markings={original_matches}|"
    f"original_unique_fourths={len(original_fourths)}|"
    f"conjugate_markings={conjugate_matches}|"
    f"conjugate_unique_fourths={len(conjugate_fourths)}",
    flush=True,
)
print(
    "Q80FOURTHFIXEDJGF73|"
    f"original_values={tuple((int(a), int(b)) for a, b, _, _ in original_values)}|"
    f"conjugate_values={tuple((int(a), int(b)) for a, b, _, _ in conjugate_values)}|"
    f"twist_x_matches={tuple(matched_x)}|"
    f"status={'PASS_HORIZONTAL_MATCH' if matched_x else 'PASS_NO_HORIZONTAL_MATCH'}",
    flush=True,
)

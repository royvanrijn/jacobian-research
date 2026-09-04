#!/usr/bin/env sage
"""Exhaustive local audit of the determinant-1236 V4 tower.

The asserted tower has common elliptic quotient

    E = 618f1: Y^2 + X*Y = X^3 - 185*X + 1401.

The known genus-two quotient is the quadratic cover with squareclass

    h = (X-4)/54,

whose Prym is 618e1.  If the other two character spaces are respectively
618a1 x 618b1 and 618c1 x 618d1, a genus-three double cover z^2=b of E and
the complementary cover z^2=h*b must have those two Prym traces, in either
order.

Every genus-three quadratic extension of F_5(E) has a representative in
L(6*O): multiply its squareclass by a square so that

    div(b) = (four branch points) + 2*Q - 6*O.

By default this script enumerates every point of P(L(6*O))(F_5), retains
exactly the classes with branch degree four, and compares both degree-two
character sums.  An optional good odd prime can be passed as the first
argument for diagnostics.  Degree two is insensitive to all quadratic
constant twists.  Thus a zero compatible count cannot be repaired by
changing a twist.

The result is a consistency audit, not an arithmetic exclusion of the
marked curve: it proves that the currently asserted quotient model and
elliptic-factor partition cannot all be used simultaneously at p=5.  The
other two pair partitions of 618a1,...,618d1 are recorded as diagnostics;
they are not substitutes for the Atkin--Lehner eigenspace calculation.  The
upstream quotient/factor identification must be reconciled before cover
lifting.
"""

import json
import sys
from itertools import product


P = 5 if len(sys.argv) == 1 else ZZ(sys.argv[1])
assert P.is_prime() and P not in (2, 3, 103)
FIELD = GF(P)
EXTENSION = GF(P**2, "z")
E_RATIONAL = EllipticCurve(QQ, [1, 0, 0, -185, 1401])
E = EllipticCurve(EXTENSION, [1, 0, 0, -185, 1401])
POLYNOMIAL_RING = PolynomialRing(FIELD, "X")
X = POLYNOMIAL_RING.gen()
WEIERSTRASS_CUBIC = X**3 - 185*X + 1401
X_MAP_DISCRIMINANT = X**2 + 4*WEIERSTRASS_CUBIC
B_POLYNOMIAL = 1944*X**6 + 441*X**4 - 90*X**2 + 9
POLE_WEIGHTS = [0, 2, 3, 4, 5, 6]

assert E_RATIONAL.has_good_reduction(P)
assert B_POLYNOMIAL.degree() == 6
assert B_POLYNOMIAL.gcd(B_POLYNOMIAL.derivative()) == 1


def projective_vectors(field, vector_dimension):
    for pivot in range(vector_dimension):
        for tail in product(field, repeat=vector_dimension-pivot-1):
            yield vector(field, [0]*pivot + [1] + list(tail))


def quadratic_character(value):
    if not value:
        return 0
    return 1 if value.is_square() else -1


def power_trace_two(curve, prime):
    return ZZ(curve.ap(prime))**2 - 2*prime


def pair_target(first, second):
    return -power_trace_two(EllipticCurve(first), P) - power_trace_two(
        EllipticCurve(second), P
    )


def pole_data(coefficients):
    nonzero_indices = [
        index for index, value in enumerate(coefficients) if value
    ]
    pole_degree = max(POLE_WEIGHTS[index] for index in nonzero_indices)
    leading_index = POLE_WEIGHTS.index(pole_degree)
    return pole_degree, coefficients[leading_index]


def coefficient_polynomials(coefficients):
    a_value = (
        coefficients[0]
        + coefficients[1]*X
        + coefficients[3]*X**2
        + coefficients[5]*X**3
    )
    b_value = coefficients[2] + coefficients[4]*X
    return a_value, b_value


def norm_from_polynomials(a_value, b_value):
    return (
        a_value**2
        - X*a_value*b_value
        - WEIERSTRASS_CUBIC*b_value**2
    )


def norm_polynomial(coefficients):
    return norm_from_polynomials(*coefficient_polynomials(coefficients))


def branch_degree(coefficients):
    # A common factor g(X) of A and B vanishes at both points above a
    # non-ramified root of g.  Its square in Norm(b) would erase both odd
    # valuations if one inspected only the norm parity.  Separate it first.
    a_value, b_value = coefficient_polynomials(coefficients)
    common = a_value.gcd(b_value).monic()
    primitive_a = a_value // common
    primitive_b = b_value // common
    primitive_norm = norm_from_polynomials(primitive_a, primitive_b)
    primitive_degree = sum(
        factor.degree()
        for factor, multiplicity in primitive_norm.factor()
        if multiplicity % 2
    )
    common_degree = sum(
        2*factor.degree()
        for factor, multiplicity in common.factor()
        if multiplicity % 2 and not X_MAP_DISCRIMINANT.mod(factor).is_zero()
    )
    pole_degree, _ = pole_data(coefficients)
    return primitive_degree + common_degree + pole_degree % 2


def character_sums(coefficients):
    lifted = [EXTENSION(value) for value in coefficients]
    pole_degree, leading = pole_data(lifted)
    if pole_degree == 0:
        sum_b = quadratic_character(leading)
    elif pole_degree % 2:
        sum_b = 0
    else:
        sum_b = quadratic_character(leading)

    product_pole_degree = pole_degree + 2
    sum_product = (
        0
        if product_pole_degree % 2
        else quadratic_character(leading / EXTENSION(54))
    )
    for point in E:
        if point.is_zero():
            continue
        x_value, y_value = point[0], point[1]
        value = (
            lifted[0]
            + lifted[1]*x_value
            + lifted[2]*y_value
            + lifted[3]*x_value**2
            + lifted[4]*x_value*y_value
            + lifted[5]*x_value**3
        )
        sum_b += quadratic_character(value)
        sum_product += quadratic_character(
            value*(x_value-4)/EXTENSION(54)
        )
    return int(sum_b), int(sum_product)


target_ab = int(pair_target("618a1", "618b1"))
target_cd = int(pair_target("618c1", "618d1"))
target_e = int(-power_trace_two(EllipticCurve("618e1"), P))
pair_partitions = {
    "ab|cd": (target_ab, target_cd),
    "ac|bd": (
        int(pair_target("618a1", "618c1")),
        int(pair_target("618b1", "618d1")),
    ),
    "ad|bc": (
        int(pair_target("618a1", "618d1")),
        int(pair_target("618b1", "618c1")),
    ),
}

assert branch_degree(vector(FIELD, [-4, 1, 0, 0, 0, 0])) == 2

# Unit-test the character convention on the already explicit B -> E cover.
known_cover_sum = 0
for point in E:
    if point.is_zero():
        known_cover_sum += quadratic_character(EXTENSION(1)/EXTENSION(54))
    else:
        known_cover_sum += quadratic_character(
            (point[0]-4)/EXTENSION(54)
        )
assert known_cover_sum == target_e

projective_count = 0
branch_four_count = 0
unfiltered_compatible_by_partition = {name: [] for name in pair_partitions}
signature_counts = {}
compatible_by_partition = {name: [] for name in pair_partitions}
for coefficients in projective_vectors(FIELD, 6):
    projective_count += 1
    signature = character_sums(coefficients)
    for name, targets in pair_partitions.items():
        if signature in (targets, tuple(reversed(targets))):
            unfiltered_compatible_by_partition[name].append(
                [int(value) for value in coefficients]
            )
    if branch_degree(coefficients) != 4:
        continue
    branch_four_count += 1
    signature_counts[signature] = signature_counts.get(signature, 0) + 1
    for name, targets in pair_partitions.items():
        if signature in (targets, tuple(reversed(targets))):
            compatible_by_partition[name].append(
                [int(value) for value in coefficients]
            )

assert projective_count == (P**6-1)//(P-1)
claimed_compatible = compatible_by_partition["ab|cd"]
if P == 5:
    assert claimed_compatible == []

if claimed_compatible:
    status = "LOCALLY_COMPATIBLE"
    conclusion = (
        "At p=%s the displayed B/E squareclass and asserted ab|cd V4 "
        "factor partition have local genus-three classes. This local "
        "diagnostic neither constructs the characteristic-zero cover nor "
        "decides a rational lift."
    ) % P
else:
    status = "EXACT_UPSTREAM_INCONSISTENCY"
    conclusion = (
        "The displayed B/E squareclass and the asserted ab|cd V4 factor "
        "partition admit no genus-three quadratic cover over F_%s in "
        "either order. The other pair partitions are diagnostic only. "
        "This is not an arithmetic exclusion of C_1236; the upstream "
        "quotient/factor identification must be reconciled before lifting."
    ) % P

payload = {
    "schema": "elkies-k3.det1236-v4-local-consistency.v1",
    "status": status,
    "prime": P,
    "field_extension_degree": 2,
    "elliptic_base": "618f1",
    "known_cover_squareclass": "(X-4)/54",
    "known_cover_prym": "618e1",
    "known_cover_character_sum": int(known_cover_sum),
    "known_cover_target": target_e,
    "claimed_complementary_pryms": [
        ["618a1", "618b1"],
        ["618c1", "618d1"],
    ],
    "degree_two_targets": {"ab": target_ab, "cd": target_cd},
    "all_pair_partitions": {
        name: {
            "left_target": targets[0],
            "right_target": targets[1],
            "compatible_classes": compatible_by_partition[name],
        }
        for name, targets in pair_partitions.items()
    },
    "unfiltered_compatible_class_counts": {
        name: len(values)
        for name, values in unfiltered_compatible_by_partition.items()
    },
    "projective_classes_enumerated": int(projective_count),
    "branch_degree_four_classes": int(branch_four_count),
    "distinct_signature_pairs": int(len(signature_counts)),
    "compatible_classes": claimed_compatible,
    "quadratic_twist_boundary": (
        "degree-two character sums are unchanged by quadratic constant twists"
    ),
    "conclusion": conclusion,
}
print(json.dumps(payload, indent=2, sort_keys=True, default=int))

#!/usr/bin/env python3
"""Exact replay of the supporting-face and branch-partition calculations."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd


# Keys encode (rho, H) as (rho_numerator, H_numerator, common_denominator).
EXPECTED_QUINTIC = {
    (5, 5, 3): (-5, -2, 1),
    (2, 2, 1): (-3, -2, -1, 0, 1),
    (3, 5, 2): (-5, -3, -1, 1),
    (3, 6, 2): (-4, -2, 0, 2),
    (2, 3, 1): (-2, -1, 0, 1),
    (3, 3, 1): (-1, 0, 1),
    (4, 10, 3): (-5, -2, 1),
    (3, 7, 2): (-3, -1, 1),
    (4, 11, 3): (-4, -1, 2),
    (5, 15, 4): (-5, -1, 3),
    (4, 12, 3): (-3, 0, 3),
    (3, 8, 2): (-2, 0, 2),
    (2, 4, 1): (-1, 0, 1, 2),
    (3, 9, 2): (-1, 1, 3),
}

EXPECTED_CUBIC = {
    (2, 2, 1): (-1, 0, 1),
    (3, 3, 2): (-3, -1, 1),
}

# Primitive signatures in increasing primitive-character order.
EXPECTED_SIGNATURES = {
    7: (1, 1, 2, 0, 1, 1),
    9: (1, 1, 2, 0, 1, 1),
    10: (1, 2, 0, 1),
    11: (1, 1, 2, 0, 0, 2, 2, 0, 1, 1),
    15: (1, 1, 0, 2, 0, 2, 1, 1),
}


def enumerate_faces(degree):
    nodes = [
        (angle, height)
        for angle in range(-degree, degree)
        for height in range(
            max(0, -2 * angle), degree + (angle < 0) - angle
        )
    ]
    negative = [point for point in nodes if point[0] < 0]
    positive = [point for point in nodes if point[0] > 0]
    top = [(i, j) for i, j in negative if i + j == degree]
    faces = {}
    for (left, high), (right, low) in product(negative, positive):
        rho = Fraction(high - low, right - left)
        intercept = high + rho * left
        if rho <= 1 or intercept >= degree:
            continue
        if all(j + rho * i > intercept for i, j in top):
            continue
        angles = tuple(i for i, j in nodes if j + rho * i == intercept)
        if len(angles) >= 3:
            faces[rho, intercept] = angles
    return faces


def check_census(faces, expected):
    encoded = {
        (Fraction(r, q), Fraction(h, q)): angles
        for (r, h, q), angles in expected.items()
    }
    assert faces == encoded, (faces, encoded)


def occupied_supports(allowed):
    for size in range(3, len(allowed) + 1):
        for occupied in combinations(allowed, size):
            if occupied[0] < 0 < occupied[-1]:
                yield occupied


def root_partitions(degree, minimum=1):
    if degree == 0:
        yield ()
        return
    for multiplicity in range(minimum, degree + 1):
        for rest in root_partitions(degree - multiplicity, multiplicity):
            yield (multiplicity, *rest)


def signature(order, exponents):
    assert sum(exponents) == 2 * order
    assert gcd(order, *exponents) == 1
    dimensions = {}
    for index in range(1, order):
        residue_sum = sum(index * exponent % order for exponent in exponents)
        assert residue_sum % order == 0
        dimensions[index] = residue_sum // order - 1
    assert all(dimension >= 0 for dimension in dimensions.values())
    assert dimensions[1] == 1
    genus = 1 - order + sum(order - gcd(order, a) for a in exponents) // 2
    assert genus == sum(dimensions.values())
    return dimensions


def check_collision(intercept, left, right, partition):
    if intercept == Fraction(3, 2):
        assert (left, right, partition) == (-3, 1, (2,))
        return "cubic repeated-root shear"
    if intercept == Fraction(5, 3):
        assert (left, right, partition) == (-5, 1, (2,))
        return "quintic repeated-quadratic shear"
    if intercept == Fraction(5, 2):
        assert (left, right, partition) == (-5, 1, (3,))
        return "quintic triple-cubic shear"
    assert intercept == 3 and (left, right, partition) == (-2, 1, (3,))
    return "height-three binomial coefficient"


def check_factor_routes(order, exponents, partition):
    dimensions = signature(order, exponents)
    factors = defaultdict(dict)
    for index, dimension in dimensions.items():
        factors[order // gcd(index, order)][index] = dimension

    routes = []
    for conductor, factor in sorted(factors.items()):
        dimension = sum(factor.values())
        if dimension == 0:
            continue
        if conductor == order and partition == (2,) and order in EXPECTED_SIGNATURES:
            assert {value + factor[order - index] for index, value in factor.items()} == {1}
            routes.append("rank-one CM")
            continue
        if dimension <= 2:
            routes.append("elliptic or surface")
            continue
        assert conductor == order, (order, exponents, factors)
        if order in (3, 4):
            assert dimension == 3
            routes.append("imaginary-quadratic threefold")
            continue

        hodge_ranks = {value + factor[order - index] for index, value in factor.items()}
        assert hodge_ranks == {2}
        assert all(value in (0, 1, 2) for value in factor.values())
        assert sum(value == 1 for value in factor.values()) <= 4
        if order in EXPECTED_SIGNATURES:
            assert tuple(factor.values()) == EXPECTED_SIGNATURES[order]
        else:
            assert order == 5 and sorted(e for e in exponents if e) == [2, 2, 2, 4]
        routes.append("rank-two CM theorem")
    return " + ".join(routes)


def audit_occupied_faces(faces):
    counts = Counter()
    for (rho, intercept), allowed in faces.items():
        for occupied in occupied_supports(allowed):
            if rho == intercept == 2:
                counts["quadratic Laurent height"] += 1
                continue

            if intercept.denominator == 1:
                angular_step = gcd(*allowed)
                left, right = occupied[0] // angular_step, occupied[-1] // angular_step
                step = 1
            else:
                left, right = occupied[0], occupied[-1]
                step = intercept.denominator
                assert step == rho.denominator
            order = intercept.numerator
            assert (right - left) % step == 0
            if intercept.denominator != 1:
                assert all((angle - left) % step == 0 for angle in occupied)
            assert left < 0 < right
            assert left >= -order and right <= order

            polynomial_degree = (right - left) // step
            for partition in root_partitions(polynomial_degree):
                exponents = (order + left, *(step * r for r in partition), order - right)
                assert sum(exponents) == 2 * order
                if max(exponents) >= order:
                    counts[check_collision(intercept, left, right, partition)] += 1
                else:
                    counts[check_factor_routes(order, exponents, partition)] += 1
    return counts


def main():
    if not __debug__:
        raise SystemExit("Run the replay without -O so its assertions are enabled.")
    for degree, expected in (
        (1, {}), (2, {}), (3, EXPECTED_CUBIC), (5, EXPECTED_QUINTIC)
    ):
        faces = enumerate_faces(degree)
        check_census(faces, expected)
        routes = audit_occupied_faces(faces)
        print(f"Degree {degree}: {len(faces)} lines, {routes.total()} support/partition cases.")
        for route, count in sorted(routes.items()):
            print(f"  {count}: {route}")

    # Deleting a single listed face must invalidate the complete census.
    damaged = dict(EXPECTED_QUINTIC)
    del damaged[(3, 5, 2)]
    try:
        check_census(enumerate_faces(5), damaged)
    except AssertionError:
        print("Negative control passed: removing the height-5/2 line was rejected.")
    else:
        raise AssertionError("The damaged census passed.")


if __name__ == "__main__":
    main()

"""Exact residue-class census and normalization and gluing data for SN1."""

from itertools import product
from math import prod

from sage.all import GF, ZZ, HyperellipticCurve, PolynomialRing, pari, power_mod


field = GF(5)
ring = PolynomialRing(field, "t")
t = ring.gen()
torsion_ring = PolynomialRing(GF(3), "X")
X = torsion_ring.gen()
last_case = (0, 2, 0, 0, 0, 1)

# Split orbit degrees, normalization order, upper bound b0.
expected_cases = {
    (0, 0, 1, 4, 3, 2): ((1,), 4536, 4),
    last_case: ((2, 2), 21, 1),
    (0, 4, 2, 2, 0, 2): ((2,), 756, 3),
    (1, 0, 0, 4, 2, 0): ((1,), 4536, 4),
    (2, 2, 1, 1, 4, 1): ((1, 1), 1512, 3),
    (2, 2, 1, 2, 1, 4): ((1, 2), 252, 2),
}


def fiber_polynomial(coefficients):
    polynomial = t**7+ring(coefficients)
    return polynomial**2+1


def finite_rank_budget(normalization, repeated_roots):
    genus = (normalization.degree()-2)//2
    if genus > 0:
        frobenius = pari(normalization).hyperellcharpoly()
        order = ZZ(frobenius(1))
        char_mod_three = torsion_ring([int(c) % 3 for c in reversed(frobenius.Vec())])
        multiplicity = 0
        while char_mod_three % (X-1) == 0:
            multiplicity += 1
            char_mod_three //= X-1
        bound = min(order.valuation(3), multiplicity)
    else:
        order = ZZ(1)
        bound = 0

    node_types = []
    for polynomial, exponent in repeated_roots.factor():
        assert exponent == 1
        degree = polynomial.degree()
        character = power_mod(normalization, (5**degree-1)//2, polynomial)
        assert character in (-1, 1)
        node_types.append((degree, 1 if character == 1 else -1))
    # Prime-to-five specialization of the infinity difference into Pic^0.
    assert order*prod(5**d-sign for d, sign in node_types) % 7 == 0
    torus = sum((5**d-sign) % 3 == 0 for d, sign in node_types)
    components = sum(sign == 1 for _, sign in node_types)
    return order, bound, node_types, bound+torus+components


def verify_census():
    candidates = {}
    total = 0
    wilder = 0
    for coefficients in product(range(5), repeat=6):
        reflection = tuple((-1)**(i+1)*c % 5 for i, c in enumerate(coefficients))
        if coefficients > reflection:
            continue
        total += 1
        fiber = fiber_polynomial(coefficients)
        repeated_roots = fiber.gcd(fiber.derivative())
        if not repeated_roots.is_squarefree():
            wilder += 1
            continue
        normalization, remainder = fiber.quo_rem(repeated_roots**2)
        assert remainder == 0 and normalization.is_squarefree()
        assert normalization.gcd(repeated_roots) == 1
        order, bound, types, total_bound = finite_rank_budget(normalization, repeated_roots)
        if total_bound >= 5:
            assert total_bound == 5
            assert all(sign == 1 for _, sign in types)
            candidates[coefficients] = (tuple(sorted(d for d, _ in types)), order, bound)
    assert total == 7875 and wilder == 616
    assert candidates == expected_cases


def jacobian_and_point(coefficients, x, y):
    fiber = fiber_polynomial(coefficients)
    repeated_roots = fiber.gcd(fiber.derivative())
    normalization = fiber//repeated_roots**2
    curve = HyperellipticCurve(normalization)
    jacobian = curve.jacobian()(field)
    plus = next(point for point in curve.points_at_infinity() if point[1] == 1)
    difference = jacobian(curve(x, y), plus)
    return jacobian, difference


def verify_order_witnesses():
    cases = [
        ((0, 0, 1, 4, 3, 2), (0, 1), 81),
        ((0, 4, 2, 2, 0, 2), (0, 1), 9),
        ((1, 0, 0, 4, 2, 0), (2, 2), 27),
        ((2, 2, 1, 1, 4, 1), (1, 2), 27),
        ((2, 2, 1, 2, 1, 4), (2, 1), 9),
    ]
    for coefficients, coordinates, exact_order in cases:
        jacobian, difference = jacobian_and_point(coefficients, *coordinates)
        orbit_degrees, order, bound = expected_cases[coefficients]
        prime_to_three = order//3**ZZ(order).valuation(3)
        witness = prime_to_three*difference
        assert exact_order*witness == jacobian(0)
        assert (exact_order//3)*witness != jacobian(0)
        remaining_bound = ZZ(order).valuation(3)-ZZ(exact_order).valuation(3)+1
        assert remaining_bound < bound
        assert remaining_bound+len(orbit_degrees)+sum(d % 2 == 0 for d in orbit_degrees) <= 4


def verify_torus_gluing():
    first = t*t+2*t+4
    second = t*t+3*t+4
    unramified_part = (t**3+t*t+2)*(t**3+4*t*t+3)
    assert fiber_polynomial(last_case) == (first*second)**2*unramified_part
    assert first.is_irreducible() and second.is_irreducible()
    assert unramified_part % first == 1 and unramified_part % second == 1
    jacobian, difference = jacobian_and_point(last_case, 0, 1)
    assert 7*difference != jacobian(0) and 21*difference == jacobian(0)

    A = (
        4*t**21+3*t**20+4*t**19+3*t**18+t**16+4*t**15+t**14+4*t**13
        +2*t**12+t**11+3*t**10+t**8+2*t**7+2*t**6+4*t**4+t**3+2*t+1
    )
    B = (
        4*t**18+3*t**17+t**16+2*t**15+4*t**14+4*t**13+t**12+3*t**11
        +2*t**10+t**8+4*t**7+2*t**4+3*t*t+3*t+4
    )
    assert A*A-B*B*unramified_part == 2*t**21
    assert (A+B)(0) == 0 and (A-B)(0) != 0
    assert A.leading_coefficient() == B.leading_coefficient() != 0
    quotient = ring.quotient(first, "a")
    a = quotient.gen()
    ratio = quotient(A+B)/quotient(A-B)
    assert ratio == a+2 and ratio**8 == 3*a != 1


if __name__ == "__main__":
    verify_census()
    verify_order_witnesses()
    verify_torus_gluing()
    print("SN1 complete census, five point orders, and final gluing obstruction verified")

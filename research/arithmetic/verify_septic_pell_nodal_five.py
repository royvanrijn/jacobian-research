"""Exact residue-class census and normalization and gluing data for SN1–SN6."""

from itertools import product
from math import prod

from sage.all import GF, ZZ, HyperellipticCurve, PolynomialRing, matrix, pari, power_mod


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


def finite_rank_budget(normalization, repeated_roots, *, check_infinity=True):
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
    if check_infinity:
        assert order*prod(5**d-sign for d, sign in node_types) % 7 == 0
    torus = sum((5**d-sign) % 3 == 0 for d, sign in node_types)
    components = sum(sign == 1 for _, sign in node_types)
    return order, bound, node_types, bound+torus+components


def verify_census():
    candidates = {}
    total = 0
    wilder = 0
    deeper_candidates = {}
    for coefficients in product(range(5), repeat=6):
        reflection = tuple((-1)**(i+1)*c % 5 for i, c in enumerate(coefficients))
        if coefficients > reflection:
            continue
        total += 1
        fiber = fiber_polynomial(coefficients)
        repeated_roots = fiber.gcd(fiber.derivative())
        if not repeated_roots.is_squarefree():
            wilder += 1
            normalization = fiber
            pencils = ring(1)
            largest = 0
            for polynomial, exponent in fiber.factor():
                largest = max(largest, exponent)
                normalization //= polynomial**(2*(exponent//2))
                if exponent == 2:
                    pencils *= polynomial
            genus = (normalization.degree()-2)//2
            assert normalization.is_squarefree()
            order, _, _, initial = finite_rank_budget(
                normalization, pencils, check_infinity=False
            )
            budget = initial+6-genus-pencils.degree()
            if budget >= 5:
                deeper_candidates[coefficients] = largest, order, budget
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
    assert deeper_candidates == {
        (0, 1, 0, 1, 1, 4): (4, 1134, 6),
        (1, 2, 4, 2, 3, 0): (4, 126, 5),
        (2, 0, 0, 0, 0, 0): (7, 126, 5),
        (2, 0, 0, 1, 0, 2): (3, 21, 5),
        (2, 0, 1, 0, 0, 4): (6, 42, 5),
        (2, 0, 3, 3, 2, 2): (3, 7, 5),
        (2, 4, 3, 4, 4, 4): (4, 126, 5),
    }


def jacobian_and_point(coefficients, x, y):
    fiber = fiber_polynomial(coefficients)
    normalization = fiber
    for polynomial, exponent in fiber.factor():
        normalization //= polynomial**(2*(exponent//2))
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


def verify_deeper_normalizations():
    for coefficients, coordinates, primary_order, exact_order in [
        ((0, 1, 0, 1, 1, 4), (0, 1), 14, 81),
        ((2, 0, 0, 0, 0, 0), (2, 2), 14, 9),
        ((2, 4, 3, 4, 4, 4), (1, 2), 14, 9),
    ]:
        jacobian, difference = jacobian_and_point(coefficients, *coordinates)
        witness = primary_order*difference
        assert exact_order*witness == jacobian(0)
        assert (exact_order//3)*witness != jacobian(0)

    for coefficients, t_coordinate in [
        ((1, 2, 4, 2, 3, 0), 1),
        ((2, 0, 1, 0, 0, 4), 0),
    ]:
        jacobian, plus = jacobian_and_point(coefficients, t_coordinate, 2)
        _, minus = jacobian_and_point(coefficients, t_coordinate, -2)
        branch_difference = 14*(plus-minus)
        assert branch_difference != jacobian(0) and 3*branch_difference == jacobian(0)

    jacobian, point = jacobian_and_point((1, 2, 4, 2, 3, 0), 1, 2)
    first = 14*point
    second = 14*jacobian(t*t+4*t+1, t+4)
    assert 3*first == 3*second == jacobian(0)
    assert first != jacobian(0) and second not in (jacobian(0), first, -first)


def verify_deeper_gluing():
    core = t*(t**5+3*t**4+2*t*t+2*t+4)
    h = t*t+t+2
    A = t**6+4*t**4+3*t*t+1
    B = t**3+t*t+3*t+3
    assert h.is_irreducible() and (A+B*(t+4)) % h == 0
    assert (A-B*(t+4)) % h != 0 and A*A-core*B*B == 2*h**3
    assert A.leading_coefficient() == B.leading_coefficient() != 0
    assert core(2) == 2 and A(2) == B(2) == 1
    local = ring.quotient(t*t-2, "b")
    b = local.gen()
    ratio = (1+b)/(1-b)
    assert ratio**6 == 1 and ratio**2 == 2+2*b != 1

    core = t**4+3*t**3+3*t+3
    A = 3*t**7+3*t**6+3*t**5+4*t**4+4*t*t+4*t+1
    B = t**4+t**3+t*t
    assert A*A-core*B*B == 4*(t-2)**14
    assert core(2) == 4 and A(2)+2*B(2) != 0
    h = t*t+4*t+1
    assert h.is_irreducible()
    local = ring.quotient(h, "a")
    a = local.gen()
    b = 2*a+2
    assert b*b == local(core) and b != 0
    ratio = (local(A)+local(B)*b)/(local(A)-local(B)*b)
    assert ratio == 2*a and ratio**8 == a+4 != 1


def verify_split_local_family():
    p = 5
    for a in (1, 2, 3):
        modulus = 3**a
        k = modulus*(3*pow(modulus, -1, 5) % 5)
        N = 4
        centers = [0, 1, 2, 3, 4, 5**k, 1+5**k]
        derivatives = []
        for c in centers:
            differences = [ZZ(c-d).valuation(p) for d in centers if d != c]
            assert set(differences) <= {0, k}
            value = sum(differences)
            assert k*N > 2*value
            derivatives.append(value)
        assert derivatives == [k, k, 0, 0, 0, k, k]
        lengths = [k]*4+[2*k*N]*3+[2*k*(N-2)]*4
        assert len(lengths)-6+1 == 6
        assert all(length > 0 and length % modulus == 0 for length in lengths)

    residue = GF(11)
    ring_11 = PolynomialRing(residue, "t")
    z = ring_11.gen()
    assert residue(5).multiplicative_order() == 5
    d = residue(5)**3
    polynomial = prod(z-c for c in [0, 1, 2, 3, 4, d, d+1])
    fiber = polynomial**2+residue(5)**24
    assert fiber.degree() == 14 and fiber.is_squarefree()
    frobenius = pari(fiber).hyperellcharpoly()
    assert ZZ(frobenius(1)) == 2643158 and 2643158 % 3 != 0


def verify_split_residue_section_profiles():
    polynomials = [
        leading*t**degree+ring(coefficients)
        for degree in (2, 5)
        for coefficients in product(range(5), repeat=degree)
        for leading in range(1, 5)
    ]
    expected = {
        (0, 0, 1, 2, 4): (4, 2),
        (0, 0, 1, 5, 1): (4, 2),
        (0, 1, 1, 1, 4): (2, 1),
        (0, 1, 1, 3, 2): (5, 2),
        (1, 1, 1, 1, 3): (2, 1),
    }
    nonempty = {}
    total = 0
    for multiplicities in product(range(8), repeat=5):
        if sum(multiplicities) != 7 or sum(e > 0 for e in multiplicities) < 3:
            continue
        orbit = [
            tuple(multiplicities[(a*i+b) % 5] for i in range(5))
            for a in range(1, 5) for b in range(5)
        ]
        if multiplicities != min(orbit):
            continue
        total += 1
        factors = [(t-r, e) for r, e in enumerate(multiplicities) if e]
        Q = prod(h**e for h, e in factors)
        rows = []
        for Z in polynomials:
            T = 2*Q+Z
            # Coprime factors would each have cube degree, whereas deg T = 7.
            if Z.gcd(T) == 1:
                continue
            if all(e % 3 == 0 for _, e in (Z*T).factor()):
                if multiplicities in ((0, 0, 1, 2, 4), (0, 0, 1, 5, 1), (0, 1, 1, 3, 2)):
                    assert Z.degree() == 5
                rows.append([(T.valuation(h)-e) % 3 for h, e in factors])
        if rows:
            nonempty[multiplicities] = len(rows), matrix(GF(3), rows).rank()
    assert total == 16 and nonempty == expected

    W = t**5-t
    Q = t*t*W
    valid = {}
    for Z in polynomials:
        T = 2*Q+Z
        if Z.gcd(T) != 1 and all(e % 3 == 0 for _, e in (Z*T).factor()):
            valid[Z] = tuple((T.valuation(t-r)-Q.valuation(t-r)) % 3 for r in range(5))
    assert valid == {
        t**5+4*t**3: (0, 1, 2, 2, 1),
        4*t**5+4*t**3: (0, 2, 1, 1, 2),
    }


if __name__ == "__main__":
    verify_census()
    verify_order_witnesses()
    verify_torus_gluing()
    verify_deeper_normalizations()
    verify_deeper_gluing()
    verify_split_local_family()
    verify_split_residue_section_profiles()
    print("SN1–SN6 exact censuses, gluing obstructions, and split local family verified")

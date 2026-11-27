"""Verify norm, prime-ratio, reduction and Selmer fixtures for TP1–TP7."""

from sage.all import (
    GF, QQ, ZZ, EllipticCurve, PolynomialRing, inverse_mod, matrix, pari, power_mod, valuation,
)


def check_seed(level, p, q):
    m = 3**level
    modulus = 3**(2*level + 3)
    residue = 1 - 3**(2*level + 2)
    assert ZZ(p).is_prime() and ZZ(q).is_prime()
    assert p % 8 == 1 and q % 8 == 5
    assert (p - residue) % modulus == (q - residue) % modulus == 0
    assert q % p == 4 and power_mod(p*inverse_mod(4, q), (q - 1)//m, q) == 1
    A = (p - q)//4
    B = -(p + q)//2
    ring = PolynomialRing(QQ, "u")
    u = ring.gen()
    radicand = A**2*u**(2*m) + B*u**m + 1
    assert B**2 - 4*A**2 == p*q
    assert radicand.is_squarefree()
    assert (A*u**m + 1)**2 - radicand == p*u**m
    assert (A*u**m - 1)**2 - radicand == q*u**m

    at_two = (-inverse_mod(B, 8)) % 8
    assert radicand(at_two) % 8 == 1
    assert radicand(1) % modulus == 3**(2*level + 2)
    assert valuation(radicand(1), p) == 1
    at_q = ZZ(GF(q)(QQ(1)/A).nth_root(m))
    assert valuation(radicand(at_q), q) == 1

    numerator_frame = [QQ(4)/p, QQ(q)/4]
    denominator_frame = [QQ(4*A*A)/p, QQ(4*A*A)/q]
    for frame in (numerator_frame, denominator_frame):
        valuations = matrix(ZZ, [[valuation(entry, fixed_prime) for entry in frame]
                                for fixed_prime in (p, q)])
        assert abs(valuations.det()) == 1

    check_dyadic_charts(A, B)

    # The full moving class subgroup has these two Smith factors.
    relation_lattice = matrix(ZZ, [[4*m, 0], [0, 4*m], [2*m, -2*m]])
    assert list(relation_lattice.smith_form()[0].diagonal()) == [2*m, 4*m]


def check_dyadic_charts(A, B):
    d = valuation(A, 2)
    A0 = QQ(A)/2**d
    p = 2*A - B
    polynomial_ring = PolynomialRing(QQ, "x")
    rational_field = polynomial_ring.fraction_field()
    x = rational_field.gen()

    def cyclic_cover(w):
        return w*(2 - w)/(p - 2*A*w)

    assert cyclic_cover(2*x)/4 == x*(1 - x)/(p - 4*A*x)
    assert 1/(2**(2 + 2*d)*cyclic_cover(p/(2*A*x))) == A0**2*x*(x - 1)/(4*A*x - p)
    other_node_coordinate = 4*A/(p*x)
    assert cyclic_cover(2/x) == other_node_coordinate/x*(x - 1)/(A*(1 - other_node_coordinate))


def check_reciprocal_seed(level, a, p, q):
    m = 3**level
    modulus = 3**(2*level + 3)
    residue = 1 - 3**(2*level + 2)
    assert a % (2*modulus) == 0 and p - q == 4*a**m
    assert ZZ(p).is_prime() and ZZ(q).is_prime()
    assert (p - residue) % modulus == (q - residue) % modulus == 0

    A = a**m
    B = -(p + q)//2
    for fixed_prime, sign in ((p, -1), (q, 1)):
        u_mod = sign*inverse_mod(a, fixed_prime)
        z_mod = power_mod(u_mod, m, fixed_prime**2)
        assert valuation(A*A*z_mod*z_mod + B*z_mod + 1, fixed_prime) == 1
    at_two = [u for u in range(1, 32, 2)
              if (A*A*pow(u, 2*m, 32) + B*pow(u, m, 32) + 1) % 32 == 4]
    assert len(at_two) == 1
    check_dyadic_charts(A, B)

    c = a**(m//3)
    h = 3*c
    polynomial_ring = PolynomialRing(QQ, "w")
    rational_field = polynomial_ring.fraction_field()
    w = rational_field.gen()
    radicand = w**6 + QQ(B)/A*w**3 + 1
    curves = []
    for sign, fixed_prime, complementary_prime in ((1, p, q), (-1, q, p)):
        b = sign*fixed_prime
        X = -fixed_prime*w/(c*(w + sign)**2)
        Y_factor = -fixed_prime/(w + sign)**3
        linear_term = h*X + b
        assert Y_factor**2*radicand == X**3 + linear_term**2
        assert linear_term == -sign*Y_factor*(w**3 + sign)
        assert Y_factor**2/X**3 == -QQ(A)/(fixed_prime*w**3)

        curve = EllipticCurve([0, h*h, 0, 2*h*b, b*b])
        assert curve.discriminant() == -432*fixed_prime**3*complementary_prime
        assert all(curve.c4() % prime != 0 for prime in (p, q))
        assert (b - h)**2 % 9 == 1
        curves.append(curve)

    partner = EllipticCurve([0, -3*h*h, 0, -18*h*q, -27*q*q])
    assert partner.is_isomorphic(curves[1].quadratic_twist(-3))
    if level == 1:
        assert all(curve.conductor() == 36*p*q for curve in curves)
        assert partner.j_invariant() in [isogeny.codomain().j_invariant()
                                         for isogeny in curves[0].isogenies_prime_degree(3)]
        check_dyadic_twists(p, q, curves)


def check_genuine_ratios():
    a, p, q = 486, 459167617, 2593
    A, B = a**3, -(p + q)//2
    Q = ZZ(1000000374033487612847617)
    assert Q.is_prime(proof=True)

    def cubic_residue(element, prime):
        return power_mod(element.numerator()*inverse_mod(element.denominator(), prime),
                         (prime - 1)//3, prime)

    assert cubic_residue(QQ(4*A*A)/p, Q) == 1
    assert cubic_residue(QQ(4*A*A)/q, Q) != 1
    for dyadic, P in ((4, ZZ(2058091146572215809781)), (20, ZZ(2057901351959280331333))):
        assert P.is_prime(proof=True) and P != Q
        assert P % 243 == Q % 243 == 1
        assert cubic_residue(QQ(4)/p, P) != 1
        assert cubic_residue(QQ(q)/4, P) == 1
        value = A*A*P**6 + B*P**3*Q**3 + Q**6
        assert value < 0 and value % 32 == dyadic
        assert value % 243 == 81
        assert value.valuation(p) == value.valuation(q) == 1


def check_dyadic_twists(p, q, curves):
    for sign, split, auxiliary_prime in ((-1, True, 23), (-1, False, 11),
                                         (1, True, 73), (1, False, 13)):
        D = ZZ(sign*p*q*auxiliary_prime)
        assert D.is_squarefree() and D % 3 == 1
        assert D % 8 == (1 if split else 5)
        for index, base in enumerate(curves):
            E = base.quadratic_twist(D)
            partner = curves[1 - index].quadratic_twist(-3*D)
            assert E.conductor() == partner.conductor() == 36*D*D
            assert (E.tamagawa_number(2), partner.tamagawa_number(2)) == (
                (3, 1) if split else (1, 3)
            )
            assert E.tamagawa_number(3) == partner.tamagawa_number(3) == 2
            assert E.root_number() == partner.root_number() == sign


def ray_dimension(discriminant, modulus, kill_two):
    variable = pari("x")
    polynomial = variable**2 - variable + (1 - discriminant)//4
    group = polynomial.bnfinit(1)
    assert group.bnfcertify() == 1
    ray = group.bnrinit(modulus, 1)
    invariants = tuple(map(int, pari("r->r.cyc")(ray)))
    dimension = sum(factor % 3 == 0 for factor in invariants)
    if kill_two:
        primes = group.idealprimedec(2)
        assert len(primes) == 2
        coordinates = tuple(map(int, ray.bnrisprincipal(primes[0], 0)))
        dimension -= any(factor % 3 == 0 and coord % 3 != 0
                         for factor, coord in zip(invariants, coordinates))
    return dimension


def check_ray_and_inertia_controls():
    for D, expected in ((-23, (1, 0)), (-11, (0, 1)), (73, (2, 0)), (13, (0, 0))):
        if D % 8 == 1:
            dimensions = ray_dimension(-3*D, 6, False), ray_dimension(D, 1, True)
        else:
            dimensions = ray_dimension(-3*D, 3, True), ray_dimension(D, 2, False)
        assert dimensions == expected
        assert dimensions[0] - dimensions[1] == (D > 0) + (1 if D % 8 == 1 else -1)

    residue_ring = PolynomialRing(GF(3), "z")
    z = residue_ring.gen()
    for c0 in (1, 2):
        assert (z**3 - z - c0).is_irreducible()
    rotation = matrix(QQ, [[0, -1], [1, -1]])
    assert rotation**3 == 1 and rotation.minimal_polynomial().is_irreducible()


if __name__ == "__main__":
    check_seed(1, 2593, 84845557)
    check_seed(2, 58321, 243701588629)
    binary_ring = PolynomialRing(GF(2), "u")
    u = binary_ring.gen()
    for level in (1, 2, 3, 4):
        assert (u**(2*3**level) + u**(3**level) + 1).is_irreducible()
    check_reciprocal_seed(1, 486, 459167617, 2593)
    check_reciprocal_seed(2, 4374, 2344061487986355189526424978686363, 403867)
    check_genuine_ratios()
    check_ray_and_inertia_controls()
    print("TP1–TP7 coefficients, charts, prime ratios, conductors and class-group fixtures verified.")

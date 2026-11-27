"""Exact prime, mirror, Dickson and isogeny degree fixtures for OPR1–OPR6."""

from sage.all import (
    CyclotomicField, GF, QQ, ZZ, PolynomialRing, companion_matrix,
    inverse_mod, matrix, pari, power_mod, primitive_root, valuation,
)

from verify_even_norm_curve_primitive_jacobian_twists import trace


def check_moving_primes():
    ell, r, a = 5, 5, 6250
    p, q = ZZ(38146972656250633751), ZZ(633751)
    A, B = a**r, -(p + q)//2
    assert p - q == 4*A and p.is_prime(proof=True) and q.is_prime(proof=True)
    assert p % 3125 == q % 3125 == (1-ell**4) % 3125
    Q = ZZ(100000000000000329838523899694676101749606900001)
    assert Q.is_prime(proof=True) and Q % 3125 == 1

    def power_residue(value, prime):
        residue = value.numerator()*inverse_mod(value.denominator(), prime)
        return power_mod(residue, (prime-1)//r, prime)

    assert power_residue(QQ(4*A*A)/p, Q) == 1
    assert power_residue(QQ(4*A*A)/q, Q) != 1
    for dyadic, P in ((4, ZZ(16000000001072118498347710176113893043556251)),
                      (20, ZZ(16000000001082997555278252349039541544006251))):
        assert P.is_prime(proof=True) and P % 3125 == 1
        assert power_residue(QQ(4)/p, P) != 1
        assert power_residue(QQ(q)/4, P) == 1
        numerator = A*A*P**(2*r) + B*P**r*Q**r + Q**(2*r)
        assert numerator < 0 and numerator % 32 == dyadic
        assert numerator % 3125 == ell**4
        assert numerator.valuation(p) == numerator.valuation(q) == 1
    even_lattice = matrix(ZZ, [[4*r, 0], [0, 4*r], [2*r, -2*r]])
    assert list(even_lattice.smith_form()[0].diagonal()) == [2*r, 4*r]


def check_mirror_discriminants():
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    for ell, level, D in ((5, 1, -11), (7, 1, 29), (5, 2, -11)):
        r = ell**level
        assert D % 4 == 1 and ZZ(D).is_squarefree()
        assert GF(ell)(D).is_square()
        cyclotomic = CyclotomicField(r)
        zeta = cyclotomic.gen()
        radicand = D*(zeta - zeta**(-1))**2
        polynomial = ring(radicand.minpoly()(x*x))
        assert polynomial.is_irreducible()
        degree = (ell-1)*ell**(level-1)
        ell_exponent = (level*(ell-1)-1)*ell**(level-1)
        assert polynomial.degree() == degree
        expected = ell**ell_exponent*abs(D)**(degree//2)
        actual = ZZ(pari(polynomial).nfdisc())
        assert abs(actual) == expected
        assert (actual > 0) == (D < 0)


def check_dickson_and_conductor():
    ring = PolynomialRing(ZZ, names=("T", "c"))
    T, c = ring.gens()
    for ell in (3, 5, 7, 11):
        dickson = [ring(2), T]
        for degree in range(2, ell+1):
            dickson.append(T*dickson[-1] - c*c*dickson[-2])
        half = (ell-1)//2
        square_factor = c**half + sum(c**(half-i)*dickson[i] for i in range(1, half+1))
        assert dickson[ell] - 2*c**ell == (T-2*c)*square_factor**2
    for ell in (5, 7):
        wild = 0
        previous = 0
        for level in (1, 2):
            r = ell**level
            cyclotomic = CyclotomicField(r)
            discriminant_exponent = valuation(cyclotomic.discriminant(), ell)
            expected_exponent = (level*(ell-1)-1)*ell**(level-1)
            assert discriminant_exponent == expected_exponent
            wild += discriminant_exponent - cyclotomic.degree() + 1
            full_half = r-1+wild
            assert full_half - previous == expected_exponent + 1
            previous = full_half

    for r, test_prime in ((5, 29), (7, 41)):
        field = GF(test_prime)
        assert test_prime % r == r-1
        assert any(trace(r, 1, parameter) != trace(r, -1, parameter)
                   for parameter in field if parameter not in (2, -field(2)))


def check_minimum_isogeny_degree():
    for ell, level, index in ((3, 2, 3), (5, 2, 5), (7, 1, 3)):
        r = ell**level
        cyclotomic = CyclotomicField(r)
        zeta = cyclotomic.gen()
        generator = primitive_root(r)
        stride = cyclotomic.degree()//index
        elements = [power_mod(generator, stride*i, r) for i in range(index)]
        assert len(set(elements)) == index
        theta = cyclotomic(1)
        for a in elements:
            theta *= zeta**a - zeta**(-a)
        assert theta.conjugate() == -theta
        b = theta/(zeta - zeta**(-1))
        assert b.conjugate() == b and b.is_integral()
        assert ell*abs(b.norm()) == ell**index
        for a in elements:
            assert theta == cyclotomic.hom([zeta**a])(theta)


def check_toric_centralizer():
    for r in (9, 25, 7):
        zeta = CyclotomicField(r).gen()
        real_polynomial = (zeta + zeta**(-1)).minpoly()
        toric_action = companion_matrix(real_polynomial)
        degree = toric_action.nrows()
        columns = []
        for i in range(degree):
            for j in range(degree):
                elementary = matrix(QQ, degree, degree, {(i, j): 1})
                commutator = toric_action*elementary - elementary*toric_action
                columns.append(commutator.list())
        commutator_matrix = matrix(QQ, columns).transpose()
        assert commutator_matrix.right_kernel().dimension() == degree


if __name__ == "__main__":
    check_moving_primes()
    check_mirror_discriminants()
    check_dickson_and_conductor()
    check_minimum_isogeny_degree()
    check_toric_centralizer()
    print("OPR1–OPR6 prime ratios, mirrors, conductors, isogeny degrees and centralizers verified.")

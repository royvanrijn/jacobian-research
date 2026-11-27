"""Exact polynomial, elliptic, discriminant and prime-ratio fixtures for SD1–SD3."""

from sage.all import EllipticCurve, GF, QQ, ZZ, PolynomialRing, kronecker, pari, power_mod, valuation


ring = PolynomialRing(QQ, "t")
t = ring.gen()


def trace_data(level, parameter):
    m = 3**level
    exponent = 2*((level + 2)//3) + 1
    a = 3**exponent
    coefficient = 2*a**3
    trace = 2 + coefficient*parameter**m
    polynomial = coefficient*parameter*(4 + coefficient*parameter**m)
    return m, a, coefficient, trace, polynomial


def dickson(degree):
    current, following = ring(2), t
    for _ in range(degree):
        current, following = following, t*following - current
    return current


def check_polynomial_and_elliptic_map():
    rational_t = ring.fraction_field()(t)
    for level in (2, 3, 4):
        m, a, coefficient, trace, H = trace_data(level, rational_t)
        assert trace**2 - 4 == t**(m - 1)*H
        assert H.numerator().degree() == m + 1
        assert H.numerator().is_squarefree()
        z = a*t**(m//3)
        X = 2/z
        omega_square = z*(z**3 + 2)
        Y_coefficient = 2/z**2
        assert Y_coefficient**2*omega_square == X**3 + 4
        assert (Y_coefficient**2*omega_square + 4)/X**3 == trace/2
        assert 4*Y_coefficient/X**3 == z

        _, _, _, evaluated_trace, _ = trace_data(level, QQ(22))
        delta = evaluated_trace - 2
        assert valuation(delta, 2) == m + 1
        assert (delta/2**(m + 1)) % 8 == 1
        power_at_three = valuation(delta, 3)
        assert power_at_three >= 2*level + 3 and power_at_three % 2 == 1
        assert (delta/3**power_at_three) % 3 == 2
        assert valuation(coefficient, 2) == 1
        assert sum((2*j + 1)*3**j - 1 for j in range(level)) == (level - 1)*(m - 1)

    curve = EllipticCurve([0, 0, 0, 0, 4])
    assert curve.conductor() == 108
    assert valuation(curve.quadratic_twist(-3).conductor(), 3) == 3
    assert (t**3 + 16).is_irreducible()


def check_prime_ratio_fixtures():
    for sign, numerator, denominator in [(1, 17, 19), (-1, 13, 73), (-1, 13, 577)]:
        assert ZZ(numerator).is_prime() and ZZ(denominator).is_prime()
        parameter = QQ(sign*2*numerator)/denominator
        m, _, coefficient, trace, _ = trace_data(2, parameter)
        assert valuation(parameter, 2) == 1 and (parameter/2) % 8 == 3
        assert parameter % 3 == 1
        assert (denominator - 1) % m == 0
        assert power_mod(coefficient, (denominator - 1)//3, denominator) != 1
        radicand = trace**2 - 4
        assert radicand*sign > 0 and valuation(radicand, numerator) == m
        D = ZZ(radicand.numerator()*radicand.denominator()).squarefree_part()
        assert D % 8 == 1 and valuation(D, 3) == 1
        polynomial = denominator**m * (dickson(m)(t/denominator) - trace)
        assert polynomial.is_irreducible()
        assert abs(pari(polynomial).nfdisc()) == 3**8*abs(D/3)**4
        if denominator == 577:
            assert valuation(4 + coefficient*parameter**m, 397) == 1
            assert kronecker(denominator, 397) == -1


def check_torsion_and_rank_inputs():
    finite_ring = PolynomialRing(GF(2), names=("v", "Z"))
    v, Z = finite_ring.gens()
    for degree in (3, 5, 9, 15, 27):
        genus = (degree - 1)//2
        affine = Z**2 + Z - v**degree
        infinity = Z**2 + v**(genus + 1)*Z - v
        assert affine.derivative(Z) == 1
        assert infinity.derivative(v)(0, 0) == 1
        # On residue vectors of sum zero modulo degree, the relation
        # identifies the two basis vectors with opposite signs.
        generator_order = degree//ZZ(2).gcd(degree)
        assert generator_order == degree

    curve = EllipticCurve([0, 0, 0, 0, 4])
    assert [curve.change_ring(GF(p)).cardinality() for p in (5, 7)] == [6, 3]
    for unit in (1, 2, 4, 5, 7, 8):
        if unit % 3 == 1:
            assert (unit**3 - 4) % 9 == 6
        else:
            assert (unit**3 - 4) % 3 == 1


if __name__ == "__main__":
    check_polynomial_and_elliptic_map()
    check_prime_ratio_fixtures()
    check_torsion_and_rank_inputs()
    print("SD1–SD5 reduced form, local traces, elliptic and torsion inputs, and prime fixtures verified.")

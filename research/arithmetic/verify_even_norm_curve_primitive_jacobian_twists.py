"""Exact quotient, cyclotomic discriminant and Frobenius fixtures for ENJ1–ENJ4."""

from sage.all import CyclotomicField, GF, QQ, PolynomialRing, ZZ, valuation

from verify_two_moving_prime_classes_even_norm_curves import check_genuine_ratios


def dickson(value, degree):
    previous, current = 2, value
    for index in range(2, degree + 1):
        previous, current = current, value*current - previous
    return current


def check_quotients():
    base = PolynomialRing(QQ, "w").fraction_field()
    w = base.gen()
    s = w + 1/w
    parameter_ring = PolynomialRing(base, "lam")
    lam = parameter_ring.gen()
    for level in (1, 2, 3):
        r = 3**level
        assert dickson(s, r) == w**r + w**(-r)
        for sign in (-1, 1):
            assert (w + sign)**2*(w**(2*r) + lam*w**r + 1)/w**(r+1) == (
                (s + 2*sign)*(dickson(s, r) + lam)
            )
        assert (r - 1)//2 - (r//3 - 1)//2 == r//3


def check_cyclotomic_contributions():
    wild_sum = 0
    preceding_exponent = 0
    for level in (1, 2, 3):
        r = 3**level
        field = CyclotomicField(r)
        root = field.gen()
        prime = field.primes_above(3)[0]
        assert prime.ramification_index() == 2*r//3
        assert prime.residue_class_degree() == 1
        discriminant_exponent = valuation(field.discriminant(), 3)
        assert discriminant_exponent == (2*level - 1)*r//3
        for j in range(1, level+1):
            difference = root**(3**(level-j)) - 1
            assert field.ideal(difference).valuation(prime) == 3**(level-j)
        wild_sum += discriminant_exponent - field.degree() + 1
        assert wild_sum == (level-2)*r + level + 2
        half_exponent = r - 1 + wild_sum
        assert half_exponent - preceding_exponent == (2*level-1)*r//3 + 1
        preceding_exponent = half_exponent


def quadratic_character(value):
    if not value:
        return 0
    return 1 if value.is_square() else -1


def trace(degree, sign, lam):
    # There are two rational points at infinity because the polynomial is monic.
    return -1 - sum(quadratic_character((s + 2*sign)*(dickson(s, degree) + lam))
                    for s in lam.parent())


def check_primitive_frobenius():
    a = 4374
    p, q = 2344061487986355189526424978686363, 403867
    assert ZZ(p).is_prime(proof=True) and ZZ(q).is_prime(proof=True)
    A, B = a**9, -(p + q)//2
    assert p - q == 4*A and (2*A*3*p*q) % 17 != 0
    field = GF(17)
    lam = field(B)/field(A)
    assert lam == 1

    assert [trace(9, sign, lam) - trace(3, sign, lam) for sign in (1, -1)] == [3, -3]


def check_simultaneous_nonisogeny():
    a = 4374
    p, q = ZZ(2344061487986355189526425125084143), ZZ(146801647)
    A, B = a**9, -(p + q)//2
    assert p.is_prime(proof=True) and q == 17*31*278561
    assert ZZ(278561).is_prime(proof=True)
    assert p - q == 4*A and p.gcd(q) == 1
    assert p % 2187 == q % 2187 == (1 - 3**6) % 2187
    for prime, degree, expected in ((53, 3, (-6, 6)), (71, 9, (3, -3))):
        field = GF(prime)
        extension = GF(prime**2, "v")
        generator = extension.multiplicative_generator()**(prime-1)
        assert generator.multiplicative_order() == prime+1
        lam = field(B)/field(A)
        assert lam == 1
        traces = tuple(trace(degree, sign, lam) - trace(degree//3, sign, lam)
                       for sign in (1, -1))
        assert traces == expected

        coefficients = dict.fromkeys(field, 0)
        for exponent, multiplier in ((degree, -degree), (degree//3, degree//3)):
            for i in range(0, prime+1, exponent):
                element = generator**i
                image = element + element.inverse()
                s = field(image.polynomial()[0])
                assert extension(s) == image
                coefficients[s] += multiplier*(-1)**i
        assert sum(coefficients.values()) == 0
        for parameter in field:
            difference = 0
            for s in field:
                weight = quadratic_character(s-2) - quadratic_character(s+2)
                difference += weight*(quadratic_character(dickson(s, degree) + parameter)
                                     - quadratic_character(dickson(s, degree//3) + parameter))
            assert difference == sum(b*quadratic_character(s+parameter)
                                     for s, b in coefficients.items())


if __name__ == "__main__":
    check_genuine_ratios()
    check_quotients()
    check_cyclotomic_contributions()
    check_primitive_frobenius()
    check_simultaneous_nonisogeny()
    print("ENJ1–ENJ4 genuine ratios, quotients, conductors and nonisogeny fixtures verified.")

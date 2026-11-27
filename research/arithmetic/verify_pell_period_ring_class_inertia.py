"""Exact residue groups and prime and class-number certificates for PI1–PI4."""

from sage.all import GF, QQ, ZZ, PolynomialRing, gcd, kronecker, valuation


def verify_dyadic_quotients():
    for radicand in (1, 2, 3):
        def scalar_class(real, radical):
            return min((real % 4, radical % 4), (-real % 4, -radical % 4))

        def multiply(first, second):
            a, b = first
            c, d = second
            return scalar_class(a*c+radicand*b*d, a*d+b*c)

        units = {
            scalar_class(a, b)
            for a in range(4)
            for b in range(4)
            if (a*a-radicand*b*b) % 2
        }
        identity = scalar_class(1, 0)
        fourth = scalar_class(1, 2)
        assert len(units) == 4
        assert multiply(fourth, fourth) == identity
        if radicand == 2:
            generator = scalar_class(1, 1)
            assert multiply(generator, generator) == fourth
        else:
            other = scalar_class(0, 1)
            assert other not in (identity, fourth)
            assert multiply(other, other) == identity
            assert all(multiply(unit, unit) == identity for unit in units)


def verify_simultaneous_prime_fixture():
    support = (3, 5)
    B = ZZ(15)
    m = ZZ(45)
    s = ZZ(479)
    c = ZZ(511)
    A = 2*(B*c-1)
    q = ZZ(11690712001)
    assert s.is_prime(proof=True) and q.is_prime(proof=True)
    assert B*c == 2**4*s+1 and gcd(s-1, m) == 1
    assert pow(q, m, 8*A) == (c*c+2*A) % (8*A)
    assert q % B == 1 and gcd(q, c*(B*c-2)*A) == 1
    n = (q**m-c*c)//A
    raw = n*(B*B*n+2)
    assert n > 0 and A*n+c*c == q**m
    assert valuation(raw, 2) == 3
    for p, ell, omitted in [(3, 29, 10), (5, 59, 15)]:
        assert ZZ(ell).is_prime(proof=True) and ell % B == B-1
        assert c % (p*p) == 1+2*p
        assert valuation(raw, p) == 1
        field = GF(ell)
        ring = PolynomialRing(field, "t")
        t = ring.gen()
        previous, chebyshev = ring(1), t
        for _ in range(2, p+1):
            previous, chebyshev = chebyshev, 2*t*chebyshev-previous
        image = {chebyshev(x) for x in field}
        assert ell-len(image) == (ell+1)*(p-1)//(2*p)
        assert (1+B*B*n) % ell == omitted and field(omitted) not in image
        assert valuation(B, p) == 1 and (B*n+c) % p != 0
    assert (1+B*B*n)**2-B*B*raw == 1
    assert (B*n+c)**2-raw == q**m
    assert raw % q != 0 and (B*n+c)**2 % q == raw % q


def verify_degree_twenty_five_fixture():
    D = ZZ(5615)
    discriminant = 4*D
    assert D == 5*1123 and ZZ(1123).is_prime(proof=True)
    assert valuation(D, 5) == valuation(15, 5) == 1
    for b in range(1, 15):
        assert (75*b-1)**2 < D*b*b+1 < (75*b)**2
    assert 1124**2-D*15**2 == 1
    assert 1124**2 < D*15**2+2 and 1124+15*74 > 3**7
    assert 542**2-D*9**2 == -11**5
    assert (542-9*4) % 11 == 0 and (542+9*4) % 11 != 0
    assert 4**2 % 11 == D % 11 and 9*pow(542, -1, 5) % 5 != 0

    scale = 10**8
    partial_sum = 0
    maximum_partial_sum = 0
    rounded_upper = 0
    for n in range(1, discriminant+1):
        character = int(kronecker(discriminant, n))
        partial_sum += character
        maximum_partial_sum = max(maximum_partial_sum, abs(partial_sum))
        rounded_upper += -(-character*scale//n)
    assert partial_sum == 0 and maximum_partial_sum == 69
    assert rounded_upper == 102972285
    upper = QQ(rounded_upper)/scale+QQ(maximum_partial_sum)/discriminant
    assert upper == QQ(23196575211)/22460000000 < 2
    assert discriminant < 150**2 and QQ(150)/7 < 25
    assert (25-1)//2 == 12 and (5-1)*25//5 == 20


if __name__ == "__main__":
    verify_dyadic_quotients()
    verify_simultaneous_prime_fixture()
    verify_degree_twenty_five_fixture()
    print("PI1–PI4 dyadic groups, simultaneous actual prime, and degree-25 data verified")

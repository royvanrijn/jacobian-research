"""Replay the weighted quartic identities and the two cyclotomic prime triangles."""

from sage.all import GF, QQ, ZZ, EllipticCurve, PolynomialRing, QuadraticField, pari


def check_identities():
    parameters = PolynomialRing(QQ, "kpol")
    coefficient_field = parameters.fraction_field()
    k = coefficient_field(parameters.gen())
    polynomial_ring = PolynomialRing(coefficient_field, "tpol")
    functions = polynomial_ring.fraction_field()
    t = functions(polynomial_ring.gen())
    w = k * (k + 1)
    m = 2 * k + 1
    quartic = -3 * t**4 + 4 * m * t**3 - 6 * w * t**2 + w**2
    assert quartic.numerator().discriminant() == -6912 * w**4
    a = t**3
    b = k**2 * (t - k - 1)**3
    c = (k + 1)**2 * (t - k)**3
    assert a**2 + b**2 + c**2 - 2 * (a * b + a * c + b * c) == w**2 * quartic
    r = t * (t - k) * (t - k - 1) / (2 * ((w + 1) * t**2 - w * m * t + w**2))
    target = 1 - 24 * (w + 1) * r**2 + 32 * m * (w - 2) * r**3 - 48 * (w + 1)**2 * r**4
    assert target == 4 * r.derivative(t)**2 * quartic


def check_prime_fixtures():
    field = QuadraticField(-3, "s")
    zeta = (field.gen() - 1) / 2
    e = EllipticCurve(field, [0, 360**2])
    t_point = e(360, -360 * 19)
    j_point = e(zeta * 360, t_point[1]) - e(zeta**2 * 360, t_point[1])
    pole = e(zeta * 360, -t_point[1])
    assert 2 * pole - t_point == -j_point
    fixtures = [
        (3, 56611, 5016, 732435231528),
        (5, 8826476898331, 732435231528, 5486777363970309410097697565862059160),
    ]
    for n, u, h, delta in fixtures:
        prime_1, prime_2, prime_3 = map(ZZ, (u, u - 10 * h, u - 9 * h))
        primes = (prime_1, prime_2, prime_3)
        assert all(q.is_prime(proof=True) and q % 3 == 1 for q in primes)
        assert len(set(primes)) == 3
        assert (prime_1 + 9 * prime_2 - 10 * prime_3) == 0
        a, b, c = (prime_1**3, 9**2 * prime_2**3, 10**2 * prime_3**3)
        assert (a - b) % 10 == (a - c) % 9 == 0
        heron = ZZ(4 * a * c - (a + c - b)**2) / (10 * 9)**2
        assert heron == 3 * delta**2
        assert all(ZZ(heron).gcd(q) == 1 for q in primes)
        point = pole + ((n + 1) // 2) * j_point
        assert 2 * point - t_point == n * j_point
        assert 180 * (point[0] - 360) / (point[1] - 360 * 19) == QQ(u) / h
        for prime, mark in zip(primes, (0, 1, -1)):
            reduction = EllipticCurve(GF(prime), [0, 360**2])
            assert (3 * n**2 + 1) * reduction(360, -360 * 19) == (2 * mark) * reduction(0, 360)
            translated = point - mark * e(0, 360)
            local_exponents = [translated[0].valuation(ideal) for ideal in field.primes_above(prime)]
            assert sorted(max(0, -v) for v in local_exponents) == [0, 2]

    rational_curve = EllipticCurve(QQ, [0, 2025]).minimal_model()
    assert tuple(pari(rational_curve).ellrank()[:2]) == (1, 1)
    print("Verified the weighted identities, two prime triangles and rank-one descent.")


if __name__ == "__main__":
    check_identities()
    check_prime_fixtures()

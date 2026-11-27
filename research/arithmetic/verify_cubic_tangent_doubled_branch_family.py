"""Replay the doubled cubic correspondence and its exact arithmetic specialization."""

from sage.all import EllipticCurve, PolynomialRing, QQ, QuadraticField


def check_rational_model():
    ring = PolynomialRing(QQ, ("m", "s", "Z"))
    m, s, Z = ring.gens()
    denominator = s * Z**2 - m
    numerator = Z * (Z + m * s)
    elimination = ((numerator**3 + m * denominator**3) * (Z**3 + m**2)
                   - m * denominator * (numerator * Z - m * denominator)**2)
    assert elimination == Z**2 * (Z**3 + m**2)**2 * (Z * (1 + m * s**3) + 2 * m * s)

    field = ring.fraction_field()
    m, s = field(m), field(s)
    w = m * s**3
    x = 2 * m * s**2 / (w - 1)
    z = -2 * m * s / (w + 1)
    branch = (w**2 - 1) * (w**2 + 4 * w - 1)
    assert x**3 + m == m * branch / (w - 1)**4
    assert z**3 + m**2 == m**2 * branch / (w + 1)**4
    assert (x**3 + m) * (z**3 + m**2) == m * (x * z - m)**2
    assert (z**2 + m * x) / (z * (x * z - m)) == s

    parameters = PolynomialRing(QQ, "u")
    u = parameters.gen()
    tangent = u**2 * ((3 + 3 * u + u**2)**2 - QQ(9) / 4 * (2 + u)**2)
    assert tangent.valuation() == 3
    assert (u * (u - 3)**2 - 4) == (u - 1)**2 * (u - 4)


def check_specialization():
    field = QuadraticField(-2807, "r")
    r = field.gen()
    s = QQ(1) / 4
    b = 6
    m = b**2
    w = m * s**3
    Y = 5 * r / 256
    assert Y**2 == (w**2 - 1) * (w**2 + 4 * w - 1)
    x = 2 * m * s**2 / (w - 1)
    z = -2 * m * s / (w + 1)
    y = b * Y / (w - 1)**2
    v = m * Y / (w + 1)**2
    assert x == -QQ(72) / 7 and z == -QQ(288) / 25

    profiles = [
        (b, (x, y), (12, -42), QQ(20) / 49, QQ(32) / 343, -QQ(624) / 1715),
        (m, (z, v), (72, -612), QQ(140) / 289, QQ(3087) / 4913, QQ(10231200) / 4721393),
    ]
    for coefficient, point, half, scale, k, expected in profiles:
        curve = EllipticCurve(field, [0, coefficient**2])
        translated = curve(*point) + curve(*half)
        w_k = k * (k + 1)
        m_k = 2 * k + 1
        t = 2 * w_k * (scale**2 * translated[0] - 4 * w_k) / (
            scale**3 * translated[1] - 4 * w_k * m_k)
        assert t == expected
        marks = (expected, expected - k, expected - k - 1)
        assert all(mark > 0 for mark in marks) or all(mark < 0 for mark in marks)

    alpha, beta = y + b, v + m
    assert alpha == 6 + 30 * r / 49 and beta == 36 + 36 * r / 125
    assert alpha.norm() == -x**3 and beta.norm() == -z**3
    for prime, expected in [(2, [(2, 3), (7, 12)]), (3, [(1, 2), (5, 4)])]:
        valuations = sorted((alpha.valuation(p), beta.valuation(p)) for p in field.primes_above(prime))
        assert valuations == expected
    assert all(not value.is_nth_power(3) for value in (alpha, beta, alpha * beta, alpha * beta**2))


def check_split_prime_classes():
    field = QuadraticField(-747375839, "r")
    w = 36 * (QQ(5) / 18)**3
    branch = (w**2 - 1) * (w**2 + 4 * w - 1)
    assert (branch / -747375839).is_square()
    group = field.class_group(proof=True)
    assert group.invariants() == (5736, 2, 2)
    for q in (2, 3):
        primes = field.primes_above(q)
        assert len(primes) == 2
        assert all(group(p).exponents()[0] % 3 == 0 for p in primes)


def check_third_quotient():
    ring = PolynomialRing(QQ, "w")
    base = ring.fraction_field()
    w = base(ring.gen())
    extension_ring = PolynomialRing(base, "Y")
    field = base.extension(extension_ring.gen()**2 - (w**2 - 1) * (w**2 + 4 * w - 1), "Y")
    w, Y = field(w), field.gen()
    X = (w + 1) / (w - 1)
    U = Y / (w - 1)**2
    assert U**2 == X**3 + X**2 - X
    alpha_over_b = U + 1
    beta_over_b_squared = Y / (w + 1)**2 + 1
    gamma = U + 2 * X - 1
    assert alpha_over_b * beta_over_b_squared * gamma == (Y / (w**2 - 1) + X)**3

    curve = EllipticCurve(QQ, [0, 1, 0, -1, 0])
    assert curve(1, 1).order() == 3
    division = curve.division_polynomial(3)
    x = division.parent().gen()
    assert division == (x - 1) * (3 * x**3 + 7 * x**2 + x + 1)
    assert not (3 * x**3 + 7 * x**2 + x + 1).roots(ring=QQ)


if __name__ == "__main__":
    check_rational_model()
    check_specialization()
    check_split_prime_classes()
    check_third_quotient()
    print("Verified the doubled cubic model, exact positive specialization, split prime classes, and third norm identity.")

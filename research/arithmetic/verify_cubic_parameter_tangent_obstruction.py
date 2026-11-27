"""Verify exact low-degree polarization arithmetic, the genus-four fixture, and the tetrahedral invariant."""

from sage.all import CyclotomicField, EllipticCurve, Matrix, PolynomialRing, QQ, QuadraticField, ZZ


def check_polarizations():
    first_row = 2
    second_row = 3
    assert [Matrix(ZZ, [[first_row, n], [n, second_row]]).det() for n in range(3)] == [6, 5, 2]
    assert [Matrix(ZZ, [[3, n], [n, 3]]).det() for n in range(3)] == [9, 8, 5]

    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    assert x**4 + 4 * x**3 - 8 * x + 4 == (x**2 + 2 * x - 2)**2
    assert not (x**2 + 2 * x - 2).roots(ring=QQ)


def check_cubic_projections():
    ring = PolynomialRing(QQ, "s")
    base = ring.fraction_field()
    s = base(ring.gen())
    extension_ring = PolynomialRing(base, "Y")
    field = base.extension(extension_ring.gen()**2 - s**9 - 36, "Y")
    s, Y = field(s), field.gen()
    curve = EllipticCurve(field, [0, 36])
    point = curve(s**3, Y)
    profiles = [
        (curve(12, -42), QQ(32) / 343, -100 * (s**3 - 12) / (343 * (3 * s**3 + 13))),
        (curve(4, 10), -QQ(343) / 375, 196 * (s**3 - 4) / (375 * (3 * s**3 + 13))),
    ]
    for half, k, expected in profiles:
        center = 2 * half
        scale = center[0] / 6
        w = k * (k + 1)
        m = 2 * k + 1
        assert scale**2 * center[0] == 4 * w
        assert scale**3 * center[1] == -4 * w * m
        translated = point + half
        actual = 2 * w * (scale**2 * translated[0] - 4 * w) / (scale**3 * translated[1] - 4 * w * m)
        assert actual == expected

    first_t = -QQ(750) / 16121
    second_t = QQ(2254) / 17625
    assert first_t == -QQ(100) * (27 - 12) / (343 * (3 * 27 + 13))
    assert second_t == QQ(196) * (27 - 4) / (375 * (3 * 27 + 13))
    assert first_t < 0 and second_t > 1 - QQ(343) / 375

    specialization = QuadraticField(2191, "r")
    edge = 6 + 3 * specialization.gen()
    assert edge.norm() == -3**9
    assert sorted(edge.valuation(prime) for prime in specialization.primes_above(3)) == [1, 8]


def check_tetrahedral_invariant():
    cyclotomic = CyclotomicField(3)
    zeta = cyclotomic.gen()
    ring = PolynomialRing(cyclotomic, "s")
    field = ring.fraction_field()
    s = field(ring.gen())

    def involution(variable):
        return (variable + 2) / (variable - 1)

    def invariant(variable):
        return (variable**12 + 264 * variable**6 + 440 * variable**3 + 24) / (variable**3 - 1)**3

    assert invariant(s) == s**3 + sum(involution(zeta**i * s)**3 for i in range(3))
    assert invariant(zeta * s) == invariant(s)
    assert invariant(involution(s)) == invariant(s)
    parameters = PolynomialRing(QQ, "lambda")
    weight = parameters.gen()
    coordinate = PolynomialRing(parameters, "u")
    u = coordinate.gen()
    branch = (u + 1)**4 + 264 * (u + 1)**2 + 440 * (u + 1) + 24 - weight * u**3
    e, d, c, b, a = branch.list()
    assert 12 * a * e - 3 * b * d + c**2 == 2916 * (weight + 24)
    at_special_level = PolynomialRing(QQ, "u")([coefficient(-24) for coefficient in branch.list()])
    u = at_special_level.parent().gen()
    assert at_special_level == (u + 1) * (u + 9)**3


if __name__ == "__main__":
    check_polarizations()
    check_cubic_projections()
    check_tetrahedral_invariant()
    print("Verified the polarization determinants, genus-four fixture, and tetrahedral branch invariant.")

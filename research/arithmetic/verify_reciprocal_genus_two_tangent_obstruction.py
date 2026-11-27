"""Replay the two reciprocal elliptic projections and their norm-class specialization."""

from sage.all import EllipticCurve, QQ, PolynomialRing, QuadraticField


def check_projections():
    ring = PolynomialRing(QQ, "s")
    fraction_field = ring.fraction_field()
    s = fraction_field(ring.gen())
    extension_ring = PolynomialRing(fraction_field, "U")
    U = extension_ring.gen()
    field = fraction_field.extension(U**2 - s**6 - 36, "Y")
    s, Y = field(s), field.gen()
    examples = [
        (QQ(-7) / 2, QQ(32) / 343, QQ(20) / 49, -100 * (s**2 - 12) / (343 * (3 * s**2 + 13))),
        (QQ(-17) / 7, QQ(3087) / 4913, QQ(140) / 289, 176400 * (2 * s**2 - 1) / (4913 * (73 * s**2 + 108))),
    ]
    for index, (r, k, scale, target) in enumerate(examples):
        w = k * (k + 1)
        m = 2 * k + 1
        assert k == -(2 * r + 1) / (r**3 * (r + 2))
        branch_value = -(1 + r)**2 / (r**3 * (r + 2))
        a = w / branch_value**2
        curve = EllipticCurve(field, [0, (4 * w)**2])
        half_center = curve(2 * w * (a - 1), 2 * w * ((a**2 + 3) * branch_value - 2 * m))
        assert 2 * half_center == curve(4 * w, -4 * w * m)
        if index == 0:
            point = curve(scale**2 * s**2, scale**3 * Y)
        else:
            point = curve(scale**2 * 36 / s**2, scale**3 * 36 * Y / s**3)
        translated = point + half_center
        actual = 2 * w * (translated[0] - 4 * w) / (translated[1] - 4 * w * m)
        assert actual == target


def check_norm_classes():
    field = QuadraticField(37, "s")
    s = field.gen()
    a = s + 6
    b = 36 * (s + 1)
    assert a.norm() == -1 and b.norm() == -36**3
    assert all(not value.is_nth_power(3) for value in (a, b, a * b, a * b**2))
    assert sorted(a.valuation(prime) for prime in field.primes_above(3)) == [0, 0]
    assert sorted(b.valuation(prime) for prime in field.primes_above(3)) == [2, 4]


if __name__ == "__main__":
    check_projections()
    check_norm_classes()
    print("Verified the reciprocal projections and the independent norm classes over Q(sqrt(37)).")

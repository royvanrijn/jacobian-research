"""Run with Sage's Python to replay the marked fibres and an exact field collision."""

from sage.all import GF, QQ, ZZ, EllipticCurve, PolynomialRing, QuadraticField


def tau(point):
    return (point[0] + 1) / point[1]


def main():
    base = EllipticCurve(QQ, [0, 1])
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    quartic = 1 + 6 * t**2 - 3 * t**4
    assert quartic.is_irreducible()

    for prime, expected_counts in [(5, (6, 36)), (7, (12, 48))]:
        counts = tuple(base.change_ring(GF(prime**degree)).cardinality() for degree in (1, 2))
        assert counts == expected_counts

    torsion = base(-1, 0)
    third = base(0, 1)
    assert 2 * torsion == base(0)
    assert 3 * third == base(0)
    assert tau(third) == tau(torsion - third) == 1
    assert tau(-third) == tau(torsion + third) == -1
    assert 5 * torsion == torsion
    assert 5 * third == -third

    p, h = ZZ(17), ZZ(6)
    d = h**4 + 6 * p**2 * h**2 - 3 * p**4
    assert d == -186843 and d.is_squarefree()
    vertices = (p - h, p, p + h)
    assert all(vertex.is_prime() for vertex in vertices)
    field = QuadraticField(d, "s")
    s = field.gen()
    curve = base.change_ring(field)
    t0 = QQ(p) / h
    x = (s / h**2 + 1 + t0**2) / (2 * t0**2)
    point = curve(x, (x + 1) / t0)
    assert tau(point) == t0
    conjugate_x = (-s / h**2 + 1 + t0**2) / (2 * t0**2)
    conjugate = curve(conjugate_x, (conjugate_x + 1) / t0)
    assert point + conjugate == curve(-1, 0)
    twisted = 2 * point - curve(-1, 0)
    assert (d * twisted[0], d * s * twisted[1]) == (311691, -154137708)

    fixtures = [
        (3, 121647372209, 82553876370, (1, 4301, 1)),
        (5, 6481558012053794672265292343537, 1613150365317644553963830999454, (23, 17, 11)),
    ]
    source_product = (p - h) * p * (p + h)
    for multiple, numerator, denominator, expected in fixtures:
        ratio = abs(QQ(tau(multiple * point)))
        assert ratio == QQ(numerator) / denominator > 1
        assert (quartic(ratio) / d).is_square()
        a, b = ratio.numerator(), ratio.denominator()
        assert b % 6 == 0
        assert tuple(mark.gcd(source_product) for mark in (a - b, a, a + b)) == expected
        assert not a.is_prime()
    print("Verified the marked fibres and the two relaxed same-field examples.")


if __name__ == "__main__":
    main()

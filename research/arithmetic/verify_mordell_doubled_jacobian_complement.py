"""Replay the complementary morphism, cubic norm identity, and third local ratios."""

from sage.all import EllipticCurve, PolynomialRing, QQ


def check_complementary_projection():
    parameter = PolynomialRing(QQ, "t")
    base = parameter.fraction_field()
    t = base(parameter.gen())
    w = t**3
    extension = PolynomialRing(base, "Y")
    field = base.extension(extension.gen()**2 - (w**2 - 1) * (w**2 + 4 * w - 1), "Y")
    t, Y = field(t), field.gen()
    w = t**3
    u, v = t - 1 / t, Y / w
    assert v**2 == (u**3 + 3 * u) * (u**3 + 3 * u + 4)
    denominator = (u + 1) * (u**2 + 3)
    X = 6 * u / denominator
    V = 3 * (2 * u**2 + 3 * u + 3) * v / denominator**2
    assert V**2 == X**3 - 3 * X**2 / 4 + 6 * X
    xi = V + 3 * X / 2 + 1
    first = 1 + Y / (w - 1)**2
    second = 1 + Y / (w + 1)**2
    cube_root = -(t + 1) * (t**2 - 3 * t - 1) / (2 * (t**2 + t + 1))
    cube_root += Y / (2 * (t - 1) * (t**2 + t - 1) * (t**2 + t + 1))
    assert xi * first / second == cube_root**3


def check_elliptic_models():
    third = EllipticCurve(QQ, [0, 1, 0, -1, 0])
    fourth = EllipticCurve(QQ, [0, -QQ(3) / 4, 0, 6, 0])
    assert third.j_invariant() == QQ(16384) / 5
    assert fourth.j_invariant() == QQ(804357) / 500
    assert third.conductor() == 20 and fourth.conductor() == 90
    assert third.ap(13) == 2 and fourth.ap(13) == -4

    point = fourth(1, QQ(5) / 2)
    assert point.order() == 3
    polynomial = fourth.division_polynomial(3)
    X = polynomial.parent().gen()
    assert polynomial == 3 * (X - 1) * (X**3 + 12 * X + 12)
    assert not (X**3 + 12 * X + 12).roots(ring=QQ)

    isogeny = third.isogeny(third(1, 1))
    target = isogeny.codomain().global_minimal_model()
    isogeny = isogeny.codomain().isomorphism_to(target) * isogeny
    for q, expected in [(2, QQ(1) / 3), (3, 1), (5, 1)]:
        tamagawa_ratio = QQ(target.tamagawa_number(q)) / third.tamagawa_number(q)
        assert tamagawa_ratio * q**isogeny.scaling_factor().valuation(q) == expected


if __name__ == "__main__":
    check_complementary_projection()
    check_elliptic_models()
    print("Verified the complementary elliptic map, its cube relation, and the third isogeny data.")

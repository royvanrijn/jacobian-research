"""Check the fixed conductors, local isogeny data, and generic-rank Frobenius bound."""

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, pari


def check_base_curves():
    expected = [(6, 972, QQ(1) / 3, QQ(1) / 3), (36, 243, QQ(1), QQ(1))]
    for coefficient, conductor, at_two, at_three in expected:
        model = EllipticCurve(QQ, [0, coefficient**2])
        assert model.conductor() == conductor
        minimal = model.global_minimal_model()
        kernel = model.isomorphism_to(minimal)(model(0, coefficient))
        isogeny = minimal.isogeny(kernel)
        target = isogeny.codomain().global_minimal_model()
        minimal_isogeny = isogeny.codomain().isomorphism_to(target) * isogeny
        differential = minimal_isogeny.scaling_factor()
        for q, expected_ratio in [(2, at_two), (3, at_three)]:
            ratio = QQ(target.tamagawa_number(q)) / minimal.tamagawa_number(q)
            assert ratio * q**differential.valuation(q) == expected_ratio
        assert model.quadratic_twist(-2807).conductor() == conductor * 2807**2


def check_frobenius():
    ring = PolynomialRing(GF(13), "s")
    s = ring.gen()
    w = 36 * s**3
    branch = (w**2 - 1) * (w**2 + 4 * w - 1)
    assert branch.degree() == 12 and branch.is_squarefree()
    polynomials = PolynomialRing(ZZ, "T")
    T = polynomials.gen()
    frobenius = polynomials([ZZ(value) for value in pari(branch).hyperellcharpoly().Vecrev()])
    first = T**2 - 5 * T + 13
    second = T**2 + 7 * T + 13
    assert frobenius == first * second * (T**2 - 2 * T + 13) * (
        T**4 - 4 * T**3 + 3 * T**2 - 52 * T + 169)
    for coefficient, polynomial in [(6, first), (36, second)]:
        curve = EllipticCurve(QQ, [0, coefficient**2])
        assert polynomial == T**2 - curve.ap(13) * T + 13
        assert frobenius.valuation(polynomial) == 1


if __name__ == "__main__":
    check_base_curves()
    check_frobenius()
    print("Verified the conductor pair, local isogeny ratios, and generic-rank Frobenius factors.")

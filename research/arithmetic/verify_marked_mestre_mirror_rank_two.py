"""Check the two marked Mestre sections and their local independence."""

from sage.all import EllipticCurve, GF, Integers, PolynomialRing, QQ, ZZ, pari


def abscissas(u):
    short = -37 * (u**4 + u**2 + 1) / (75 * u**2 * (u**2 + 1))
    return short - 3, u**2 * short - 3


def radicand(x):
    return x**3 + (3 * x + 17)**2


def check_twist():
    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    u = field(ring.gen())
    first, second = abscissas(u)
    A = u**4 + u**2 + 1
    B = u**2 + 1
    twist = -75 * 37 * B * (37**2 * A**3 + 75**3 * u**4 * B**2)
    assert radicand(first) * (75**2 * u**3 * B**2)**2 == twist
    assert radicand(second) == u**6 * radicand(first)

    reduced = ring(twist).change_ring(GF(7))
    assert reduced.degree() == 14 and reduced.is_squarefree()
    polynomials = PolynomialRing(ZZ, "T")
    T = polynomials.gen()
    frobenius = polynomials([ZZ(n) for n in pari(reduced).hyperellcharpoly().Vecrev()])
    elliptic_factor = T**2 + T + 7
    assert frobenius == elliptic_factor**2 * (T**4 + T**3 - 2 * T**2 + 7 * T + 49)**2
    model = EllipticCurve(QQ, [0, 9, 0, 102, 289])
    assert model.ap(7) == -1 and model.j_invariant() == QQ(108000000) / 63869
    division = model.division_polynomial(3)
    x = division.parent().gen()
    assert division == 3 * x * (x**3 + 12 * x**2 + 204 * x + 1156)
    assert (x**3 + 12 * x**2 + 204 * x + 1156).is_irreducible()


def check_local_conditions():
    x, _ = abscissas(QQ(1))
    assert x == -QQ(187) / 50
    residue = Integers(9)
    ring = PolynomialRing(residue, "I")
    I = ring.gen()
    quadratic = ring.quotient(I**2 + 1)
    i = quadratic.gen()
    assert residue(radicand(x)) == residue(-4)
    assert residue(3 * x + 17) == residue(2)
    assert (2 + 2 * i)**8 == 1 and (2 - 2 * i)**8 == 1

    assert abscissas(GF(17)(2)) == (GF(17)(11), GF(17)(2))
    examples = [(7, 2, 3, (2, 2)), (19, 13, 3, (7, 11))]
    for prime, parameter, ordinate, expected in examples:
        field = GF(prime)
        u, y = field(parameter), field(ordinate)
        x1, x2 = abscissas(u)
        assert y**2 == radicand(x1)
        assert (u**3 * y)**2 == radicand(x2)
        evaluations = (y + 3 * x1 + 17, u**3 * y + 3 * x2 + 17)
        assert tuple(int(a**((prime - 1) // 3)) for a in evaluations) == expected

    parameter = QQ(271) / 433
    for modulus, expected in [(2, 1), (27, 1), (17, 2), (7, 2), (19, 13)]:
        assert Integers(modulus)(parameter) == expected
    assert radicand(abscissas(parameter)[0]) < 0


if __name__ == "__main__":
    check_twist()
    check_local_conditions()
    print("Verified two Mestre sections, generic rank two, and the local cube and independence residues.")

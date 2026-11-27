"""Replay the local cube conditions and an exact pair of mirror class groups."""

from sage.all import GF, Integers, PolynomialRing, QQ, QuadraticField, pari


def check_mirror_fixture():
    s = -QQ(5) / 21
    w = 225 * s**3
    assert w == -QQ(3125) / 1029
    assert w.valuation(2) == 0 and w.valuation(7) < 0
    residue = Integers(9)
    assert residue(1 / w) == residue(3)
    radicand = (w**2 - 1) * (w**2 + 4 * w - 1)
    d = -565358149646
    assert (radicand / d).is_square()
    imaginary = QuadraticField(d, "a")
    mirror = QuadraticField(-3 * d, "b")
    assert imaginary.class_group(proof=True).invariants() == (66690, 2, 2, 2, 2)
    assert mirror.class_group(proof=True).invariants() == (6, 2, 2, 2, 2)
    Y = QQ(radicand / d).sqrt() * imaginary.gen()
    gamma = Y / (w - 1)**2 + (w + 3) / (w - 1)
    assert gamma.norm() == (-2 / (w - 1))**3
    assert not gamma.is_nth_power(3)
    assert GF(7)(2)**2 != 1
    assert residue(gamma[0]) - residue(gamma[1]) == residue(8)

    ring = PolynomialRing(QQ, "T")
    T = ring.gen()
    cardano = T**3 + 6 / (w - 1) * T - 2 * (w + 3) / (w - 1)
    monic = 2077**3 * cardano(T / 2077)
    assert monic == T**3 - 6411699 * T - 163929302
    assert monic.is_irreducible()
    assert pari(monic).nfdisc() == mirror.discriminant()


def check_varying_fixture():
    t = QQ(10)
    b = t**6 / 4 - 1
    s = 70 / (35 * t**4 + 22)
    w = b**2 * s**3
    assert b.valuation(3) == 1 and s.valuation(3) == -1
    assert w.valuation(2) == 0 and w.valuation(5) > 0 and w.valuation(7) > 0
    assert Integers(9)(1 / w) == Integers(9)(3)
    assert QQ(1) / 3 < w < QQ(1) / 2


if __name__ == "__main__":
    check_mirror_fixture()
    check_varying_fixture()
    print("Verified the cube condition at three, both class groups, and the fundamental reflected cubic discriminant.")

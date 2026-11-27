"""Replay the varying inverse parameter identities and a local specialization."""

from sage.all import GF, PolynomialRing, QQ, QuadraticField


def coefficient_from_branch(r):
    k = -(2 * r + 1) / (r**3 * (r + 2))
    return 4 * k * (k + 1)


def check_parameter_identities():
    ring = PolynomialRing(QQ, "t")
    field = ring.fraction_field()
    t = field(ring.gen())
    A = t**3 / 2
    B = (A**2 + 1) / (2 * A)
    b = A**2 - 1
    s = 70 / (35 * t**4 + 2)
    w = b**2 * s**3
    assert B**2 - 1 == b**2 / t**6
    assert w == (1 - 4 / t**6)**2 / (2 * (1 + 2 / (35 * t**4))**3)
    for Z in (A, B):
        r = 2 * Z / (3 - Z)
        scale = -(r + 1) / r**2
        assert coefficient_from_branch(r) == (Z**2 - 1) * scale**3


def check_specialization():
    t = QQ(2)
    b = t**6 / 4 - 1
    s = 70 / (35 * t**4 + 2)
    w = b**2 * s**3
    assert b == 15 and s == QQ(35) / 281
    assert w.valuation(2) == 0 and w.valuation(5) > 0 and w.valuation(7) > 0
    assert QQ(1) / 3 < w < QQ(1) / 2
    radicand = (w**2 - 1) * (w**2 + 4 * w - 1)
    numerator_part = radicand.numerator().squarefree_part()
    denominator_part = radicand.denominator().squarefree_part()
    d = (numerator_part * denominator_part).squarefree_part()
    field = QuadraticField(d, "r")
    r = field.gen()
    Y = QQ(radicand / d).sqrt() * r
    X = (w + 1) / (w - 1)
    U = Y / (w - 1)**2
    gamma = U + 2 * X - 1
    assert U**2 == X**3 + X**2 - X
    assert gamma.norm() == -(X - 1)**3
    assert not gamma.is_nth_power(3)
    assert GF(7)(-2)**2 != 1 and w.valuation(7) > 0


if __name__ == "__main__":
    check_parameter_identities()
    check_specialization()
    print("Verified varying inverse profile identities and the t=2 virtual product specialization.")

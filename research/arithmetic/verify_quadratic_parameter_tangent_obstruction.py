"""Replay the genus-three parameter identities and one exact specialization."""

from sage.all import PolynomialRing, QQ, QuadraticField


def check_parameter_identities():
    ring = PolynomialRing(QQ, "z")
    field = ring.fraction_field()
    z = field(ring.gen())

    def heron(variable):
        return 1 + 6 * variable - 3 * variable**2

    assert heron(2 - z) == heron(z)
    assert heron(-z / (6 * z + 1)) == heron(z) / (6 * z + 1)**2
    assert heron(2 + 1 / (3 * z)) == -heron(z) / (3 * z**2)
    assert heron(QQ(2)) == 1
    assert heron(-QQ(1) / 6) == -QQ(1) / 12


def check_specialization():
    field = QuadraticField(193, "r")
    r = field.gen()
    first_parameter = QQ(7) / 5
    second_parameter = QQ(1) / 5
    assert first_parameter**2 + second_parameter**2 == 2

    def edge_norm(parameter):
        a, b = parameter.numerator(), parameter.denominator()
        radicand = b**4 + 6 * a**2 * b**2 - 3 * a**4
        assert radicand == 4 * 193
        element = (2 * a**3 - 3 * a**2 * b - b**3 + 2 * b * r) / 2
        assert element.norm() == (a * (a - b))**3
        return element

    first = edge_norm(first_parameter)
    second = edge_norm(second_parameter)
    assert first == -87 + 5 * r and second == -69 + 5 * r
    vals = [(first.valuation(prime), second.valuation(prime)) for prime in field.primes_above(2)]
    assert sorted(vals) == [(1, 5), (2, 1)]
    assert all(not a.is_nth_power(3) for a in (first, second, first * second, first * second**2))
    assert field.class_number(proof=True) == 1


if __name__ == "__main__":
    check_parameter_identities()
    check_specialization()
    print("Verified the Möbius identities and the independent dyadic norm classes over Q(sqrt(193)).")

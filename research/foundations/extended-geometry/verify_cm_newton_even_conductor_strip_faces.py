"""Replay arithmetic certificates for the order-fourteen and order-eighteen faces."""

from math import gcd

from sage.all import (
    CyclotomicField,
    GF,
    PolynomialRing,
    QQ,
    euler_phi,
    lcm,
    matrix,
)


def relative_polynomial(prime, order, first, second):
    coefficient = CyclotomicField(order)
    zeta = coefficient.gen()
    relative = PolynomialRing(coefficient, "X")
    x = relative.gen()

    field = GF(prime**2, "u")
    generator = field.multiplicative_generator()
    logarithms = {}
    power = field.one()
    for index in range(prime**2 - 1):
        logarithms[power] = index
        power *= generator
    traces = []
    for elements, stride in ((map(field, range(prime)), prime+1), (iter(field), 1)):
        counts = [0] * order
        for z in elements:
            value = z**(order-1) * (z-1)**first * (z-3)**second
            if value:
                log = logarithms[value]
                assert log % stride == 0
                counts[(log // stride) % order] += 1
        traces.append(-sum(count * zeta**i for i, count in enumerate(counts)))
    trace, second_power_trace = traces
    determinant = (trace**2 - second_power_trace) / 2
    polynomial = relative.one()
    for i in primitive_characters(order):
        conjugate = coefficient.hom([zeta**i], coefficient)
        polynomial *= x**2 - conjugate(trace)*x + conjugate(determinant)
    return PolynomialRing(QQ, "X")(polynomial), trace**2 - 4*determinant


def primitive_characters(order):
    return tuple(i for i in range(1, order) if gcd(order, i) == 1)


def distinct_powers(polynomial, ell, exponent):
    reduced = polynomial.change_ring(GF(ell))
    ring = reduced.parent().quotient(reduced)
    root = ring.gen()
    power = root**exponent
    rows = [(power * root**j).list() for j in range(polynomial.degree())]
    powered = matrix(GF(ell), rows).charpoly()
    assert powered.gcd(powered.derivative()).degree() == 0


def verify_arithmetic():
    exponent = lcm([j for j in range(1, 41473) if euler_phi(j) <= 144])
    records = (
        (14, 3, 3, 43, 109, 43**4 * 1583**2),
        (14, 3, 3, 127, 109, 13**2 * 71**2 * 127**4),
        (18, 2, 7, 19, 103, 3**4 * 19**4 * 71**2),
        (18, 2, 7, 109, 103, 19**2 * 109**4 * 163**2),
    )
    for order, first, second, prime, ell, wanted_norm in records:
        polynomial, discriminant = relative_polynomial(prime, order, first, second)
        assert polynomial.degree() == 12 and polynomial.is_irreducible()
        assert discriminant.norm() == wanted_norm
        distinct_powers(polynomial, ell, exponent)
        if order == 14 and prime == 43:
            field = GF(1583)
            values = map(field, discriminant.list())
            expression = PolynomialRing(field, "T")(list(values))
            zeros = {int(value) for value in field
                     if value and value.multiplicative_order() == 14 and expression(value) == 0}
            assert zeros == {261, 279}
        print(f"n={order}, p={prime}: irreducible and distinct under all powers, "
              f"discriminant norm {wanted_norm}.")


def verify_path(order, first, second, p, expected_path):
    infinity = order - first - second + 1
    branches = (order - 1, first, second, infinity)
    assert sum(branches) == 2*order and min(branches) > 0 and max(branches) < order
    dimensions = {i: sum(i*e % order for e in branches)//order - 1
                  for i in primitive_characters(order)}
    assert [i for i, d in dimensions.items() if d == 1] == [1, order-1]
    path = [order-1]
    while path[-1] != 1:
        path.append(p*path[-1] % order)
        assert len(path) <= len(dimensions)
    assert tuple(path) == expected_path
    assert all(dimensions[i] > 0 for i in path)
    assert dimensions[p % order] == 0
    assert {pow(p, j, order) for j in range(len(dimensions))} == set(dimensions)

    r = len(path)-1
    assert ((order-1)*p**r-1) % order == 0
    m = ((order-1)*p**r-1)//order
    z = PolynomialRing(GF(p), "z").gen()
    selected = (((z-1)**first * (z-3)**second)**m)[m]
    assert selected != 0
    print(f"n={order}, p={p}, m={m}: path {tuple(path)}, selected entry {selected}.")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    verify_arithmetic()
    verify_path(14, 3, 3, 5, (13, 9, 3, 1))
    verify_path(18, 2, 7, 29, (17, 7, 5, 1))


if __name__ == "__main__":
    main()

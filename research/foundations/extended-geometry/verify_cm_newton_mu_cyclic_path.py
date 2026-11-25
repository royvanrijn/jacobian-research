"""Replay the fixed degree-five cyclic-cover calculation with SageMath."""

from math import comb

from sage.all import (CyclotomicField, GF, PolynomialRing, QQ, euler_phi,
                      lcm, matrix)


def relative_frobenius(p, zeta):
    """Recover relative trace and determinant from two character sums."""
    field = GF(p**2, "u")
    generator = field.multiplicative_generator()
    logs = {}
    value = field.one()
    for index in range(p**2 - 1):
        logs[value] = index
        value *= generator

    traces = []
    for elements, stride in ((map(field, range(p)), p + 1), (iter(field), 1)):
        counts = [0] * 5
        for z in elements:
            value = (z * (z - 1) * (z - 3))**2
            if value:
                index = logs[value]
                assert index % stride == 0
                counts[(index // stride) % 5] += 1
        traces.append(-sum(count * zeta**i for i, count in enumerate(counts)))
    trace, second_power_trace = traces
    return trace, (trace**2 - second_power_trace) / 2


def verify_distinct_powers(polynomial, ell, exponent):
    reduced = polynomial.change_ring(GF(ell))
    ring = reduced.parent().quotient(reduced)
    root = ring.gen()
    power = root**exponent
    values = [(power * root**j).list() for j in range(8)]
    powered_polynomial = matrix(GF(ell), values).charpoly()
    assert powered_polynomial.gcd(powered_polynomial.derivative()).degree() == 0
    return powered_polynomial


def main():
    field = CyclotomicField(5)
    zeta = field.gen()
    relative = PolynomialRing(field, "x")
    x = relative.gen()
    rationals = PolynomialRing(QQ, "x")
    exponent = lcm([j for j in range(1, 8193) if euler_phi(j) <= 64])
    records = (
        (11, 103, -zeta**2 - 2*zeta, 7*zeta**3 + 9*zeta**2 - 3*zeta + 3,
         5**2 * 11**2 * 29**2),
        (61, 137, 3*zeta**3 + 2*zeta**2 + 2, 24*zeta**3 + 16*zeta**2 + 39*zeta + 72,
         5**2 * 61**2 * 181**2),
    )
    for p, ell, expected_trace, expected_determinant, expected_norm in records:
        trace, determinant = relative_frobenius(p, zeta)
        assert (trace, determinant) == (expected_trace, expected_determinant)
        assert (trace**2 - 4*determinant).norm() == expected_norm
        polynomial = relative.one()
        for i in range(1, 5):
            conjugate = field.hom([zeta**i], field)
            polynomial *= x**2 - conjugate(trace)*x + conjugate(determinant)
        polynomial = rationals(polynomial)
        assert polynomial.is_irreducible()
        distinctness = verify_distinct_powers(polynomial, ell, exponent)
        print(f"p={p}: irreducible degree eight, norm discriminant {expected_norm}, "
              f"separable power polynomial modulo {ell}: {distinctness}")

    p = 13
    k = (4*p**2 - 1) // 5
    entry = (-1)**k * sum(comb(2*k, i) * comb(2*k, 3*k-i) * pow(3, i-k, p)
                         for i in range(k, 2*k + 1)) % p
    assert entry == 11
    print(f"The selected directed entry at p={p}, k={k}, is {entry}.")


if __name__ == "__main__":
    main()

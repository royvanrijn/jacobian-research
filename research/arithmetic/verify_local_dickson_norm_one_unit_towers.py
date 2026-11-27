"""Exact unit representatives and discriminant fixtures for LD1–LD2."""

from sage.all import GF, QQ, NumberField, PolynomialRing, lcm, pari, valuation


ring = PolynomialRing(QQ, "x")
x = ring.gen()


def verify_representatives():
    unramified = NumberField(x**2 + 1, "i")
    i = unramified.gen()
    unramified_unit = (1 + 3*i)/(1 - 3*i)
    assert unramified_unit.norm() == 1
    assert unramified_unit.trace() == -QQ(8)/5

    noncyclotomic = NumberField(x**2 - 3, "r")
    ramified_unit = -2 - noncyclotomic.gen()
    assert ramified_unit.norm() == 1
    assert ramified_unit.trace() == -4
    assert valuation(ramified_unit.trace() + 2, 3) == 0
    assert valuation(ramified_unit.trace() - 2, 3) == 1

    cyclotomic = NumberField(x**2 + 3, "s")
    pi = cyclotomic.gen()
    zeta = (-1 + pi)/2
    eta = (1 + pi**3)/(1 - pi**3)
    assert zeta**3 == 1 and zeta != 1
    assert eta.norm() == 1 and eta.trace() == -QQ(13)/7
    assert (zeta*eta).trace() == QQ(11)/7
    assert (zeta**2*eta).trace() == QQ(2)/7
    cubic = x**3 + 2*x**2 + x + QQ(1)/7
    assert cubic.change_ring(GF(3)).is_irreducible()


def dickson_polynomials(maximum_degree):
    polynomials = [ring(2), x]
    for degree in range(2, maximum_degree + 1):
        polynomials.append(x*polynomials[degree - 1] - polynomials[degree - 2])
    return polynomials


def field_discriminant_valuation(polynomial):
    denominator = lcm([coefficient.denominator() for coefficient in polynomial])
    integral_polynomial = denominator**polynomial.degree()*polynomial(x/denominator)
    return valuation(pari(integral_polynomial).nfdisc(), 3)


def verify_discriminants_and_conductors(dickson):
    for n in range(1, 4):
        m = 3**n
        split_or_unramified = n*m
        depth_one_ramified = n*m + (m - 1)//2
        pure_cyclotomic = ((2*n - 1)*m - 3)//2
        examples = [
            (QQ(17)/4, split_or_unramified),
            (-QQ(8)/5, split_or_unramified),
            (QQ(4), depth_one_ramified),
            (-QQ(13)/7, pure_cyclotomic),
            (QQ(11)/7, depth_one_ramified),
            (QQ(2)/7, depth_one_ramified),
        ]
        for trace, expected in examples:
            polynomial = dickson[m] - trace
            assert polynomial.is_irreducible()
            assert field_discriminant_valuation(polynomial) == expected

        full_base = (2*n - 1)*(m//3)
        real_base = (full_base - 1)//2
        assert (2*n + 1)*(m//3) - full_base == 2*(m//3)
        assert (2*n + 2)*(m//3) - 2*real_base == m + 1
        if n > 1:
            assert n*(m//3) - real_base == (m//3 + 1)//2


def verify_imprimitive_types(dickson):
    cyclotomic = NumberField(x**2 + 3, "p")
    pi = cyclotomic.gen()
    zeta = (-1 + pi)/2
    eta = (1 + pi**3)/(1 - pi**3)

    for n in (2, 3):
        m = 3**n
        for s in range(1, n + 1):
            r = n - s
            for j in (0, 1):
                gamma = zeta**j * eta**(3**s)
                polynomial = dickson[m] - gamma.trace()
                if j == 0:
                    central = ((2*r - 1)*3**r - 3)//2 if r else 0
                    expected = [(3**r, central)]
                    for degree_exponent in range(r, n):
                        degree = 3**degree_exponent
                        discriminant = ((2*degree_exponent + 1)*degree - 1)//2
                        expected.extend([(degree, discriminant)] * 2)
                else:
                    expected = [(m, ((2*n + 1)*m - 1)//2)]
                actual = [
                    (factor.degree(), field_discriminant_valuation(factor))
                    for factor, multiplicity in polynomial.factor()
                    for _ in range(multiplicity)
                ]
                assert sorted(actual) == sorted(expected)
                local_factors = pari(polynomial).factorpadic(3, 32)
                assert sorted(int(factor.poldegree()) for factor in local_factors[0]) == sorted(
                    degree for degree, _ in expected
                )

            if r:
                multiplier = 4
                translation_stride = 3**s
                subgroup_order = 3**(n - 1)
                for lift in range(3**r):
                    translation = 1 + lift*translation_stride
                    part = translation*((multiplier**subgroup_order - 1)//3)
                    assert valuation(part, 3) == n - 1
                    orbit = {0}
                    value = translation
                    while value not in orbit:
                        orbit.add(value)
                        value = (multiplier*value + translation) % m
                    assert len(orbit) == m


if __name__ == "__main__":
    dickson = dickson_polynomials(27)
    verify_representatives()
    verify_discriminants_and_conductors(dickson)
    verify_imprimitive_types(dickson)
    print("LD1–LD3 unit, discriminant, conductor, factor and affine action fixtures verified.")

"""Check the resultant and finite-field certificates for the shifted windows."""

from argparse import ArgumentParser

from sage.all import GF, PolynomialRing, QQ, lcm, matrix
from sympy import Poly, degree, isprime, resultant, symbols


r, x, t = symbols("r x t")
d = -373248 * r**5 * (12*r + 1)**2
a = 161243136*r**5 + 42456960*r**4 + 3203712*r**3 + 51624*r**2 - 612*r + 1
b = -216 * (559872*r**3 + 201528*r**2 + 23412*r + 881)
c = -3888 * (5832*r**2 + 540*r - 1)
u = a**2 - 12*b*d
j = c*d**2
k = a**3


def certificate(p):
    reduced_u = Poly(u, r, modulus=p)
    j_reduced = Poly(j, r, modulus=p).rem(reduced_u)
    k_reduced = Poly(k, r, modulus=p).rem(reduced_u)
    reduced_v = sum(
        coeff * j_reduced**power * k_reduced**(10 - power)
        for (power,), coeff in reduced_u.terms()
    ).rem(reduced_u)
    return reduced_u, reduced_v, reduced_u.gcd(reduced_v)


def multiplication_determinant(cubic, numerator, denominator):
    cubic_ring = cubic.parent()
    coefficients = cubic_ring.base_ring()
    x = cubic_ring.gen()

    def multiplication_matrix(multiplier):
        columns = [(multiplier*x**i) % cubic for i in range(3)]
        return matrix(coefficients, 3, 3, lambda i, j: columns[j][i])

    denominator_matrix = multiplication_matrix(denominator)
    numerator_matrix = multiplication_matrix(numerator)
    v0, v1, v2, v3 = [(i*denominator_matrix - numerator_matrix).det()
                       for i in range(4)]
    d = (v3 - 3*v2 + 3*v1 - v0) / 6
    a = (v2 - 2*v1 + v0) / 2 - 3*d
    b = v1 - v0 - d - a
    return d, a, b, v0


def second_window_polynomials():
    coefficients = PolynomialRing(QQ, "h")
    h = coefficients.gen()
    cubic_ring = PolynomialRing(coefficients, "x")
    x = cubic_ring.gen()
    product = (240*h*h - 40*h + 3) / 128
    cubic = x**3 - x*x + h*x - product
    d, a, b, c = multiplication_determinant(
        cubic, 1 - 4*h + 2*x - 3*x*x, x**7 * (1 - 2*x)**2,
    )
    u = 3*a**4 - 40*b*d*a*a + 240*(b*d)**2 - 128*c*d*d*a
    assert u.degree() == 53
    return u, b*d, a*a


def third_shift_at_seven():
    sqsum, total, product = symbols("h s r")
    delta = 4*sqsum - total**2
    third_relation = (5*delta**3 + 24*sqsum*delta**2
                      + 128*(sqsum**2 - total*product)*delta
                      + 1024*(sqsum**3 - 2*total*sqsum*product + product**2))
    assert Poly(third_relation.subs(total, 1), sqsum, product, modulus=7) == (
        Poly(2*(sqsum*product + product**2 + product + 1), sqsum, product, modulus=7)
    )
    assert Poly(third_relation.subs(total, 0), sqsum, product, modulus=7) == (
        Poly(2*product**2, sqsum, product, modulus=7)
    )
    polynomials = PolynomialRing(GF(7), "r")
    r = polynomials.gen()
    rational_functions = polynomials.fraction_field()
    h = -r - 1 - 1/rational_functions(r)
    cubic_ring = PolynomialRing(rational_functions, "x")
    x = cubic_ring.gen()
    cubic = x**3 - x*x + h*x - r
    d, a, b, c = multiplication_determinant(
        cubic, 1 - 4*h + 2*x - 3*x*x, x**9 * (1 - 2*x)**2,
    )
    next_h = b*d/a**2
    next_r = c*d*d/a**3
    required = (next_h*next_r + next_r**2 + next_r + 1).numerator().monic()
    discriminant = h*h - 4*h**3 - 4*r - 27*r*r + 18*h*r
    excluded = (r * discriminant.numerator() * a.numerator() * d.numerator()
                * next_r.denominator())
    while (factor := required.gcd(excluded)).degree() > 0:
        required = required // factor
    assert required.degree() == 144

    numerator = next_r.numerator() % required
    denominator = next_r.denominator() % required
    iterate = polynomials.zero()
    denominator_power = polynomials.one()
    for coefficient in reversed(required.list()):
        iterate = (iterate*numerator + coefficient*denominator_power) % required
        denominator_power = denominator_power*denominator % required
    assert required.gcd(iterate).is_one()
    assert pow(7, 60, 29**2) == 226
    assert pow(7, 60, 25*9) == 1


def second_certificate(p, polynomials):
    field = PolynomialRing(GF(p), "h")
    u, b, a = map(field, polynomials)
    assert not u.is_zero()
    a %= u
    b %= u
    # Use the unreduced degree of U in the homogenized substitution.
    v = sum(
        coeff * (b**power % u) * (a**(53 - power) % u)
        for power, coeff in enumerate(u.list())
    ) % u
    return u, u.gcd(v)


def parameter_field(p, factor):
    degree = factor.degree()
    if degree == 1:
        return GF(p), -factor[0]
    field = GF(p**degree, name="z", modulus=factor)
    return field, field.gen()


def cubic_frobenius_periods(p, field, degree, h, r, shift, total=1):
    ring = PolynomialRing(field, "x")
    x = ring.gen()
    cubic = x**3 - total*x*x + h*x - r
    denominator = x**(2*shift+3) * (total-2*x)**2
    if (not cubic.gcd(cubic.derivative()).is_one()
            or not cubic.gcd(denominator).is_one()):
        return []
    value = ((total*total-4*h+2*total*x-3*x*x)
             * denominator.inverse_mod(cubic)) % cubic
    if total:
        trace = sum(((value*x**i) % cubic)[i] for i in range(3))
        if not trace:
            return []
        expected = value * (1/trace)
    else:
        if not value.gcd(cubic).is_one():
            return []
        inverse = value.inverse_mod(cubic)
    relative_period = lcm([factor.degree() for factor, _ in cubic.factor()])
    period = degree * relative_period
    image = x
    periods = []
    for exponent in range(1, period+1):
        image = image**p % cubic
        if total:
            compatible = image == expected
        else:
            compatible = ((image*inverse) % cubic).degree() == 0
        if compatible:
            periods.append((period, exponent % period))
    return periods


def root_frobenius_periods(p, parameter_polynomial, shift):
    """Find exponent residues allowed by compatibility on every cubic root."""
    periods = []
    for factor, _ in parameter_polynomial.factor():
        field, parameter = parameter_field(p, factor)
        if shift == 1:
            h, r = 1/field(12), parameter
        else:
            h = parameter
            r = (240*h*h - 40*h + 3)/128
        periods.extend(cubic_frobenius_periods(p, field, factor.degree(), h, r, shift))
    return sorted(set(periods))


def first_parameter_polynomial(p):
    _, _, common = certificate(p)
    ring = PolynomialRing(GF(p), "r")
    return ring([int(value) for value in reversed(common.all_coeffs())])


def fixed_shift_certificate(p, shift):
    field = GF(p)
    ring = PolynomialRing(field, "h,r")
    h, r = ring.gens()
    cubic_ring = PolynomialRing(ring, "x")
    x = cubic_ring.gen()

    def relation(total, pairs, product):
        delta = 4*pairs - total*total
        numerators = [1, pairs]
        for _ in range(shift-1):
            third = numerators[-3] if len(numerators) >= 3 else 0
            numerators.append(pairs*numerators[-1] - total*product*numerators[-2]
                              + product*product*third)
        coefficient = field(1)
        result = ring.zero()
        for i, value in enumerate(numerators):
            result += coefficient*delta**(shift-i)*value
            if i < shift:
                coefficient *= field(8*(shift-i))/field(2*(shift-i)-1)
        return result

    def determinant(total, pairs, product):
        cubic = x**3 - total*x*x + pairs*x - product
        numerator = total*total - 4*pairs + 2*total*x - 3*x*x
        return multiplication_determinant(cubic, numerator, x**(2*shift+3) * (total-2*x)**2)

    # The s=0 chart permits r=1, and h=0 is the known cube-root point.
    d, a, _, c = determinant(0, h, 1)
    zero_sum = ring.ideal([relation(0, h, 1), a])
    excluded = h*(-4*h**3-27)*d*c
    admissible_zero_sum = zero_sum.saturation(ring.ideal(excluded))[0]
    zero_sum_basis = admissible_zero_sum.groebner_basis()
    assert len(zero_sum_basis) == 1 and zero_sum_basis[0].degree(r) == 0
    if shift == 4 or p == 11:
        assert admissible_zero_sum.is_one()
    elif shift == 5 and p == 13:
        assert zero_sum_basis[0] == h**3 - 6
    zero_sum_periods = []
    for factor, _ in zero_sum_basis[0].univariate_polynomial().factor():
        base, pairs = parameter_field(p, factor)
        zero_sum_periods.extend(cubic_frobenius_periods(
            p, base, factor.degree(), pairs, base.one(), shift, total=0,
        ))

    d, a, b, c = determinant(1, h, r)
    required = relation(1, h, r)
    image = relation(a, b*d, c*d*d)
    discriminant = h*h-4*h**3-4*r-27*r*r+18*h*r
    both = ring.ideal([required, image]).saturation(
        ring.ideal(r*discriminant*d*a)
    )[0]
    lex_ring = PolynomialRing(field, "h,r", order="lex")
    basis = both.transformed_basis("fglm", other_ring=lex_ring)
    eliminants = [value for value in basis if value.degree(lex_ring.gen(0)) == 0]
    assert len(eliminants) == 1
    parameter_polynomial = eliminants[0].univariate_polynomial().monic()
    periods = []
    for factor, _ in parameter_polynomial.factor():
        base, product = parameter_field(p, factor)
        local_ring = PolynomialRing(base, "h")
        z = local_ring.gen()

        def specialize(value):
            return sum(coeff*z**i*product**j for (i, j), coeff in value.dict().items())

        fiber = specialize(required).gcd(specialize(image))
        excluded = specialize(d*a*discriminant)
        while (cut := fiber.gcd(excluded)).degree() > 0:
            fiber = fiber // cut
        if fiber.degree() == 0:
            continue

        candidate_image = (specialize(c*d*d) *
                           specialize(a**3).inverse_mod(fiber)) % fiber
        same_conjugacy_class = local_ring.zero()
        for coefficient in reversed(factor.list()):
            same_conjugacy_class = (same_conjugacy_class*candidate_image
                                    + coefficient) % fiber
        for h_factor, _ in fiber.gcd(same_conjugacy_class).factor():
            assert h_factor.degree() == 1
            periods.extend(cubic_frobenius_periods(
                p, base, factor.degree(), -h_factor[0], product, shift,
            ))
    return parameter_polynomial.degree(), sorted(set(periods + zero_sum_periods))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("primes", type=int, nargs="*")
    args = parser.parse_args()
    for p in args.primes:
        if p < 5 or not isprime(p):
            parser.error(f"{p} must be a prime at least 5")

    actual = Poly(resultant(
        12*x**3 - 12*x**2 + x - 12*r,
        3*t*x**5*(1 - 2*x)**2 - (2 + 6*x - 9*x**2), x
    ), t)
    expected = Poly(1152 * (d*t**3 + a*t**2 + b*t + c), t)
    assert actual == expected
    assert degree(u, r) == 10
    for p in (5, 7, 17):
        reduced_u, reduced_v, common = certificate(p)
        assert common.is_one
        if p == 5:
            assert reduced_u.degree() == 10
            assert reduced_v == Poly(-r**9 + r**7 + 2*r**6 + r**4 - r**2 + 2,
                                     r, modulus=5)

    second_polynomials = second_window_polynomials()
    for p in (19, 29, 37, 41):
        reduced_u, common = second_certificate(p, second_polynomials)
        assert common.is_one()
        assert reduced_u.degree() == 53
    print("First shift: resultant and coprimality at 5, 7 and 17 pass.")
    print("Second shift: degree 53 and coprimality at 19, 29, 37 and 41 pass.")
    third_shift_at_seven()
    print("Third shift at 7: degree 144 and Frobenius compatibility gcd 1 pass.")
    first_all = (11, 29, 59, 83, 359)
    second_all = (13, 23, 43, 53, 79, 89, 293, 313)
    for p in first_all:
        common = first_parameter_polynomial(p)
        assert common.degree() > 0 and not root_frobenius_periods(p, common, 1)
    for p in second_all:
        _, common = second_certificate(p, second_polynomials)
        assert common.degree() > 0 and not root_frobenius_periods(p, common, 2)
    assert root_frobenius_periods(13, first_parameter_polynomial(13), 1) == [(3, 1)]
    assert root_frobenius_periods(281, first_parameter_polynomial(281), 1) == [(2, 1)]
    for p, residues in [(7, [(3, 0)]), (17, [(3, 1), (6, 5)]),
                        (73, [(3, 0)]), (211, [(2, 0)])]:
        _, common = second_certificate(p, second_polynomials)
        assert root_frobenius_periods(p, common, 2) == residues
    assert pow(13, 60, 25*9) == 1
    assert pow(13, 60, 7) == 1
    print("Individual cubic roots: additional prime and exponent filters pass.")
    for p, expected in [(11, (298, [(3, 0)])), (13, (293, [(2, 0)])),
                        (17, (298, []))]:
        assert fixed_shift_certificate(p, 4) == expected
    assert pow(17, 1260, 27*49*25) == 1
    print("Fourth shift: certificates at 11, 13 and 17 pass.")
    assert fixed_shift_certificate(11, 5) == (533, [])
    assert fixed_shift_certificate(13, 5) == (532, [(4, 0)])
    assert pow(11, 10080, 23**2) == 404
    assert pow(11, 1260, 27*49*25) == 1
    print("Fifth shift: certificates at 11 and 13 pass.")
    for p in args.primes:
        _, _, common = certificate(p)
        print(f"p={p}, first shift: gcd={common.as_expr()}; candidate periods "
              f"{root_frobenius_periods(p, first_parameter_polynomial(p), 1)}")
        if p >= 7 and p != 11:
            _, common = second_certificate(p, second_polynomials)
            print(f"p={p}, second shift: gcd={common}; candidate periods "
                  f"{root_frobenius_periods(p, common, 2)}")


if __name__ == "__main__":
    main()

"""Exact polynomial, local, independence, and reduction data for CP1–CP4."""

from sage.all import CRT_list, EllipticCurve, GF, QQ, ZZ, PolynomialRing, QuadraticField, matrix


ring = PolynomialRing(QQ, "t")
t = ring.gen()
P = (
    QQ(55296)/1153**2*t**4
    - QQ(2670347)/(4*1153**2)*t**2
    + QQ(1037316177403)/(2**17*1153**2)
)
Z = t**2/1153 - t/64 + QQ(1378561)/(2**14*1153)
H = 48*t**2/1153 - t/4 - QQ(591617)/(2**10*1153)
F = P*P + 1
conditions = {
    2**12: 1105,
    3**5: 1,
    5: 2,
    13: 6,
    17: 1,
    19: 2,
    29: 2,
    31: 1,
    1153**2: 657,
    2677: 1,
}


def verify_identities():
    parameter_ring = PolynomialRing(QQ, "v")
    v = parameter_ring.gen()
    formula_ring = PolynomialRing(parameter_ring.fraction_field(), "u")
    u = formula_ring.gen()
    s = u - 3/(8*v)
    z_squared = (v - 24*v*v)/48
    numerator = (
        s + 3 + 6*v*(s+1)**2 + 12*v*v*(s+1)*(s*s+3)
        + 8*v**3*(s*s+3)**2 - 8*z_squared*(s*s+3)
    )
    assert numerator == numerator(-u)
    w = QQ(1)/48
    v = 1/(24*(1+2*w*w))
    z = w*v
    s = t - 3/(8*v)
    assert Z == z*(s*s+3)
    assert H == (s+1)/2 + v*(s*s+3)
    assert P == P(-t)
    assert P.degree() == 4 and F.degree() == 8 and F.is_squarefree()
    assert Z.degree() == H.degree() == 2
    assert 2*P*Z+Z*Z-1 == H**3
    assert H.resultant(Z*Z+1) == QQ(5*13*17*29*2677)/(2**4*1153**4)


def verify_integral_slice():
    moduli = [ZZ(m) for m in conditions]
    residues = [ZZ(conditions[m]) for m in conditions]
    moment = CRT_list(residues, moduli)
    period = ZZ(1)
    for modulus in moduli:
        period *= modulus
    assert P((moment+period*t)/128).denominator() == 1
    assert [a.valuation(2) for a in P((1105+2**12*t)/128).list()] == [0, 1, 9, 21, 31]
    p = ZZ(1153)
    argument = (657+p*p*t)/128
    for polynomial in [P, Z, Z(-t), H, H(-t)]:
        assert polynomial(argument).denominator() % p != 0
    assert H(argument)[0] % p == 792 and H(-argument)[0] % p == 718
    for q in (5, 17, 29, 2677):
        argument = QQ(conditions[q])/128
        assert H(argument) % q != 0 and H(-argument) % q != 0
    argument = QQ(conditions[13])/128
    assert H(argument) % 13 == 0 and H(-argument) % 13 == 5
    assert (Z(argument)**2+1) % 13 == 11


def verify_three_and_cubic_residues():
    modulus = 3**5
    argument = QQ(conditions[modulus])/128
    p = ZZ(P(argument) % modulus)
    r = 100
    assert p == 87 and r % 3 == 1 and (r*r-p*p-1) % modulus == 0
    unit = p+r
    assert unit % 9 == 7
    first = (1-Z(argument)*unit) % modulus
    second = (1-Z(-argument)*unit) % modulus
    assert first == 223 and second == 27
    assert first*unit**2 % 9 == 1 and (second/27) % 9 == 1

    rows = []
    for q, root, omega, exponents in [
        (13, 11, 3, [2, 1, 0]),
        (19, 5, 7, [2, 2, 0]),
        (31, 6, 25, [0, 1, 1]),
    ]:
        residue_field = GF(q)
        argument = QQ(conditions[q])/128
        p = residue_field(P(argument))
        r = residue_field(root)
        assert r != 0 and r*r == p*p+1
        omega = residue_field(omega)
        assert omega != 1 and omega**3 == 1
        unit = p+r
        values = [unit, 1-residue_field(Z(argument))*unit, 1-residue_field(Z(-argument))*unit]
        for value, exponent in zip(values, exponents):
            assert value != 0 and value**((q-1)//3) == omega**exponent
        rows.append(exponents)
    assert matrix(GF(3), rows).det() == 2


def verify_good_reduction():
    trace_sums = []
    for degree, expected_count in [(1, 12), (2, 140), (3, 1332)]:
        field = GF(11**degree, "a")
        polynomial = F.change_ring(field)
        assert polynomial.degree() == 8 and polynomial.is_squarefree()
        assert polynomial.leading_coefficient().is_square()
        count = 2
        for x in field:
            value = polynomial(x)
            if value == 0:
                count += 1
            elif value.is_square():
                count += 2
        assert count == expected_count
        trace_sums.append(11**degree+1-count)

    first = -trace_sums[0]
    second = -(trace_sums[1]+first*trace_sums[0])//2
    third = -(trace_sums[2]+first*trace_sums[1]+second*trace_sums[0])//3
    assert (first, second, third) == (0, 9, 0)
    assert ZZ(1+second+11*second+11**3).valuation(3) == 2


def verify_local_affine_orbits():
    for dimension in (1, 2, 3):
        pure_exponent = (3**dimension-1)//2
        mixed_exponent = 3 + 7*(3**(dimension-1)-1)//2
        assert mixed_exponent-pure_exponent == 2*3**(dimension-1)
    assert tuple((3**dimension-1)//2 for dimension in (2, 3)) == (4, 13)
    assert 3 + 2*2 == 7


def verify_family_squareclasses():
    gaussian = QuadraticField(-1, "i")
    i = gaussian.gen()
    parameter_ring = PolynomialRing(gaussian, "w")
    w = parameter_ring.gen()
    fraction = parameter_ring.fraction_field()
    quartic_ring = PolynomialRing(fraction, "t")
    u = quartic_ring.gen()
    v = 1/(24*(1+2*w*w))
    z = w*v
    s = u - 3/(8*v)
    quartic = (
        s+3+6*v*(s+1)**2+12*v*v*(s+1)*(s*s+3)
        + 8*v**3*(s*s+3)**2-8*z*z*(s*s+3)
    )/(16*z)
    assert quartic == quartic(-u)
    a, b, c = quartic[4], quartic[2], quartic[0]
    assert a == v*v/(2*w)
    constant = (c-i)/a
    discriminant = (b*b-4*a*(c-i))/(a*a)
    assert constant == -19440*(w*w+i*w-QQ(1)/5)*(w*w-i*w/3+QQ(1)/3)**3
    assert discriminant == 147456*(w-3*i/4)*(w+i/4)**3*(w*w+QQ(1)/2)**2
    first = w**4+QQ(3)/5*w*w+QQ(1)/25
    second = w**4+QQ(7)/9*w*w+QQ(1)/9
    assert (w*w+i*w-QQ(1)/5)*(w*w-i*w-QQ(1)/5) == first
    assert (w*w-i*w/3+QQ(1)/3)*(w*w+i*w/3+QQ(1)/3) == second
    assert (w-3*i/4)*(w+3*i/4) == w*w+QQ(9)/16
    assert (w+i/4)*(w-i/4) == w*w+QQ(1)/16

    C = gaussian(constant(QQ(1)/48))
    D = gaussian(discriminant(QQ(1)/48))
    values = [C, D, C.conjugate(), D.conjugate()]
    places = [
        gaussian.ideal(p, i-r)
        for p, r in [(101, 10), (101, -10), (1297, 36), (1297, -36)]
    ]
    rows = [[gaussian.ideal(value).valuation(place) % 2 for value in values] for place in places]
    assert matrix(GF(2), rows).rank() == 4


def verify_rank_zero_descents():
    square_mod_eight = {0, 1, 4}
    for u in range(8):
        for v in range(8):
            if (u | v) & 1:
                assert (-u**4+15*u*u*v*v-25*v**4) % 8 not in square_mod_eight
            if u & 1:
                assert (2*u**4-20*u*u*v*v+32*v**4) % 8 not in square_mod_eight
            if u & 1 and v & 1:
                assert (2*u**4-5*u*u*v*v+2*v**4) % 8 == 7
    for A, B, primes, order, points in [
        (10, 9, (5, 7), 8, [(0, 0), (-1, 0), (-9, 0), (3, 12), (3, -12), (-3, 6), (-3, -6)]),
        (15, 25, (3, 7), 4, [(0, 0), (5, 25), (5, -25)]),
    ]:
        curve = EllipticCurve([0, A, 0, B, 0])
        for prime in primes:
            assert curve.discriminant() % prime != 0
            assert curve.change_ring(GF(prime)).cardinality() == order
        assert len(set(points)) == order-1
        assert all(y*y == x*(x*x+A*x+B) for x, y in points)


def verify_squarefree_slice():
    refined = {25 if modulus == 5 else modulus: residue for modulus, residue in conditions.items()}
    moduli = [ZZ(m) for m in refined]
    start = CRT_list(list(refined.values()), moduli)
    period = ZZ(1)
    for modulus in moduli:
        period *= modulus
    integral_quartic = P((start+period*t)/128)
    assert integral_quartic.denominator() == 1
    assert integral_quartic.change_ring(GF(2)) == 1
    polynomial = integral_quartic**2+1
    for prime in period.prime_divisors():
        if prime == 2:
            continue
        modulus = 25 if prime == 5 else prime
        if prime == 5:
            assert polynomial[0] % modulus == 15
        else:
            assert polynomial[0] % modulus != 0
        assert all(coefficient % modulus == 0 for coefficient in polynomial.list()[1:])
    assert set(polynomial.leading_coefficient().numerator().prime_divisors()) <= set(period.prime_divisors())
    assert {k: sum(m*n for r, m in [(0, 5), (2, 2), (4, 1)]
                   for s, n in [(0, 5), (2, 2), (4, 1)] if r+s == k)
            for k in range(0, 9, 2)} == {0: 25, 2: 20, 4: 14, 6: 4, 8: 1}


if __name__ == "__main__":
    verify_identities()
    verify_integral_slice()
    verify_three_and_cubic_residues()
    verify_good_reduction()
    verify_local_affine_orbits()
    verify_family_squareclasses()
    verify_rank_zero_descents()
    verify_squarefree_slice()
    print("CP1–CP4 exact certificates passed")

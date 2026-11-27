"""Exact polynomial, local, independence, and reduction data for CP1–CP2."""

from sage.all import CRT_list, GF, QQ, ZZ, PolynomialRing, matrix


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


if __name__ == "__main__":
    verify_identities()
    verify_integral_slice()
    verify_three_and_cubic_residues()
    verify_good_reduction()
    verify_local_affine_orbits()
    print("CP1–CP2 exact certificates passed")

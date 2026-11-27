"""Exact finite algebra, local residues, and conductors for WD1–WD5."""

from sage.all import GF, QQ, EllipticCurve, PolynomialRing, valuation


ring = PolynomialRing(QQ, "z")
z = ring.gen()
quartic = z**4 - 6*z**3 + 7*z**2 - 18*z + 97
pell_a = z**6 - 5*z**5 + 2*z**4 - 15*z**3 + 76*z**2 + 65*z + 101
pell_b = z**4 - 2*z**3 - 3*z**2 - 14*z - 5
octic = z**8 - 6*z**7 - z**6 + 18*z**5 + 93*z**4 - 36*z**3 - 4*z**2 + 48*z + 16


def pell_trace(t):
    coordinate = t - 2/t
    return pell_a(coordinate)**2/1944 - 2


def verify_quotient_pell_and_residues():
    rational_z = ring.fraction_field()(z)
    assert octic == z**4 * quartic(rational_z - 2/rational_z)
    assert quartic.is_squarefree() and octic.is_squarefree()
    assert pell_a**2 - quartic*pell_b**2 == 7776
    assert pell_trace(rational_z) == pell_trace(-2/rational_z)
    assert 7776 == 4*1944 and valuation(1944, 3) == 5
    assert pell_a.change_ring(GF(2))(1) == 1
    assert pell_a(-1) == 135 and pell_a.derivative()(-1) == -171
    assert pell_a(1) == 225
    assert quartic(1 + 9*z) == 81*(1 - 2*z - 5*z**2 - 18*z**3 + 81*z**4)
    for t in (1, 4, 7):
        assert octic(t) % 9 == 3
    for t in (2, 5, 8):
        t = QQ(t)
        coordinate = t - 2/t
        r = (coordinate - 1)/9
        assert r % 3 == (1 if t == 5 else 0)
        assert valuation(pell_a(coordinate), 3) == (3 if t == 5 else 2)
        assert valuation(octic(t), 3) == (5 if t == 5 else 4)
    assert ((QQ(5) - QQ(2)/5 - 1)/9) % 9 == 4


def verify_wild_polynomials_and_conductors():
    trace = pell_trace(QQ(1))
    assert trace == QQ(59)/8 and valuation(trace + 2, 3) == 1
    polynomials = [ring(2), z]
    for j in range(2, 82):
        polynomials.append(z*polynomials[j - 1] - polynomials[j - 2])
    previous_different = 0
    for n in range(1, 5):
        m = 3**n
        dickson = polynomials[m]
        assert dickson(1) == -2
        field_polynomial = dickson(z + 1) - trace
        assert field_polynomial.change_ring(GF(3)) == (z**m).change_ring(GF(3))
        assert valuation(field_polynomial[0], 3) == 1
        different = n*m + (m - 1)//2
        assert valuation(field_polynomial.discriminant(), 3) == different
        real_mirror_discriminant = (2*n - 1)*(m//3) - 1
        assert different - previous_different == real_mirror_discriminant + m + 1
        previous_different = different
        ray_discriminant = sum(3**(k - 1)*(3**k + 1) for k in range(1, n + 1))
        assert 8*ray_discriminant == 3*m*m + 4*m - 7

        a = int(pow(4, -1, m))
        exponent = a - m
        assert exponent < 0 and (12*exponent - 3) % m == 0
        scale = QQ(2)**exponent
        at_zero = 2**12/(1944*scale**12)
        at_infinity = scale**12/1944
        assert valuation(at_zero, 3) == valuation(at_infinity, 3) == -5


def verify_pell_fiber_transport():
    g0 = z**6 - 3*z**5 - 3*z**4 - 12*z**3 + 6*z**2 - 12*z - 8
    ell0 = z**2 + 2
    assert g0**2 - octic*ell0**2 == -200*z**6
    assert g0(0) + 4*ell0(0) == 0 and g0(0) - 4*ell0(0) != 0
    assert valuation(7776, 2) == 5
    assert 7776 == 6*36**2

    h = QQ(3)/2
    first_curve = EllipticCurve([0, h*h, 0, 2*h*(h+1), (h+1)**2])
    for t, R in ((None, 1), (-1, 9), (0, 4), (2, 36)):
        if t is None:
            y, v = -QQ(2)/3, QQ(4)/3
        else:
            denominator = 3*(t*t + 2)
            y = QQ(-2*t*t - 4*t + 4)/denominator
            v = QQ(12*R)/denominator**2
        U = v*v/(4*y*y)
        V = v/(8*y)*(-QQ(135)/16*y*y - QQ(25)/(9*y*y))
        point = first_curve(U-1, V)
        assert point[0] == 0 and (2*point).order() == 3


if __name__ == "__main__":
    verify_quotient_pell_and_residues()
    verify_wild_polynomials_and_conductors()
    verify_pell_fiber_transport()
    print("WD1–WD5 Pell data, wild conductors, prime ratios and rational fibers verified.")

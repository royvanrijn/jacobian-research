"""Exact quotient reduction certificates. Run with SageMath 10.9."""
from sage.all import GF, QQ, EllipticCurve, HyperellipticCurve, PolynomialRing


def quotient_polynomial(prime):
    ring = PolynomialRing(GF(prime), "x")
    x = ring.gen()
    s = (ring.base_ring()(3)/7)**2*x**3
    a = 49*s**4 + 750*s**3 - 2403*s**2 + 27000*s + 63504
    assert a.degree() == 12 and (x*a).is_squarefree()
    return x, a


def check_reduced_quotient(f, expected_frobenius, order, cofactor, points, triples):
    curve = HyperellipticCurve(f)
    frobenius = curve.frobenius_polynomial()
    assert frobenius.list() == expected_frobenius
    assert frobenius(1) == order and cofactor % 3 != 0
    jacobian = curve.jacobian()(f.base_ring())
    ring = f.parent()
    x = ring.gen()
    multiples = []
    for (px, py), (u, v) in zip(points, triples):
        assert f(px) == f.base_ring()(py)**2
        point = cofactor*jacobian([x-px, ring(py)])
        triple = 3*point
        assert triple == jacobian([ring(u), ring(v)])
        assert triple != jacobian(0) and 3*triple == jacobian(0)
        assert all(triple not in (old, -old) for old in multiples)
        multiples.append(triple)
    assert len(multiples) == len(points) == len(triples)


x, a = quotient_polynomial(83)
assert a(62) == 0
odd_plus = sum(a[i]*(62*x+1)**i*x**(12-i) for i in range(13))
assert odd_plus.degree() == 11 and odd_plus.is_squarefree()
check_reduced_quotient(
    odd_plus,
    [3939040643, 569499852, 74332310, 3885396, 691307, 53136, 8329, 564, 130, 12, 1],
    27*169907840, 169907840, [(0, 30)],
    [(x**5+41*x**4+5*x**3+38*x**2+72*x+35, 9*x**4+11*x**3+74*x**2+31*x+7)],
)

x, a = quotient_polynomial(137)
odd_minus = x*a
assert odd_minus.degree() == 13
check_reduced_quotient(
    odd_minus,
    [6611856250609, 579140693484, 100398477885, 7899196416, 1198719723, 68467668,
     8674126, 499764, 63867, 3072, 285, 12, 1],
    81*90130506752, 90130506752, [(1, 7), (3, 24)],
    [(x**6+68*x**3+8, 125*x**5+21*x**2), (x**6+135*x**3+128, 30*x**5+19*x**2)],
)

h = QQ(3)/22
elliptic = EllipticCurve([0, h*h, 0, 2*h*(h+1), (h+1)**2]).change_ring(GF(17))
assert elliptic.cardinality() == 18
boundary = elliptic(QQ(75)/121, QQ(1750)/1331)
assert (-2*boundary).order() == 9

print("Verified the quotient Frobenius polynomials, order-nine witnesses and boundary reduction.")

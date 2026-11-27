"""Verify the Tate curve, the two boundary marks, and the fixed local slices."""
from sage.all import EllipticCurve, FunctionField, GF, Matrix, PolynomialRing, QQ, ZZ


v = QQ(13)/2
b = (v*v-1)**3*(v*v+3)/(256*v*v)
c = -1664
assert b == QQ(813074625)/2768896
E = EllipticCurve([1, b, b, 0, 0])
T = E(0, 0)
assert 2*T == E(-b, 0) and T.order() == 4
z = -4/(3+v*v)
xr = b*z
yr = (b*(z+1)*(z+2)/(v*z)-xr-b)/2
R = E(xr, yr)
assert R.order() == 3
slope = (3*xr*xr+2*b*xr-yr)/(2*yr+xr+b)
assert slope == QQ(23625)/1664


def F(t):
    return (b*t*t+t-1)**2-4*b*t**3


ring = PolynomialRing(QQ, "t")
t = ring.gen()
assert F(t).discriminant() == 256*b**4*(1-16*b)
assert F(QQ(1)/17) < 0 < F(-QQ(1)/17)
K = FunctionField(QQ, "t")
t = K.gen()
polynomials = PolynomialRing(K, "X")
xvar = polynomials.gen()
C = K.extension(xvar*xvar + (b*t*t+t-1)*xvar/(t*t) + b/t, "x")
x = C.gen()
t = C(t)
y = t*x*(x+b)
xbar = b/(t*x)
ybar = t*xbar*(xbar+b)
x2 = -b*x/(x+b)
y2 = -b*b*(y+x+b)/(x+b)**2
assert y*y+(x+b)*y == x*x*(x+b)
assert y2*y2+(x2+b)*y2 == x2*x2*(x2+b)
assert y2/(x2*(x2+b)) == 1/(b*t)
assert (y*y/(x+b))*(ybar*ybar/(xbar+b)) == b**3


def marked_values(t, x, b, xr, yr, slope):
    y = t*x*(x+b)
    x2 = -b*x/(x+b)
    y2 = -b*b*(y+x+b)/(x+b)**2
    xbar = b/(t*x)
    ybar = t*xbar*(xbar+b)
    xbar2 = -b*xbar/(xbar+b)
    ybar2 = -b*b*(ybar+xbar+b)/(xbar+b)**2
    line = y2-yr-slope*(x2-xr)
    linebar = ybar2-yr-slope*(xbar2-xr)
    return (y*y/(x+b))**2/b**3, (y2*y2/(x2+b))**2/b**3, line/linebar


assert all(u*v == 1 for u, v in zip(
    marked_values(t, x, C(b), xr, yr, slope),
    marked_values(t, xbar, C(b), xr, yr, slope),
))

ramified = [(2, 30336, 9, 6), (3, -1664, 3, 6), (5, 1, 3, 10),
            (11, 1, 3, 22), (13, 10686, 1, 13), (181, 1, 1, 181)]
for p, center, valuation, local_twist in ramified:
    center = QQ(center)
    radicand = F(center)
    assert radicand.valuation(p) == valuation and valuation % 2
    ratio = center/c
    assert ratio.valuation(p) == 0
    if p == 3:
        assert ratio == 1
    else:
        prime_to_three = (p-1)//3**ZZ(p-1).valuation(3)
        assert GF(p)(ratio)**prime_to_three == 1
    square = radicand/local_twist
    vp = square.valuation(p)
    assert vp % 2 == 0
    unit = square/QQ(p)**vp
    assert unit % 8 == 1 if p == 2 else GF(p)(unit).is_square()
    expected_conductor_exponent = 6 if p == 2 else 2
    assert E.quadratic_twist(local_twist).conductor().valuation(p) == expected_conductor_exponent

assert E.conductor() == 2*3*5*11*13*157*173*181
assert set(b.numerator().prime_factors()+b.denominator().prime_factors()) == {
    p for p, _, _, _ in ramified
}

for p, u, radical, expected in [
    (157, 1, 46, [75, 39, 141]), (173, 2, 37, [149, 10, 7]),
]:
    field = GF(p)
    assert ZZ(p-1).valuation(3) <= 1
    t = field(c*u**3)
    root = field(radical)
    bq = field(b)
    assert root != 0 and root**2 == (bq*t*t+t-1)**2-4*bq*t**3
    x = (1-t-bq*t*t+root)/(2*t*t)
    assert list(marked_values(t, x, bq, field(xr), field(yr), field(slope))) == [
        field(value) for value in expected
    ]
    assert all(expected)

rows = []
for p, u, radical, zeta, expected in [
    (7, 3, 1, 2, [0, 2, 2]), (31, 1, 9, 25, [2, 0, 2]),
    (43, 4, 11, 36, [1, 1, 1]),
]:
    field = GF(p)
    assert ZZ(p-1).valuation(3) == 1
    t = field(c*u**3)
    root = field(radical)
    bq = field(b)
    assert root != 0 and root**2 == (bq*t*t+t-1)**2-4*bq*t**3
    x = (1-t-bq*t*t+root)/(2*t*t)
    vals = marked_values(t, x, bq, field(xr), field(yr), field(slope))
    zeta = field(zeta)
    assert zeta.multiplicative_order() == 3
    assert [a**((p-1)//3) for a in vals] == [zeta**i for i in expected]
    rows.append(expected)
assert Matrix(GF(3), rows).rank() == 3

print("Verified the mixed three-power construction, ramified local centers, and independence rows.")

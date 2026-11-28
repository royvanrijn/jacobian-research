"""Exact algebra for the standalone genus-three construction; run with sage -python.

This checks the maps and marked points, not the written proof of Selmer transfer.
All inputs are the rational constants below; there is no search or random seed.
"""
from sage.all import EllipticCurve, GF, PolynomialRing, QQ

ring = PolynomialRing(QQ, "t")
t = ring.gen()
F = t**8 - 6*t**7 - t**6 + 18*t**5 + 93*t**4 - 36*t**3 - 4*t**2 + 48*t + 16
F7 = F.change_ring(GF(7))
assert F7.gcd(F7.derivative()) == 1
# A repeated-root control must fail the same squarefreeness criterion.
bad = F7 * (F7.parent().gen() - 1)**2
assert bad.gcd(bad.derivative()).degree() > 0

field = ring.fraction_field()
t = field(t)
d0 = 3*(t*t + 2)
x = (t*t - 4*t - 2)/d0
y = (-2*t*t - 4*t + 4)/d0
w_coefficient = 12/d0**2
assert (y - QQ(2)/3)/(x - QQ(1)/3) == t
assert 6*x*x + 3*y*y == 2
assert F*w_coefficient**2 == -QQ(9)/2*x**4 - QQ(117)/16*y**4 + QQ(59)/18

rows = [
    (y, w_coefficient, -QQ(135)/16, QQ(3)/2, QQ(25)/9, QQ(3)/2),
    (x, w_coefficient/2, -QQ(135)/16, QQ(39)/8, QQ(1)/144, -QQ(3)/4),
    (x/y, -w_coefficient/(3*y*y), QQ(25)/9, QQ(59)/18, QQ(1)/144, -QQ(7)/6),
]
for z, v_coefficient, a, b, c, h in rows:
    assert a*c*(b*b - 4*a*c) != 0 and c.is_square()
    v_square = F*v_coefficient**2
    assert v_square == a*z**4 + b*z*z + c
    U = v_square/(4*z*z)
    V_coefficient = v_coefficient*(a*z*z - c/(z*z))/(8*z)
    X = U - 1
    assert F*V_coefficient**2 == X**3 + (h*X + h + 1)**2
    # X on C factors through the degree-two t-map; X on E has degree two.
    # Its reduced rational degree in t is therefore exactly deg(pi).
    assert max(X.numerator().degree(), X.denominator().degree()) == 8
    E = EllipticCurve([0, h*h, 0, 2*h*(h+1), (h+1)**2])
    assert E.discriminant() == 16*(h+1)**3*(4*h**3 - 27*(h+1)) != 0
    P = E(0, -h-1)
    assert not P.is_zero() and (3*P).is_zero()

# Verify each character before imposing the conic relations: signs alone suffice.
ambient = PolynomialRing(QQ, names=("x", "y", "w")).fraction_field()
xx, yy, ww = ambient.gens()
characters = [(1, -1), (-1, 1), (-1, -1)]
for (z, v), (_, _, a, b, c, _), character in zip(
        [(yy, ww), (xx, ww/2), (xx/yy, -ww/(3*yy*yy))], rows, characters):
    U = v*v/(4*z*z)
    V = v*(a*z*z - c/(z*z))/(8*z)
    for images, sign in zip([(-xx, yy, ww), (xx, -yy, ww)], character):
        sign_map = ambient.hom(images, ambient)
        assert sign_map(U) == U and sign_map(V) == sign*V

print("PASS: squarefree octic and repeated-root control; conic inverse and identity;")
print("three smooth quartics; three degree-eight elliptic maps and order-three points;")
print("three distinct nontrivial sign characters. No Selmer or Sha computation.")

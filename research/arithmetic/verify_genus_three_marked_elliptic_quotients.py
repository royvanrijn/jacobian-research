"""Exact quotient, local residue, and conductor checks for HG1–HG4."""
from sage.all import GF, QQ, EllipticCurve, Matrix, NumberField, PolynomialRing

ring = PolynomialRing(QQ, "T")
T = ring.gen()
F = T**8 - 6*T**7 - T**6 + 18*T**5 + 93*T**4 - 36*T**3 - 4*T**2 + 48*T + 16
assert F.is_squarefree() and F.change_ring(GF(7)).is_squarefree()
assert F.is_irreducible()

rational_functions = ring.fraction_field()
t = rational_functions(T)
denominator = 3 * (t*t + 2)
x = (t*t - 4*t - 2) / denominator
y = (-2*t*t - 4*t + 4) / denominator
w_coefficient = 12 / denominator**2
w_square = F * w_coefficient**2
assert 6*x*x + 3*y*y == 2
assert w_square == -QQ(9)/2*x**4 - QQ(117)/16*y**4 + QQ(59)/18

rows = [
    (y, w_coefficient, -QQ(135)/16, QQ(3)/2, QQ(25)/9, QQ(3)/2),
    (x, w_coefficient/2, -QQ(135)/16, QQ(39)/8, QQ(1)/144, -QQ(3)/4),
    (x/y, -w_coefficient/(3*y*y), QQ(25)/9, QQ(59)/18, QQ(1)/144, -QQ(7)/6),
]
points = []
curves = []
for coordinate, v_coefficient, a, b, c, h in rows:
    v_square = F * v_coefficient**2
    assert v_square == a*coordinate**4 + b*coordinate**2 + c
    U = v_square / (4*coordinate**2)
    V_coefficient = v_coefficient * (a*coordinate**2 - c/coordinate**2) / (8*coordinate)
    X = U - 1
    L = 1 + h*U
    assert F*V_coefficient**2 == X**3 + L**2
    assert max(X.numerator().degree(), X.denominator().degree()) == 8
    E = EllipticCurve([0, h*h, 0, 2*h*(h+1), (h+1)**2])
    psi = ring(E.division_polynomial(3).list())
    assert psi % T == 0 and (psi // T).is_irreducible()
    points.append((X, V_coefficient, L))
    curves.append(E)

assert [E.change_ring(GF(11)).trace_of_frobenius() for E in curves] == [6, -6, 0]
assert [E.conductor() for E in curves] == [90, 90, 30]
assert [E.quadratic_twist(3).conductor() for E in curves] == [720, 720, 720]
# HG2 imports only the exponent at 3 from these twists.
assert [E.conductor().valuation(2) for E in curves] == [1, 1, 1]
assert [E.quadratic_twist(3).conductor().valuation(3) for E in curves] == [2, 2, 2]
assert [E.quadratic_twist(3).conductor().valuation(2) for E in curves] == [4, 4, 4]
# At t=1, D=129 is locally split at 2 and 5 and locally 3 at 3.
# This is a conductor fixture, not the full CRT independence slice.
assert 129 % 8 == 1 and GF(5)(129).is_square()
assert (QQ(129)/3) % 3 == 1
assert [E.quadratic_twist(129).conductor() for E in curves] == [10*129**2]*3

assert F(1) == 129 and F(1) % 9 == 3
assert points[0][0] == (-QQ(135)*y*y/16 - QQ(5)/2 + QQ(25)/(9*y*y))/4
assert points[1][0] == (-1215*x**4 + 126*x*x + 1)/(576*x*x)

exponents = []
for q, argument, w, zeta, targets in [
    (13, 1, 6, 3, (9, 10, 4)),
    (19, 4, 7, 7, (15, 3, 17)),
    (31, 1, 6, 5, (20, 10, 8)),
]:
    field = GF(q)
    w = field(w)
    assert w**2 == field(w_square(argument))
    r = w / field(w_coefficient(argument))
    residues = [r * field(Y(argument)) + field(L(argument)) for _, Y, L in points]
    assert tuple(residues) == targets
    x3 = field(zeta)
    assert x3 != 1 and x3**3 == 1
    exponents.append([(field(1), x3, x3**2).index(value**((q-1)//3)) for value in residues])
assert exponents == [[2, 1, 2], [2, 1, 1], [1, 1, 0]]
assert Matrix(GF(3), exponents).rank() == 3

K = NumberField(T*T - 129, "r")
r = K.gen()
assert K.class_group(proof=True).invariants() == ()
unit = 16855 - 1484*r
fundamental = K.unit_group(proof=True).fundamental_units()[0]
assert unit in (fundamental, -fundamental, 1/fundamental, -1/fundamental)
_, V, L = points[2]
alpha3 = V(1)*r + L(1)
assert alpha3 == -QQ(251)/50 - QQ(217)/250*r
assert alpha3 == unit * (-QQ(397)/10 - QQ(7)/2*r)**3

print("Verified the genus-three maps, independent cube residues, conductors, and the unit-class fixture.")

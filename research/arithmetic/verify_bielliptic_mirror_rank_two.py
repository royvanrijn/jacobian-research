"""Exact identity, specialization, and conductor checks for BT1–BT4."""
from sage.all import CRT_list, GF, QQ, EllipticCurve, PolynomialRing, ZZ

ring = PolynomialRing(QQ, "t")
t = ring.gen()
A, B = 61, 19
G = -A * (t**2 + 1) * (B * (t**2 + 1)**2 - A**2 * t**2)
assert G.degree() == 6 and G.is_squarefree()
assert max(f.degree() for f, _ in G.factor()) <= 4

field = ring.fraction_field()
tq = field(t)
X1 = -B * (tq**2 + 1) / A
X2 = X1 / tq**2
F = t * (t**2 + A * t + B)
for index, coordinate in enumerate((X1, X2)):
    assert max(coordinate.numerator().degree(), coordinate.denominator().degree()) == 2
    assert F(coordinate) == B**2 * G / (A**4 * tq**(6 * index))
    x = coordinate - 1
    assert (1 + 8 * coordinate)**2 - F(coordinate) == -x**3

curve = EllipticCurve([0, 64, 0, 144, 81])
assert ring(curve.division_polynomial(3).list()) == t * (3*t**3 + 256*t**2 + 864*t + 972)
assert (3*t**3 + 256*t**2 + 864*t + 972).is_irreducible()
assert curve.conductor() == 1140
assert curve.quadratic_twist(2).conductor().valuation(2) == 6
assert 2**6 * 3 * 5 * 19 == 285 * 8**2

local_y_square = QQ(F(X1(1)))
assert local_y_square / 3**6 == QQ(3610) / 226981
assert (local_y_square / 3**6) % 9 == 1
assert 1 + 8 * X1(1) == 1 + 8 * X2(1) == -QQ(243) / 61
assert G(1) % 16 == 2
assert GF(5)(G(0)) == 1 and GF(19)(G(2)) == 7 and GF(19)(7).is_square()

for q, residue, y1, expected in [(7, 1, 3, (4, 4)), (13, 3, 2, (9, 1))]:
    finite = GF(q)
    X1q, X2q = finite(X1(residue)), finite(X2(residue))
    y1 = finite(y1)
    y2 = y1 / finite(residue)**3
    finite_polynomial = F.change_ring(finite)
    assert y1**2 == finite_polynomial(X1q) and y2**2 == finite_polynomial(X2q)
    assert tuple((y + 1 + 8 * X)**((q - 1) // 3)
                 for X, y in [(X1q, y1), (X2q, y2)]) == expected

moduli = [8, 3**8, 7, 13, 5, 19]
residues = [1, 1, 1, 3, 0, 2]
modulus = ZZ(1)
for q in moduli:
    modulus *= q
t0 = CRT_list(residues, moduli)
b = (1 - 21 * t0) % modulus + modulus
a = (t0 - b) % modulus + modulus
U, V = a + b, 21 * (a + b) + b
specialization = QQ(U) / V
assert 0 < specialization < QQ(1) / 20 and G(specialization) < 0
for q, residue in zip(moduli, residues):
    assert V % q == 1 and U % q == residue

for m in [36, 1225, 225]:
    factors = [m*t**3 - 1, m*t**3 + 1, m**2*t**6 + 4*m*t**3 - 1]
    polynomial = factors[0] * factors[1] * factors[2]
    assert polynomial.degree() == 12 and polynomial.is_squarefree()
    assert max(f.degree() for f in factors) <= 6

print("Verified the bielliptic maps, local cube and independence residues, and counting form hypotheses.")

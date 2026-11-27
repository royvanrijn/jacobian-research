"""Verify boundary coefficients and the higher three-power local tests."""
from pathlib import Path
from runpy import run_path

from sage.all import GF, LaurentSeriesRing, Matrix, NumberField, PolynomialRing, QQ, Qp, ZZ

data = run_path(str(Path(__file__).with_name("verify_marked_octic_mirror_class_rank_four.py")))
F, marks, beta = (data[name] for name in ("F", "points", "beta"))

series_ring = LaurentSeriesRing(QQ, "s", default_prec=50)
s = series_ring.gen()
expected = [
    [(-3, -QQ(7)/11), (-3, QQ(21)/11), (0, QQ(15)/11)],
    [(-3, QQ(7)/11), (-3, -QQ(21)/11), (0, -QQ(5)/11)],
    [(-3, QQ(7)/66), (-3, QQ(7)/22), (0, -QQ(5)/11)],
    [(-3, -QQ(7)/66), (-3, -QQ(7)/22), (0, QQ(15)/11)],
]
for index, (argument, sign) in enumerate([(s, 1), (s, -1), (1/s, 1), (1/s, -1)]):
    value = F(argument)
    v = value.valuation()
    leading = value[v]
    radical = sign*leading.sqrt()*s**(v//2)*(value/(s**v*leading)).sqrt()
    for (coefficient, rational), (order, cube_root) in zip(marks, expected[index]):
        evaluation = coefficient(argument)*radical+rational(argument)
        assert evaluation.valuation() == order
        assert evaluation[order] == cube_root**3

for n in (1, 2, 3):
    m = 3**n
    a = (3*ZZ(4).inverse_mod(m)) % m
    c = QQ(3)/7*11**a
    assert a % 3 == 0 and (4*a-3) % m == 0
    for p in (2, 3, 5, 7, 11):
        direction = 1 if p == 5 else -1
        t = c*QQ(p)**(direction*5*m)
        field = Qp(p, prec=80*m+100)
        radical = field(F(t)).sqrt()
        for sign in (1, -1):
            for coefficient, rational in marks:
                value = field(rational(t))+sign*field(coefficient(t))*radical
                assert value.valuation() % 3 == 0 and value.precision_relative() >= 2
                if p == 3:
                    assert ZZ(value.unit_part().lift()) % 9 in (1, 8)
            b = field(beta[1](t))+sign*field(beta[0](t))*radical
            gamma = b*b/11**3
            assert gamma.valuation() % m == 0 and gamma.precision_relative() > n+1
            if p == 3:
                modulus = 3**(n+1)
                assert ZZ(gamma.unit_part().lift()) % modulus in (1, modulus-1)

rows = []
for p, u, root, zeta, expected_row in [
    (13, 7, 1, 3, [0, 2, 2, 1]), (31, 1, 5, 25, [0, 1, 2, 2]),
    (43, 10, 14, 36, [1, 1, 0, 0]), (61, 4, 17, 47, [0, 2, 1, 0]),
]:
    k = GF(p)
    assert ZZ(p-1).valuation(3) == 1
    t = k(3*u**3)/7
    root, zeta = k(root), k(zeta)
    assert root != 0 and root**2 == F.change_ring(k)(t)
    assert zeta.multiplicative_order() == 3

    def evaluate(rational):
        return rational.numerator().change_ring(k)(t)/rational.denominator().change_ring(k)(t)

    for (coefficient, rational), exponent in zip(marks+[beta], expected_row):
        evaluation = evaluate(rational)+root*evaluate(coefficient)
        assert evaluation**((p-1)//3) == zeta**exponent
    rows.append(expected_row)
assert Matrix(GF(3), rows).rank() == 4

m = 9
c = QQ(3)/7*11**3
constant = (QQ(98)/3)**4/11**6
for p, v in [(73, m), (37, -m)]:
    t = c*QQ(p)**v
    field = Qp(p, prec=160)
    root = field(F(t)).sqrt()
    chosen = [field(beta[1](t))+sign*field(beta[0](t))*root for sign in (1, -1)]
    dominant = min(chosen, key=lambda value: value.valuation())
    gamma = dominant*dominant/11**3
    assert gamma.valuation() == -4*m
    actual = gamma.unit_part().residue()
    expected = GF(p)(constant if v > 0 else 1)
    assert actual**((p-1)//m) == expected**((p-1)//m)
    assert (expected**((p-1)//m)).multiplicative_order() == (m if v > 0 else 1)

# A mixed exponent instance of YR4 with local fifteenth powers at three and five.
m = 15
c = QQ(3)/7*11**12
for p in (3, 5):
    t = c*QQ(p)**(-5*m)
    local = Qp(p, prec=1000)
    root = local(F(t)).sqrt()
    for sign in (1, -1):
        evaluation = local(beta[1](t))+sign*local(beta[0](t))*root
        gamma = evaluation**2/11**3
        assert gamma.valuation() % m == 0 and gamma.precision_relative() >= 3
        assert ZZ(gamma.unit_part().lift())**(p-1) % p**2 == 1

mixed_rows = []
for p, u, radical, zeta, row in [
    (13, 7, 1, 3, [0, 2, 2, 1]), (19, 4, 3, 7, [1, 2, 2, 2]),
    (31, 3, 6, 25, [2, 1, 2, 1]), (37, 2, 10, 26, [2, 2, 2, 2]),
]:
    k = GF(p)
    t = k(c*u**m)
    root = k(radical)
    assert root != 0 and root**2 == F.change_ring(k)(t)
    zeta = k(zeta)
    for (coefficient, rational), exponent in zip(marks+[beta], row):
        value = evaluate(rational)+root*evaluate(coefficient)
        assert value != 0 and value**((p-1)//3) == zeta**exponent
    if p == 31:
        gamma = (evaluate(beta[1])+root*evaluate(beta[0]))**2/k(11**3)
        assert gamma**6 == k(16)
    mixed_rows.append(row)
assert Matrix(GF(3), mixed_rows).rank() == 4

# For a rational local radical, the Dickson algebra has the cyclotomic discriminant.
ring = PolynomialRing(QQ, "Z")
Z = ring.gen()
previous, dickson = ring(2), Z
for d in range(2, 16):
    previous, dickson = dickson, Z*dickson-previous
    if d not in (5, 9, 15):
        continue
    local_model = dickson-(QQ(2)**d+QQ(2)**(-d))
    factors = local_model.factor()
    assert all(multiplicity == 1 for _, multiplicity in factors)
    assert sorted(f.degree() for f, _ in factors) == {
        5: [1, 4], 9: [1, 2, 6], 15: [1, 2, 4, 8]
    }[d]
    discriminant = 1
    for factor, _ in factors:
        if factor.degree() != 1:
            discriminant *= abs(NumberField(factor, "a").discriminant())
    assert discriminant == {5: 5**3, 9: 3**10, 15: 3**5*5**9}[d]

print("Verified the boundary cube coefficients, higher three-adic powers, and prime residue tests.")

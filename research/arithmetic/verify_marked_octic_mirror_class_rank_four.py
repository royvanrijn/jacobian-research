"""Exact chart, norm, local cube and conductor checks for MO1–MO4."""
from sage.all import GF, QQ, EllipticCurve, Matrix, PolynomialRing

ring = PolynomialRing(QQ, "T")
T = ring.gen()
F = 49*T**8 + 750*T**6 - 2403*T**4 + 27000*T**2 + 63504
assert F.is_squarefree() and F.change_ring(GF(13)).is_squarefree()
base_change = F(QQ(3)/7*T**3)
assert base_change.degree() == 24 and base_change.is_squarefree()

function_field = ring.fraction_field()
t = function_field(T)
x = (t + 6/t)/2
y = (6/t - t)/6
w_coefficient = 108 / (121*(6*t)**2)
w_square = F * w_coefficient**2
assert x*x - 9*y*y == 6
assert w_square == (5778*x**4 + 103518*y**4 - 116883)/14641

rows = [
    (y, -QQ(11)/54*w_coefficient, QQ(196)/121, QQ(214)/121, QQ(125)/484, -QQ(16)/11),
    (x, QQ(11)/6*w_coefficient, QQ(196)/121, -QQ(426)/121, -QQ(7875)/484, QQ(24)/11),
    (x/y, 11*w_coefficient/(9*y*y), QQ(125)/484, QQ(1443)/242, -QQ(7875)/484, QQ(3)/22),
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
    curve = EllipticCurve([0, h*h, 0, 2*h*(h+1), (h+1)**2])
    division = ring(curve.division_polynomial(3).list())
    assert division % T == 0 and (division // T).is_irreducible()
    points.append((V_coefficient, L))
    curves.append(curve)

assert [E.change_ring(GF(13)).trace_of_frobenius() for E in curves] == [-4, 2, -4]
assert [E.change_ring(GF(23)).trace_of_frobenius() for E in curves] == [6, 6, 0]
assert [E.conductor() for E in curves] == [1540, 13860, 6930]
for p in [5, 7, 11]:
    assert all(E.quadratic_twist(p).conductor().valuation(p) == 2 for E in curves)

P = 7*t**4 + QQ(375)/7*t*t + 252
beta_coefficient = QQ(7)/18/(t*t)
beta = beta_coefficient, beta_coefficient*P
assert P*P - F == 11*(QQ(198)/7*t*t)**2
assert beta[1]**2 - F*beta[0]**2 == 11**3

radical_at_zero = 252
_, quotient_v, _, _, _, _ = rows[2]
z0 = (x/y)(0)
v0 = quotient_v(0) * radical_at_zero
x0 = v0*v0/(4*z0*z0) - 1
y0 = v0*(QQ(125)/484*z0*z0 + QQ(7875)/(484*z0*z0))/(8*z0)
boundary_point = curves[2](x0, y0)
assert (x0, y0) == (QQ(75)/121, QQ(1750)/1331)
assert 6*boundary_point != curves[2](0)
assert [curves[2].change_ring(GF(p)).cardinality() for p in [13, 23]] == [18, 24]


def assert_split(prime, center):
    value = F(center)
    valuation = value.valuation(prime)
    assert valuation % 2 == 0
    unit = value / QQ(prime)**valuation
    if prime == 2:
        assert unit % 8 == 1
    else:
        assert GF(prime)(unit).is_square()


def split_valuation_data(prime, center, evaluations):
    return [((F(center)*irrational(center)**2).valuation(prime)//2,
             rational(center).valuation(prime))
            for irrational, rational in evaluations]


def assert_cube_residues(center, evaluations, scales, expected):
    assert_split(3, center)
    for (irrational, rational), exponent, residues in zip(evaluations, scales, expected):
        assert exponent % 3 == 0
        squared_term = F(center) * (3**exponent * irrational(center))**2
        rational_term = 3**exponent * rational(center)
        assert (squared_term % 9, rational_term % 9) == residues
        assert any((k*k) % 9 == residues[0] and (k+residues[1]) % 9 in (1, 8)
                   for k in range(9))


for p, center, expected in [
    (2, -QQ(1)/2, [(-6, 1), (-6, -1), (1, 0)]),
    (11, -QQ(1), [(-3, -3), (-3, -3), (-3, -3)]),
]:
    assert_split(p, center)
    assert split_valuation_data(p, center, points) == expected
for p, center in [(5, -QQ(1)/2), (7, -QQ(1))]:
    assert F(center).valuation(p) % 2 == 1
assert_cube_residues(-QQ(2)/3, points, [6, 3, 0], [(1, 0), (1, 0), (7, 4)])

four_functions = points + [beta]
for p, center, expected in [
    (2, QQ(3)/56, [(-12, -4), (-12, -5), (1, 0), (-7, -7)]),
    (7, QQ(3)/2401, [(-9, -6), (-9, -6), (1, 0), (-6, -6)]),
]:
    assert_split(p, center)
    assert split_valuation_data(p, center, four_functions) == expected
for p, center in [(5, -QQ(1)/2), (11, QQ(19))]:
    assert F(center).valuation(p) % 2 == 1
assert_cube_residues(QQ(1)/63, four_functions, [9, 6, 0, 6],
                     [(1, 0), (1, 0), (7, 4), (7, 5)])


def character_row(prime, argument, w, root, expected, evaluations):
    field = GF(prime)
    w = field(w)
    assert w*w == field(w_square(argument))
    radical = w / field(w_coefficient(argument))
    residues = [radical*field(irrational(argument)) + field(rational(argument))
                for irrational, rational in evaluations]
    assert tuple(residues) == expected
    zeta = field(root)
    assert zeta != 1 and zeta**3 == 1
    return [(field(1), zeta, zeta**2).index(value**((prime-1)//3)) for value in residues]


direct = [
    (13, 10, 2, 3, (10, 4, 7)),
    (19, 3, 1, 7, (13, 15, 2)),
    (31, 2, 14, 5, (10, 12, 5)),
]
direct_rows = [character_row(*row, points) for row in direct]
assert direct_rows == [[1, 2, 2], [2, 2, 1], [1, 2, 1]]
assert Matrix(GF(3), direct_rows).rank() == 3

changed = [
    (19, 4, 3, 1, 7, (13, 15, 2, 4)),
    (31, 3, 16, 14, 5, (14, 22, 11, 21)),
    (37, 3, 1, 3, 10, (9, 4, 7, 35)),
    (43, 21, 13, 13, 6, (3, 26, 39, 35)),
]
changed_rows = []
for q, u, argument, w, zeta, expected in changed:
    assert GF(q)(3*u**3) / GF(q)(7) == argument
    changed_rows.append(character_row(q, argument, w, zeta, expected, four_functions))
assert changed_rows == [[2, 2, 1, 2], [2, 1, 1, 1], [2, 1, 1, 2], [2, 1, 0, 0]]
assert Matrix(GF(3), changed_rows).rank() == 4

print("Verified the octic maps, fourth norm, local cube conditions, independence, and conductor inputs.")

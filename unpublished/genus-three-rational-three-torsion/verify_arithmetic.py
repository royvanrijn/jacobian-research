"""Exact fixed-curve arithmetic and marked-norm tangent certificate (Sage 10.9).

No height search or random input is used. The final matrix calculation concerns
only the displayed marked-norm scheme, not rational or geometrically simple
deformations. The written arguments connect these finite certificates to HG6.
"""
from sage.all import (EllipticCurve, FunctionField, GF, PolynomialRing, QQ,
                      gcd, matrix, prod, vector)

ring = PolynomialRing(QQ, "t")
t = ring.gen()
F = t**8-6*t**7-t**6+18*t**5+93*t**4-36*t**3-4*t**2+48*t+16
NORM_ROWS = [
    (t*t+5*t+7,
     t**6+2*t**5-13*t**4-52*t**3-49*t*t+38*t+22, (t+1)**2, QQ(-300)),
    (t*t+2,
     t**6-3*t**5-3*t**4-12*t**3+6*t*t-12*t-8, t*t, QQ(-200)),
    (t*t-4*t-2,
     t**6-7*t**5+5*t**4+20*t**3+50*t*t-52*t+40, (t-2)**2, QQ(24)),
]
for ell, g, q, c in NORM_ROWS:
    assert g*g-F*ell*ell == c*q**3
    assert gcd(g, ell*q) == 1

curves = []
for h in [QQ(3)/2, -QQ(3)/4, -QQ(7)/6]:
    E = EllipticCurve([0, h*h, 0, 2*h*(h+1), (h+1)**2])
    # Disable the rank database: obtain proved rank bounds from descent.
    assert E.rank(proof=True, use_database=False) == 0
    assert E.rank_bounds() == (0, 0)
    curves.append(E)
assert [E.torsion_order() for E in curves] == [6, 6, 12]
assert [E.conductor() for E in curves] == [90, 90, 30]
assert all(E.j_invariant() != 0 for E in curves)
assert prod(E.conductor() for E in curves) == 243000
assert {tuple(P) for P in curves[0].torsion_points()} == {
    (0, 1, 0), (-1, 0, 1), (0, QQ(-5)/2, 1), (0, QQ(5)/2, 1),
    (5, -15, 1), (5, 15, 1)}
counts = []
for p, expected in [(7, 432), (11, 1296), (13, 3888), (17, 3456), (19, 13824)]:
    assert F.change_ring(GF(p)).is_squarefree()
    count = prod(E.change_ring(GF(p)).cardinality() for E in curves)
    assert count == expected
    counts.append(count)
assert gcd(counts) == 432

# Exhaust the first elliptic target, including the pole fibre z=0.
for U, expected in [(0, set()), (1, {QQ(-2)/3, QQ(2)/3}), (6, set())]:
    residual = -QQ(135)/16*t**4+(QQ(3)/2-4*U)*t*t+QQ(25)/9
    assert {a for a, _ in residual.roots(QQ)} == expected
assert not (t*t-QQ(1)/3).roots(QQ)  # z=0 has no rational conic preimage.
for a, b in [(-1, 9), (0, 4), (2, 36)]:
    assert F(a) == b*b
d = 3*(t*t+2)
xx = (t*t-4*t-2)/d
yy = (-2*t*t-4*t+4)/d
fibres = set()
for y0 in [QQ(-2)/3, QQ(2)/3]:
    fibres.update(a for a, _ in (yy-y0).numerator().roots(QQ))
assert fibres == {-1, 0, 2}
assert yy.numerator().leading_coefficient()/yy.denominator().leading_coefficient() == -QQ(2)/3


def function_field(k):
    K = FunctionField(k, "t")
    tt = K.gen()
    S = PolynomialRing(K, "R")
    L = K.extension(S.gen()**2-K(F(tt)), "R")
    rr, tt = L.gen(), L(tt)
    infs = L.places_infinite()
    plus = next(P for P in infs if (rr/tt**4-1).valuation(P) > 0)
    minus = next(P for P in infs if P != plus)
    return L, tt, rr, plus, minus


L, tt, rr, plus, minus = function_field(QQ)
P0 = L.maximal_order().ideal(tt, rr-4).place()
# A polynomial Pell certificate: div(B+A R)=12(infinity_- - infinity_+).
A = -t**8+8*t**7-4*t**6-52*t**5-115*t**4+296*t**3+746*t*t+764*t-139
B = (-t**12+11*t**11-23*t**10-74*t**9-3*t**8+1053*t**7+621*t**6
     -2916*t**5-8697*t**4-3923*t**3+487*t*t+278*t-4444)
assert (L(A(tt))*rr+L(B(tt))).divisor() == 12*(minus.divisor()-plus.divisor())
ell, g, _, _ = NORM_ROWS[1]
assert (L(ell(tt))*rr+L(g(tt))).divisor() == 6*(P0.divisor()-plus.divisor())

# Certify the *whole* Sylow 2 group at 7 by producing 16 distinct elements.
L7, u, v, plus7, minus7 = function_field(GF(7))
G7 = L7.jacobian(model="hess", base_div=3*plus7.divisor()).group()
zero = G7.zero()
I7 = G7.point(minus7.divisor()-plus7.divisor())
Q7 = L7.maximal_order().ideal(u*u+1, v-u+2).place()
P07 = L7.maximal_order().ideal(u, v-4).place()
generators = [27*I7, 27*G7.point(P07.divisor()-plus7.divisor()),
              27*G7.point(Q7.divisor()-2*plus7.divisor())]
two_group = [zero]
for P in generators:
    assert 4*P == zero
    old = list(two_group)
    Q = P
    while Q not in old:
        two_group.extend(Q+R for R in old)
        Q += P
assert len(two_group) == 16
assert sum(2*P == zero for P in two_group) == 4
assert 4*(3*I7) == zero and 2*(3*I7) != zero

# The three norm classes 2[(a,b)-infinity_+] are independent modulo 7.
three_group = [zero]
for a, b in [(-1, 9), (0, 4), (2, 36)]:
    P = L7.maximal_order().ideal(u-a, v-b).place()
    D = 2*G7.point(P.divisor()-plus7.divisor())
    assert D != zero and 3*D == zero
    assert D not in three_group
    old = list(three_group)
    three_group = old+[D+R for R in old]+[2*D+R for R in old]
assert len(three_group) == 27

# A nonzero x-T image obstructs halving T=3[(0,4)-infinity_+] at 17.
R17 = PolynomialRing(GF(17), "z")
z = R17.gen()
q17 = z*z+5
assert q17.is_irreducible() and R17(F.list()) % q17 == 0
k17 = GF(17**2, "a17", modulus=q17)
theta = k17.gen()
assert (-theta)**18 == k17(5)  # Norm from F_(17^2) to F_17.
assert not GF(17)(5).is_square()
assert all(not (k17(c)*(-theta)).is_square() for c in GF(17) if c)
# Square positive control for the same residue-field test.
assert (theta*theta).is_square()
print("PASS: rank zero; counts and gcd; exactly eight rational points;")
print("Sylow 2 at 7=(Z/4)^2; global Sylow 2=Z/2+Z/4 by the written 17-obstruction;")
print("J(Q)=Z/3+Z/6+Z/12 (order 216); conductor 243000; three independent norm rows.")

# First derivatives of g_i^2-f*l_i^2-c_i*q_i^3, with monic degrees 8,2,6,2.
# Variables: 8 f coefficients, then (2 l,6 g,2 q,1 c) for each of 3 rows.
def column(blocks):
    return vector(QQ, [poly[j] for poly in blocks for j in range(12)])


columns = [column([-ell*ell*t**j for ell, g, q, c in NORM_ROWS]) for j in range(8)]
for i, (ell, g, q, c) in enumerate(NORM_ROWS):
    derivatives = ([-2*F*ell, -2*F*ell*t]+[2*g*t**j for j in range(6)]
                   +[-3*c*q*q, -3*c*q*q*t, -q**3])
    for derivative in derivatives:
        blocks = [ring(0)]*3
        blocks[i] = derivative
        columns.append(column(blocks))
M = matrix(columns).transpose()
assert M.dimensions() == (36, 41) and M.rank() == 36
kernel = M.right_kernel().basis_matrix()
f_basis = kernel[:, :8].row_space().basis_matrix()
assert kernel.nrows() == f_basis.nrows() == 5

# Permit the base involution itself to move: [a,b;c,-a], and its multiplier.
for a, b, c in [(2, 2, 1), (-1, 2, 1), (0, -2, 1)]:
    X, Z = a*t+b, c*t-a
    transformed = sum(F[j]*X**j*Z**(8-j) for j in range(9))
    multiplier = transformed//F
    assert transformed == multiplier*F and multiplier.degree() == 0
    involution_columns = []
    for row in f_basis:
        df = sum(row[j]*t**j for j in range(8))
        involution_columns.append(sum(row[j]*X**j*Z**(8-j) for j in range(8))-multiplier*df)
    da = db = dc = ring(0)
    for j in range(9):
        if j:
            da += F[j]*j*t*X**(j-1)*Z**(8-j)
            db += F[j]*j*X**(j-1)*Z**(8-j)
        if j < 8:
            da -= F[j]*(8-j)*X**j*Z**(7-j)
            dc += F[j]*(8-j)*t*X**j*Z**(7-j)
    involution_columns += [da, db, dc, -F]
    B = matrix(QQ, [[poly[j] for poly in involution_columns] for j in range(9)])
    assert B.right_kernel().basis_matrix()[:, :5].rank() == 3
print("PASS: marked-norm Jacobian rank 36/41; tangent dimension 5;")
print("each known involution persists in only a 3-dimensional tangent subspace.")

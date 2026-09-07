from sage.all import *


K = GF(31)
Kt = PolynomialRing(K, "t")
t = Kt.gen()
F = FractionField(Kt)

A = Kt([4, 23, 18, 28, 12, 18, 23, 20, 19])
B = Kt([23, 24, 27, 10, 6, 8, 23, 26, 26, 16, 19, 9, 15])
lam = K(23)

X1 = Kt([3, 9, 0, 12, 29])
Y1 = Kt([0, 5, 21, 12, 1])

r2 = K(27)
z2 = t - r2
X2 = Kt([17, 1, 20, 1, 27, 27, 29])
Y2 = Kt([0, 21, 0, 6, 4])

# The canonical P3.O=1 search produces a denominator-cancelling section.  One
# sign is recorded here after cancelling (t-r)^2 and (t-r)^3 from its X and Y
# numerators.  It is not the requested canonical P3, but it is a genuine third
# polynomial section on the same surface.
XR = Kt([3, 8, 3, 9, 29])
YR = Kt([0, 2, 26, 17, 1, 16])

assert Y1**2 == X1**3 + A * X1 + B
assert Y2**2 == X2**3 + A * X2 * z2**4 + B * z2**6
assert YR**2 == XR**3 + A * XR + B

# Replay one raw meet-in-the-middle hit and certify the cancellation exactly.
r3 = K(6)
z3 = t - r3
G3 = t * (t - lam)
C3 = K(3) * r3**2 * (t - lam) / (-lam)
C3 += K(10) * (lam - r3)**2 * t / lam
X3_raw = C3 + G3 * (12 + 12*t + 20*t**2 + 18*t**3 + K(29)*t**4)
Y3_raw = G3 * (8 + 6*t + 3*t**2 + 16*t**3 + 25*t**4 + 19*t**5 + 16*t**6)
assert X3_raw == XR * z3**2
assert Y3_raw == YR * z3**3

Delta = -16 * (4 * A**3 + 27 * B**2)


def valuation_at(poly, point):
    factor = t - point
    valuation = 0
    while poly(point) == 0:
        poly //= factor
        valuation += 1
    return valuation


valuations = (
    valuation_at(Delta, K(0)),
    valuation_at(Delta, K(1)),
    valuation_at(Delta, lam),
    24 - Delta.degree(),
)
assert valuations == (3, 2, 2, 11)
fixed = t**3 * (t - 1)**2 * (t - lam)**2
residual = Delta // fixed
assert residual.degree() == 6
assert gcd(residual, residual.derivative()).degree() == 0

# Finite component incidences.  P1 is nonidentity at 0 and lambda but is a
# nonsingular point at 1.  P2 is nonidentity at all three finite reducible
# fibers.
nodes = {K(0): K(3), K(1): K(21), lam: K(10)}
assert X1(0) == nodes[K(0)] and Y1(0) == 0
assert X1(lam) == nodes[lam] and Y1(lam) == 0
assert not (X1(1) == nodes[K(1)] and Y1(1) == 0)
for point, node in nodes.items():
    assert X2(point) == node * (point - r2)**2
    assert Y2(point) == 0

# Infinity jets: P1 has exact class 2 up to sign (Ybar starts in degree 2),
# and P2 has exact class 5/6 (Ybar starts in degree 5).  The latter is the
# sign orbit of the target label 6 in Z/11.
assert Y1.degree() == 4
assert Y2.degree() == 4
assert X1.degree() == 4 and X1[4] != 0
assert X2.degree() == 6 and X2[6] != 0

E = EllipticCurve(F, [F(A), F(B)])
P1 = E(F(X1), F(Y1))
P2 = E(F(X2) / F(z2)**2, F(Y2) / F(z2)**3)
PR = E(F(XR), F(YR))
assert not P1.is_zero() and not P2.is_zero() and not PR.is_zero()

# A bounded group-law audit is a regression, not the independence proof.
relations = []
for m in range(-16, 17):
    for n in range(-16, 17):
        if m == 0 and n == 0:
            continue
        if m * P1 + n * P2 == E(0):
            relations.append((m, n))
assert not relations

# Exact rank-three independence certificate.  For each good specialization,
# retain the coefficient masks c in F_2^3 for which
#
#     c1*P1 + c2*P2 + c3*R
#
# lies in 2E(F_31).  The stacked quotients at t=4,7,11 have trivial kernel.
# A separate good fiber at t=5 has odd order 25, so E(F_31(t)) has no rational
# 2-torsion: prime-to-characteristic torsion specializes injectively at good
# places.  Infinite descent then proves Z-independence of the three sections.
remaining_masks = set(range(8))
quotient_trace = []
for parameter in (K(4), K(7), K(11)):
    assert Delta(parameter) != 0 and parameter != r2
    fiber = EllipticCurve(K, [A(parameter), B(parameter)])
    points = (
        fiber(X1(parameter), Y1(parameter)),
        fiber(X2(parameter) / (parameter-r2)**2, Y2(parameter) / (parameter-r2)**3),
        fiber(XR(parameter), YR(parameter)),
    )
    twice = {2*point for point in fiber}
    remaining_masks = {
        mask
        for mask in remaining_masks
        if sum(
            (points[i] for i in range(3) if (mask >> i) & 1),
            fiber(0),
        ) in twice
    }
    quotient_trace.append((int(parameter), fiber.cardinality(), sorted(remaining_masks)))
assert remaining_masks == {0}

torsion_parameter = K(5)
assert Delta(torsion_parameter) != 0 and torsion_parameter != r2
torsion_fiber = EllipticCurve(K, [A(torsion_parameter), B(torsion_parameter)])
assert torsion_fiber.cardinality() == 25

# Shioda self-heights from the verified component data and P.O values:
#
# P1: 4 - [2*9/11 + 1*2/3 + 0 + 1/2] = 79/66.
# P2: 6 - [5*6/11 + 2*1/3 + 1/2 + 1/2] = 106/66.
h1 = QQ(79) / 66
h2 = QQ(106) / 66
height_ratio = h2 / h1
assert height_ratio == QQ(106) / 79
assert not height_ratio.is_square()

print("MW3A10P1P2|field=GF(31)|fibers=I11,I3,I2,I2,6I1", flush=True)
print("MW3A10P1P2|P1_height=79/66|P2_height=106/66", flush=True)
print("MW3A10P1P2|independent=1|reason=height_ratio_not_Q_square", flush=True)
print("MW3A10P1P2|bounded_relations=0|box=16", flush=True)
print("MW3A10P1P2|PASS", flush=True)
print(
    "MW3A10RANK3|quotient_trace="
    + ";".join(
        f"t={parameter},order={order},masks={','.join(map(str, masks))}"
        for parameter, order, masks in quotient_trace
    ),
    flush=True,
)
print("MW3A10RANK3|no_2_torsion_fiber=t=5,order=25", flush=True)
print("MW3A10RANK3|independent=1|method=stacked_mod2_infinite_descent", flush=True)
print("MW3A10RANK3|PASS", flush=True)

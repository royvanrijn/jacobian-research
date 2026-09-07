from sage.all import *


K = GF(31)
Kt = PolynomialRing(K, "t")
t = Kt.gen()
F = FractionField(Kt)

A = Kt([20, 23, 11, 10, 8, 8])
B = Kt([8, 20, 17, 0, 19, 27, 6, 28, 1])

X1 = Kt([11, 27, 29])
Y1 = Kt([3, 2, 19, 6, 1])

r2 = K(8)
z2 = t-r2
X2 = Kt([5, 21, 15, 27, 22])
Y2 = Kt([0, 0, 12, 11, 22, 0, 16, 1])

assert Y1**2 == X1**3 + A*X1 + B
assert Y2**2 == X2**3 + A*X2*z2**4 + B*z2**6
assert X2(r2) != 0

Delta = -16*(4*A**3 + 27*B**2)


def valuation_at(poly, point):
    factor = t-point
    valuation = 0
    while poly != 0 and poly(point) == 0:
        poly //= factor
        valuation += 1
    return valuation


lam, mu = K(10), K(23)
valuations = (
    valuation_at(Delta, K(0)),
    valuation_at(Delta, K(1)),
    valuation_at(Delta, lam),
    valuation_at(Delta, mu),
    24-Delta.degree(),
)
assert valuations == (4, 4, 2, 2, 8)
fixed = t**4*(t-1)**4*(t-lam)**2*(t-mu)**2
residual = Delta // fixed
assert residual.degree() == 4
assert gcd(residual, residual.derivative()).degree() == 0
assert all(residual(point) != 0 for point in (K(0), K(1), lam, mu))

nodes = {K(0): K(18), K(1): K(5), lam: K(16), mu: K(14)}
assert not (X1(0) == nodes[K(0)] and Y1(0) == 0)
assert X1(1) == nodes[K(1)] and Y1(1) == 0
assert not (X1(lam) == nodes[lam] and Y1(lam) == 0)
assert not (X1(mu) == nodes[mu] and Y1(mu) == 0)
for point, node in nodes.items():
    assert X2(point) == node*(point-r2)**2
    assert Y2(point) == 0

E = EllipticCurve(F, [F(A), F(B)])
P1 = E(F(X1), F(Y1))
P2 = E(F(X2)/F(z2)**2, F(Y2)/F(z2)**3)
assert not P1.is_zero() and not P2.is_zero()


def meets_node(section, point, node):
    if section.is_zero():
        return False
    x_coordinate = F(section[0])
    y_coordinate = F(section[1])
    if (
        Kt(x_coordinate.denominator())(point) == 0
        or Kt(y_coordinate.denominator())(point) == 0
    ):
        return False
    return x_coordinate(point) == node and y_coordinate(point) == 0


# The numerator square is genuine, but the deeper I4 test rejects it as the
# canonical P2.  At both I4 fibers, 2*P2 is on the identity component, so P2
# has class 2 at both; the target finite profile is (1,2,1,1).
assert not meets_node(2*P2, K(0), nodes[K(0)])
assert not meets_node(2*P2, K(1), nodes[K(1)])

# Regression only; exact independence is certified below.
relations = []
for m in range(-16, 17):
    for n in range(-16, 17):
        if (m, n) != (0, 0) and m*P1+n*P2 == E(0):
            relations.append((m, n))
assert not relations

# Stack E_t(F_31)/2E_t(F_31) at good fibers until the two coefficient masks
# have trivial common kernel.  A good odd-order fiber rules out rational
# 2-torsion.  Infinite descent then proves P1,P2 are Z-independent.
remaining_masks = set(range(4))
quotient_trace = []
odd_fiber = None
bad = {K(0), K(1), lam, mu, r2}
for parameter in K:
    if parameter in bad or Delta(parameter) == 0:
        continue
    fiber = EllipticCurve(K, [A(parameter), B(parameter)])
    points = (
        fiber(X1(parameter), Y1(parameter)),
        fiber(
            X2(parameter)/(parameter-r2)**2,
            Y2(parameter)/(parameter-r2)**3,
        ),
    )
    twice = {2*point for point in fiber}
    remaining_masks = {
        mask for mask in remaining_masks
        if sum(
            (points[index] for index in range(2) if (mask >> index) & 1),
            fiber(0),
        ) in twice
    }
    quotient_trace.append(
        (int(parameter), fiber.cardinality(), sorted(remaining_masks))
    )
    if odd_fiber is None and fiber.cardinality() % 2:
        odd_fiber = (int(parameter), fiber.cardinality())
    if remaining_masks == {0} and odd_fiber is not None:
        break

assert remaining_masks == {0}
assert odd_fiber is not None

print("E6MW3WRONGP2|field=GF(31)|fibers=IV*,I4,I4,I2,I2,4I1", flush=True)
print("E6MW3WRONGP2|delta_valuations=8,4,4,2,2|residual_squarefree=1", flush=True)
print("E6MW3WRONGP2|finite_profile=2,2,1,1|target=1,2,1,1", flush=True)
print("E6MW3WRONGP2|bounded_relations=0|box=16", flush=True)
print(
    "E6MW3WRONGP2|quotient_trace="
    + ";".join(
        f"t={parameter},order={order},masks={','.join(map(str, masks))}"
        for parameter, order, masks in quotient_trace
    ),
    flush=True,
)
print(
    f"E6MW3WRONGP2|no_2_torsion_fiber=t={odd_fiber[0]},order={odd_fiber[1]}",
    flush=True,
)
print("E6MW3WRONGP2|independent=1|method=stacked_mod2_infinite_descent", flush=True)
print("E6MW3WRONGP2|status=REJECTED_COMPONENT_PROFILE", flush=True)

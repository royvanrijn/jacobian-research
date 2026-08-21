from sage.all import *


K = GF(11)
T = PolynomialRing(K, "t")
t = T.gen()
F = FractionField(T)

A = T([6, 9, 2, 0, 0, 0])
B = T([10, 6, 10, 8, 7, 6, 10, 7, 1])

X1 = T([10, 3, 1])
Y1 = T([5, 9, 9, 9, 1])

r2 = K(10)
z2 = t-r2
X2 = T([3, 6, 3, 9, 2])
Y2 = T([0, 7, 9, 3, 8, 4, 1, 1])

X3 = T([3, 6, 1])
Y3 = T([0, 3, 2, 2, 10])

assert Y1**2 == X1**3+A*X1+B
assert Y2**2 == X2**3+A*X2*z2**4+B*z2**6
assert Y3**2 == X3**3+A*X3+B
assert X2(r2) != 0

Delta = -16*(4*A**3+27*B**2)


def valuation_at(poly, point):
    factor = t-point
    valuation = 0
    while poly != 0 and poly(point) == 0:
        poly //= factor
        valuation += 1
    return valuation


lam, mu = K(3), K(9)
valuations = (
    valuation_at(Delta, K(0)),
    valuation_at(Delta, K(1)),
    valuation_at(Delta, lam),
    valuation_at(Delta, mu),
    24-Delta.degree(),
)
assert valuations == (4, 4, 2, 2, 8)
fixed = t**4*(t-1)**4*(t-lam)**2*(t-mu)**2
residual = Delta//fixed
assert residual.degree() == 4
assert gcd(residual, residual.derivative()).degree() == 0
assert all(residual(point) != 0 for point in (K(0), K(1), lam, mu))


def node_at(point):
    nodes = [
        node for node in K
        if A(point) == -3*node**2
        and B(point) == 2*node**3
        and B.derivative(t)(point)+node*A.derivative(t)(point) == 0
    ]
    assert len(nodes) == 1
    return nodes[0]


points = (K(0), K(1), lam, mu)
nodes = tuple(node_at(point) for point in points)
E = EllipticCurve(F, [F(A), F(B)])
P1 = E(F(X1), F(Y1))
P2 = E(F(X2)/F(z2)**2, F(Y2)/F(z2)**3)
P3 = E(F(X3), F(Y3))


def meets_node(section, point, node):
    if section.is_zero():
        return False
    x_coordinate = F(section[0])
    y_coordinate = F(section[1])
    if (
        T(x_coordinate.denominator())(point) == 0
        or T(y_coordinate.denominator())(point) == 0
    ):
        return False
    return x_coordinate(point) == node and y_coordinate(point) == 0


def meets_ivstar_singular(section):
    if section.is_zero():
        return False
    x_coordinate = F(section[0])
    y_coordinate = F(section[1])
    return (
        T(x_coordinate.numerator()).degree()
        - T(x_coordinate.denominator()).degree() < 4
        and T(y_coordinate.numerator()).degree()
        - T(y_coordinate.denominator()).degree() < 6
    )


# These exact group-law tests orient the component classes.  In particular,
# P2 and P3 have the same I4(0) orientation, whereas the target basis requires
# opposite orientations.
assert tuple(meets_node(P1, point, node) for point, node in zip(points, nodes)) == (False, True, False, False)
assert tuple(meets_node(P2, point, node) for point, node in zip(points, nodes)) == (True, True, True, True)
assert tuple(meets_node(P3, point, node) for point, node in zip(points, nodes)) == (True, False, False, False)
assert tuple(meets_node(2*P2, point, node) for point, node in zip(points, nodes)) == (True, False, False, False)
assert tuple(meets_node(2*P3, point, node) for point, node in zip(points, nodes)) == (True, False, False, False)
assert meets_node(P2+P3, points[0], nodes[0])
assert all(meets_ivstar_singular(point) for point in (P1, P2, P3))
assert not meets_ivstar_singular(P1-P2)
assert not meets_ivstar_singular(P1+P3)
assert not meets_ivstar_singular(P2+P3)


def section_O_intersection(section):
    x_coordinate = F(section[0])
    numerator = T(x_coordinate.numerator())
    denominator = T(x_coordinate.denominator())
    finite_order = denominator.degree()
    infinity_order = max(0, numerator.degree()-denominator.degree()-4)
    assert finite_order % 2 == 0 and infinity_order % 2 == 0
    return ZZ((finite_order+infinity_order)//2)


def component_add(left, right):
    return tuple(
        (a+b) % modulus
        for a, b, modulus in zip(left, right, (3, 4, 4, 2, 2))
    )


def local_self(label):
    e6, i4a, i4b, i2a, i2b = label
    result = QQ(0) if e6 == 0 else QQ(4)/3
    for modulus, value in ((4, i4a), (4, i4b), (2, i2a), (2, i2b)):
        if value:
            result += QQ(value*(modulus-value))/modulus
    return result


def height(section, label):
    if section.is_zero():
        return QQ(0)
    return QQ(4)+2*section_O_intersection(section)-local_self(label)


labels = (
    (1, 0, 1, 0, 0),
    (1, 1, 2, 1, 1),
    (2, 1, 0, 0, 0),
)
sections = (P1, P2, P3)
diagonal = tuple(height(section, label) for section, label in zip(sections, labels))
H = matrix(QQ, 3, 3)
for i in range(3):
    H[i, i] = diagonal[i]
    for j in range(i):
        sum_height = height(sections[i]+sections[j], component_add(labels[i], labels[j]))
        H[i, j] = H[j, i] = (sum_height-diagonal[i]-diagonal[j])/2

target = matrix(QQ, [[23, -10, -8], [-10, 23, 1], [-8, 1, 23]])/12
assert H == matrix(QQ, [[23, -22, 16], [-22, 23, -17], [16, -17, 23]])/12
assert H.det() == QQ(13)/48
assert target.det() == QQ(79)/16
assert H != target

# Stacked specialization in E_t(F_11)/2E_t(F_11), followed by infinite
# descent using one odd-order good fiber, proves Z-independence.
remaining_masks = set(range(8))
quotient_trace = []
odd_fiber = None
bad = set(points) | {r2}
for parameter in K:
    if parameter in bad or Delta(parameter) == 0:
        continue
    fiber = EllipticCurve(K, [A(parameter), B(parameter)])
    specialized = (
        fiber(X1(parameter), Y1(parameter)),
        fiber(X2(parameter)/(parameter-r2)**2, Y2(parameter)/(parameter-r2)**3),
        fiber(X3(parameter), Y3(parameter)),
    )
    twice = {2*point for point in fiber}
    remaining_masks = {
        mask for mask in remaining_masks
        if sum(
            (specialized[index] for index in range(3) if (mask >> index) & 1),
            fiber(0),
        ) in twice
    }
    quotient_trace.append((int(parameter), fiber.cardinality(), sorted(remaining_masks)))
    if odd_fiber is None and fiber.cardinality() % 2:
        odd_fiber = (int(parameter), fiber.cardinality())
    if remaining_masks == {0} and odd_fiber is not None:
        break

assert remaining_masks == {0}
assert odd_fiber is not None

print("E6WRONGRANK3|field=GF(11)|fibers=IV*,I4,I4,I2,I2,4I1", flush=True)
print("E6WRONGRANK3|profiles=P1:1,0,1,0,0;P2:1,1,2,1,1;P3:2,1,0,0,0", flush=True)
print(
    "E6WRONGRANK3|O_intersections="
    + ",".join(
        f"{name}:{section_O_intersection(point)}"
        for name, point in (
            ("P1", P1), ("P2", P2), ("P3", P3),
            ("P1+P2", P1+P2), ("P1+P3", P1+P3),
            ("P2+P3", P2+P3), ("P1-P2", P1-P2),
        )
    ),
    flush=True,
)
print("E6WRONGRANK3|height_gram=(1/12)*[23,-22,16;-22,23,-17;16,-17,23]", flush=True)
print("E6WRONGRANK3|height_det=13/48|target_det=79/16", flush=True)
print(
    "E6WRONGRANK3|quotient_trace="
    + ";".join(
        f"t={parameter},order={order},masks={','.join(map(str, masks))}"
        for parameter, order, masks in quotient_trace
    ),
    flush=True,
)
print(f"E6WRONGRANK3|no_2_torsion_fiber=t={odd_fiber[0]},order={odd_fiber[1]}", flush=True)
print("E6WRONGRANK3|independent=1|status=REJECTED_HEIGHT_LATTICE", flush=True)

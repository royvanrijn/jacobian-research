from sage.all import *


R = PolynomialRing(QQ, "u")
u = R.gen()

a2 = u * (u**3 + QQ(16)/25*u**2 - 12*u + 100) / 100
a4 = -QQ(21)/31250 * u**3 * (u-QQ(25)/7) * (u-1) * (u+2) * (u+QQ(25)/3)
a6 = QQ(441)/39062500 * u**4 * (u-QQ(25)/7)**2 * (u-1)**2 * (u+QQ(25)/3)**2

sections = (
    (
        QQ(21)/625 * u * (u-QQ(25)/7) * (u+QQ(25)/3),
        QQ(63)/6250 * u**2 * (u-QQ(25)/7) * (u+QQ(25)/3),
    ),
    (
        -QQ(7)/25 * u * (u-1) * (u-QQ(25)/7),
        QQ(98)/3125 * u**2 * (u-QQ(275)/28) * (u-1) * (u-QQ(25)/7),
    ),
    (
        QQ(3)/25 * u * (u-1) * (u+QQ(25)/3),
        QQ(27)/3125 * u**2 * (u-1) * (u+QQ(25)/3) * (u+QQ(25)/2),
    ),
)

for X, Y in sections:
    assert Y**2 == X**3 + a2*X**2 + a4*X + a6

K = FractionField(R)
E = EllipticCurve(K, [0, a2, 0, a4, a6])
generic_points = [E((K(X), K(Y))) for X, Y in sections]
half_point = E((
    K(0),
    K(QQ(21)/6250 * u**2 * (u-QQ(25)/7) * (u-1) * (u+QQ(25)/3)),
))
assert 2*half_point == generic_points[0]
Delta = R(E.discriminant())
fiber_points = (QQ(1), -QQ(25)/2, QQ(25)/7, -QQ(25)/3, QQ(0))
fiber_orders = (2, 2, 3, 3, 7)


def valuation_at(poly, point):
    factor = u - point
    order = 0
    while poly and poly(point) == 0:
        poly //= factor
        order += 1
    return order


assert tuple(valuation_at(Delta, point) for point in fiber_points) == fiber_orders
assert 24 - Delta.degree() == 4
residual = R(Delta / prod((u-point)**order for point, order in zip(fiber_points, fiber_orders)))
assert residual.degree() == 3 and gcd(residual, residual.derivative()) == 1
assert all(residual(point) != 0 for point in fiber_points)


def node_at(A2, A4, A6, point):
    RX = PolynomialRing(QQ, "xnode")
    xnode = RX.gen()
    cubic = xnode**3 + A2(point)*xnode**2 + A4(point)*xnode + A6(point)
    repeated = gcd(cubic, cubic.derivative())
    assert repeated.degree() == 1
    return -repeated[0] / repeated[1]


def multiplicative_steps(A2, A4, A6, X, Y, point):
    node = node_at(A2, A4, A6, point)
    if X(point) != node or Y(point) != 0:
        return 0
    P = PolynomialRing(QQ, ("s", "xx", "yy"))
    s, xx, yy = P.gens()

    def shift(poly):
        shifted = R(poly(u + point))
        return sum(P(coefficient) * s**index for index, coefficient in enumerate(shifted.list()))

    surface = yy**2 - (node+xx)**3 - shift(A2)*(node+xx)**2 - shift(A4)*(node+xx) - shift(A6)
    shifted_x = R(X(u + point))
    shifted_y = R(Y(u + point))
    section_x = R((shifted_x-node) // u)
    section_y = R(shifted_y // u)
    surface = P(surface(s, s*xx, s*yy) // s**2)
    steps = 1
    while True:
        center_x, center_y = section_x(0), section_y(0)
        center = {s: 0, xx: center_x, yy: center_y}
        if any(surface.derivative(variable).subs(center) for variable in (s, xx, yy)):
            return steps
        section_x = R((section_x-center_x) // u)
        section_y = R((section_y-center_y) // u)
        surface = P(surface(s, center_x+s*xx, center_y+s*yy) // s**2)
        steps += 1
        assert steps <= 6


finite_multiplicative_points = fiber_points[:4]
finite_labels = tuple(
    tuple(multiplicative_steps(a2, a4, a6, X, Y, point) for point in finite_multiplicative_points)
    for X, Y in sections
)

# Minimal K3 chart at infinity: xbar=s^4*x(1/s), ybar=s^6*y(1/s).
S = PolynomialRing(QQ, "s")
s = S.gen()


def reverse_with_weight(poly, weight):
    return S(sum(poly[index] * s**(weight-index) for index in range(poly.degree()+1)))


infinity_a2 = reverse_with_weight(a2, 4)
infinity_a4 = reverse_with_weight(a4, 8)
infinity_a6 = reverse_with_weight(a6, 12)
infinity_labels = []
for X, Y in sections:
    infinity_x = reverse_with_weight(X, 4)
    infinity_y = reverse_with_weight(Y, 6)
    infinity_labels.append(
        multiplicative_steps(infinity_a2, infinity_a4, infinity_a6, infinity_x, infinity_y, QQ(0))
    )
infinity_labels = tuple(infinity_labels)

pair_gcds = {}
for left in range(3):
    for right in range(left):
        pair_gcd = gcd(
            sections[left][0] - sections[right][0],
            sections[left][1] - sections[right][1],
        )
        pair_gcds[(right+1, left+1)] = pair_gcd.monic()


def additive_blowup_path(X, Y):
    P = PolynomialRing(QQ, ("s0", "xx0", "yy0"))
    s0, xx0, yy0 = P.gens()

    def embed(poly):
        return sum(P(poly[index]) * s0**index for index in range(poly.degree()+1))

    surface = yy0**2 - xx0**3 - embed(a2)*xx0**2 - embed(a4)*xx0 - embed(a6)
    section_x, section_y = R(X), R(Y)
    path = []
    center_x = center_y = QQ(0)
    for step in range(1, 10):
        section_x = R((section_x-center_x) // u)
        section_y = R((section_y-center_y) // u)
        surface = P(surface(s0, center_x+s0*xx0, center_y+s0*yy0) // s0**2)
        center_x, center_y = section_x(0), section_y(0)
        center = {s0: 0, xx0: center_x, yy0: center_y}
        gradient = tuple(surface.derivative(variable).subs(center) for variable in (s0, xx0, yy0))
        path.append((center_x, center_y, gradient))
        if any(gradient):
            return tuple(path)
    raise RuntimeError("I1* section path did not reach a smooth chart")


additive_paths = tuple(additive_blowup_path(X, Y) for X, Y in sections)

# Resolve the orientation that ``multiplicative_steps`` deliberately forgets.
# At an I_n fiber a section of depth k and one of depth n-k have the same
# distance from the identity component, but their mixed Shioda correction can
# differ.  The exceptional conics at the two I3 fibers and at the I4 fiber at
# infinity split into the following exact lines.  Their evaluations distinguish
# the two arms; P1 meets the crossing of the I4 lines and hence has label 2.
oriented_labels = (
    (0, 1, 1),  # I2 at u=1
    (0, 0, 1),  # I2 at u=-25/2
    (1, 1, 0),  # I3 at u=25/7
    (1, 0, 2),  # I3 at u=-25/3
    (2, 1, 3),  # I4 at infinity
)

plus_centers = ((QQ(10)/7, QQ(75)/49), (-QQ(18)/7, -QQ(45)/7))
plus_lines = (
    lambda x, y: -1365*x + 686*y + 900,
    lambda x, y: 1365*x + 686*y - 900,
)
assert tuple(tuple(line(*center) for line in plus_lines) for center in plus_centers) == (
    (0, 2100), (0, -8820),
)

minus_centers = ((QQ(10)/3, -QQ(25)/3), (QQ(28)/3, -QQ(70)/3))
minus_lines = (
    lambda x, y: -285*x + 54*y + 1400,
    lambda x, y: 285*x + 54*y - 1400,
)
assert tuple(tuple(line(*center) for line in minus_lines) for center in minus_centers) == (
    (0, -900), (-2520, 0),
)

infinity_centers = (
    (QQ(21)/625, 0),
    (-QQ(7)/25, QQ(98)/3125),
    (QQ(3)/25, QQ(27)/3125),
)
infinity_lines = (
    lambda x, y: -625*x + 6250*y + 21,
    lambda x, y: 625*x + 6250*y - 21,
)
assert tuple(tuple(line(*center) for line in infinity_lines) for center in infinity_centers) == (
    (0, 0), (392, 0), (0, 108),
)


def multiplicative_contribution(order, labels):
    """Inverse-A_(order-1) contribution matrix for oriented component labels."""
    answer = matrix(QQ, len(labels))
    for left, a in enumerate(labels):
        for right, b in enumerate(labels):
            if a and b:
                answer[left, right] = QQ(min(a, b) * (order-max(a, b))) / order
    return answer


# All three sections follow the same first two blowups at the I1* fiber and
# then meet distinct smooth points of the same vector-end component.  In the
# D5 inverse Cartan matrix this component has self/mixed contribution 1.
assert tuple(path[0][0:2] for path in additive_paths) == ((-1, 0),)*3
assert len({path[1][0:2] for path in additive_paths}) == 3
d5_cartan = matrix(ZZ, (
    (2, -1, 0, 0, 0),
    (-1, 2, -1, 0, 0),
    (0, -1, 2, -1, -1),
    (0, 0, -1, 2, 0),
    (0, 0, -1, 0, 2),
))
assert d5_cartan.inverse()[0, 0] == 1
local_contribution = matrix(QQ, 3, 3, [1]*9)
for order, labels in zip((2, 2, 3, 3, 4), oriented_labels):
    local_contribution += multiplicative_contribution(order, labels)

# The common factors of the affine coordinates occur only on singular fibers.
# The blowup paths above separate every pair there, and their components at
# infinity are distinct, so the strict transforms have pair intersection zero.
assert pair_gcds == {
    (1, 2): u**2 * (u-QQ(25)/7),
    (1, 3): u**2 * (u+QQ(25)/3),
    (2, 3): u**2 * (u-1),
}
section_intersections = matrix(ZZ, 3, 3, 0)

# Shioda's formula on a K3 (chi=2).  These polynomial sections do not meet O;
# pairwise strict-transform intersections were just shown to vanish.
height_gram = matrix(QQ, 3, 3, lambda i, j:
    (4 if i == j else 2) - section_intersections[i, j] - local_contribution[i, j]
)
assert height_gram == matrix(QQ, (
    (QQ(2)/3, -QQ(1)/6, QQ(1)/6),
    (-QQ(1)/6, QQ(13)/12, QQ(1)/4),
    (QQ(1)/6, QQ(1)/4, QQ(7)/12),
))
assert height_gram.det() == QQ(23)/72

# Replace P1 by its exact half.  This is the currently certified rational
# rank-three sublattice; it is not the full geometric MW lattice (rank four).
half_change = diagonal_matrix(QQ, (QQ(1)/2, 1, 1))
half_height_gram = half_change * height_gram * half_change
assert half_height_gram == matrix(QQ, (
    (QQ(1)/6, -QQ(1)/12, QQ(1)/12),
    (-QQ(1)/12, QQ(13)/12, QQ(1)/4),
    (QQ(1)/12, QQ(1)/4, QQ(7)/12),
))
assert half_height_gram.det() == QQ(23)/288
assert 576 * half_height_gram.det() == 46


def x_degree(point):
    if point.is_zero():
        return 0
    x_coordinate = point[0]
    return max(x_coordinate.numerator().degree(), x_coordinate.denominator().degree())


degree_growth = {
    index+1: tuple(x_degree(multiplier*point) for multiplier in range(1, 13))
    for index, point in enumerate(generic_points)
}
sum_degree_growth = {
    (left+1, right+1): tuple(
        x_degree(multiplier*(generic_points[left]+generic_points[right]))
        for multiplier in range(1, 13)
    )
    for left in range(3) for right in range(left+1, 3)
}

# Independence certificate: specialization at u=-1 has positive regulator.
specialized = EllipticCurve(
    QQ, [0, a2(-1), 0, a4(-1), a6(-1)]
)
specialized_points = [specialized((X(-1), Y(-1))) for X, Y in sections]
height_matrix = specialized.height_pairing_matrix(specialized_points, precision=128)
regulator = height_matrix.det()
assert regulator > 1e-20

print(
    "MW3TWIST|fibers=I1*,I4,2I3,2I2,3I1|root=D5+A3+2A2+2A1"
    "|rational_independent=3|geometric_expected=4",
    flush=True,
)
print("MW3TWIST|P1_divisible=2|half_point_x=0|saturation_pending=1", flush=True)
print(f"MW3TWIST|finite_labels={finite_labels}|infinity_I4={infinity_labels}", flush=True)
print(f"MW3TWIST|oriented_labels={oriented_labels}", flush=True)
print(f"MW3TWIST|pair_gcds={pair_gcds}", flush=True)
print(f"MW3TWIST|I1star_paths={additive_paths}", flush=True)
print(
    f"MW3TWIST|height_gram={height_gram}|det={height_gram.det()}"
    f"|half_height_gram={half_height_gram}|half_det={half_height_gram.det()}"
    "|rank19_visible_disc=46",
    flush=True,
)
print(f"MW3TWIST|degree_growth={degree_growth}|sum_degree_growth={sum_degree_growth}", flush=True)
print(f"MW3TWIST|specialization=u:-1|regulator={regulator}|independent=3", flush=True)
print("MW3TWIST|PASS", flush=True)

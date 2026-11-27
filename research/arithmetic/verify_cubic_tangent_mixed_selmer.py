"""Replay the local mixed isogeny pairing for the marked triple (53,167,281)."""
from itertools import product

from sage.all import EllipticCurve, PolynomialRing, Qp, Qq, ZZ


P = 167
H = 114
D = 10179717
B = ZZ((H // 2) ** 2 * D)
X0 = P * (P - H)
Y0 = (2 * P**3 - 3 * P**2 * H - H**3) // 2
X1 = P**2 - H**2
Y1 = -P**3 + 3 * P * H**2
PRECISION = 48


def preimage_coordinates(b, x, y, s, root):
    """Return the dual-isogeny preimage with Cardano cube root root."""
    assert (root**3 - (y + s)).valuation() >= 24
    x_image = 3 * (x + root**2 + (x / root) ** 2)
    y_image = 27 * y * x_image**3 / (x_image**3 + 216 * b)
    assert (y_image**2 - x_image**3 + 27 * b).valuation() >= 18
    assert ((x_image**3 - 108 * b) / (9 * x_image**2) - x).valuation() >= 12
    return x_image, y_image


def has_tame_cube_root(a, q):
    field = a.parent()
    start = next(field(i) for i in range(q) if field(i).residue() ** 3 == a.residue())
    for _ in range(7):
        start -= (start**3 - a) / (3 * start**2)
    assert (start**3 - a).valuation() >= 24
    return start


def tame_term(q):
    base = Qp(q, prec=PRECISION)
    if q == 2:
        t = PolynomialRing(base, "t").gen()
        field = base.extension(t**2 - t + (1 - D) // 4, names="omega")
        sqrt_d = 2 * field.gen() - 1
    else:
        field = base
        sqrt_d = base(D).sqrt()
        if sqrt_d.residue() != 7:
            sqrt_d = -sqrt_d
    alpha = field(Y0) + (H // 2) * sqrt_d
    beta = (sqrt_d - P**2 - H**2) / 2
    delta = -(P - H) * beta / alpha
    sqrt_minus_three = field(-3).sqrt()
    if q == 19 and sqrt_minus_three.residue() != 4:
        sqrt_minus_three = -sqrt_minus_three
    zeta = (sqrt_minus_three - 1) / 2
    generator = zeta.residue()
    root = has_tame_cube_root(alpha, q)
    values = []
    for j, (x, y, theta) in enumerate(((X0, Y0, root), (X1, Y1, root * delta))):
        s = 2**j * (H // 2) * sqrt_d
        _, image_y = preimage_coordinates(4**j * B, field(x), field(y), s, theta)
        values.append(image_y + 3 * sqrt_minus_three * s)
    u0, u1 = values
    v0, v1 = (int(u.valuation()) for u in values)
    residue_power = (generator.parent().cardinality() - 1) // 3
    normalized = tuple((u / field(q)**v).residue() ** residue_power
                       for u, v in zip(values, (v0, v1)))

    def index(root_of_unity):
        return next(i for i in range(3) if generator**i == root_of_unity)

    if q == 2:
        assert (v0, v1) == (0, 4)
        epsilon_squared = field.gen().residue() ** 2
        return 2 * index(epsilon_squared ** (-v0 - v1)) % 3

    assert q == 19 and (v0, v1) == (2, 0) and tuple(map(int, normalized)) == (7, 1)
    gamma_pair = normalized[1] / field(7).residue() ** v1
    epsilon_pair = field(11).residue() ** (-v0 - v1)
    return (index(epsilon_pair) - index(gamma_pair)) % 3


def over_base(element, base, minimum_precision=12):
    coefficients = element.polynomial().list()
    assert all(coefficient.valuation() >= minimum_precision for coefficient in coefficients[1:])
    return base(coefficients[0])


def cube_root_over_ramified_quadratic(a, uniformizer):
    """Use v(3)=2. Starting at precision above four certifies Hensel lifting."""
    field = a.parent()
    valuation = int(a.valuation())
    assert valuation % 3 == 0
    if valuation:
        return uniformizer ** (valuation // 3) * cube_root_over_ramified_quadratic(
            a / uniformizer**valuation, uniformizer)
    residue = field.residue_field()
    candidates = [field.teichmuller(a.residue() ** 9)]
    for depth in range(1, 5):
        next_candidates = []
        for candidate in candidates:
            for digit in residue:
                updated = candidate + field.teichmuller(digit) * uniformizer**depth
                if (updated**3 - a).valuation() >= depth + 3:
                    next_candidates.append(updated)
        candidates = next_candidates
        for candidate in candidates:
            if (candidate**3 - a).valuation() > 4:
                for _ in range(7):
                    candidate -= (candidate**3 - a) / (3 * candidate**2)
                assert (candidate**3 - a).valuation() >= 30
                return candidate
    raise AssertionError("No certified local cube root")


def cube_in_unramified_cubic(a):
    valuation = int(a.valuation())
    assert a.precision_absolute() >= valuation + 2
    if valuation % 3:
        return False
    unit = a / a.parent()(3)**valuation
    return (unit**26 - 1).valuation() >= 2


def wild_term():
    base = Qp(3, prec=PRECISION)
    unramified = Qq(ZZ(27), names="mu", prec=PRECISION)
    t = PolynomialRing(unramified, "t").gen()
    large = unramified.extension(t**6 - 3*t**4 + 3*t**2 + 3, names="v")
    v = large.gen()
    pi = (v**2 - 3) / (v**2 - 1)
    r = pi - 1
    zeta = (v * pi - 1) / 2
    quadratic = unramified.extension(t**2 + 3, names="lambda")
    lam = quadratic.gen()
    w = (-base(B) / 27).sqrt()
    if w.residue() != 1:
        w = -w

    original = EllipticCurve(large, [0, B])
    reflected = EllipticCurve(large, [0, -27 * B])
    e = [EllipticCurve(unramified, [0, 4**j * B]) for j in range(3)]
    ep = [EllipticCurve(unramified, [0, -27 * 4**j * B]) for j in range(3)]

    def lifted(j, point, root=None):
        x, y = map(quadratic, (point[0], point[1]))
        s = 3 * 2**j * lam * w
        if root is None:
            root = cube_root_over_ramified_quadratic(y + s, lam)
        x, y = preimage_coordinates(4**j * B, x, y, s, root)
        return ep[j](over_base(x, unramified), over_base(y, unramified))

    def kummer(j, point):
        s = -27 * 2**j * w
        if point.is_zero():
            return unramified(1)
        if point[0] == 0:
            return unramified(2 * s) if (point[1] - s).valuation() >= 15 else unramified(1 / (2 * s))
        return point[1] + s

    conjugates = []
    for i in range(3):
        pi_i = 1 + zeta**i * r
        square_root = (1 + B / pi_i**6).sqrt()
        if square_root.residue() != 1:
            square_root = -square_root
        conjugates.append(-original(pi_i**2, pi_i**3 * square_root))

    trace_points = []
    for j in range(3):
        trace_point = original(0)
        for i, point in enumerate(conjugates):
            trace_point += original(zeta**(-i*j) * point[0], point[1])
        trace_points.append(e[j](over_base(r**(2*j) * trace_point[0], unramified),
                                 over_base(2**j * trace_point[1], unramified)))
    base_lifts = [lifted(j, point) for j, point in enumerate(trace_points)]
    torsion = [ep[j](0, -27 * 2**j * w) for j in range(3)]
    x, y = conjugates[0][0:2]
    target = reflected((x**3 + 4*B) / x**2, y * (x**3 - 8*B) / x**3)
    target_shift = reflected(0, -27*w)
    assert all(((target + sign * target_shift)[0] - target[0]).valuation() < 15
               for sign in (-1, 1))

    compatible_lifts = []
    for offsets in product(range(3), repeat=3):
        lifts = [base_lifts[j] + offsets[j] * torsion[j] for j in range(3)]
        row = reflected(0)
        for j, point in enumerate(lifts):
            row += reflected(r**(-2*j) * large(point[0]), large(2)**(-j) * large(point[1]))
        if row.is_zero() or (row[0] - target[0]).valuation() < 15:
            continue
        if (row[1] - target[1]).valuation() < 15:
            continue
        if cube_in_unramified_cubic(kummer(0, lifts[0]) * kummer(1, lifts[1]) * kummer(2, lifts[2])):
            compatible_lifts.append(lifts)
    assert len(compatible_lifts) == 3

    alpha = Y0 + 3 * lam * w
    beta = (3 * lam * w / (H // 2) - P**2 - H**2) / 2
    delta = -(P - H) * beta / alpha
    root = cube_root_over_ramified_quadratic(alpha, lam)
    rational_lift = lifted(1, e[1](X1, Y1), root * delta)
    assert all(cube_in_unramified_cubic(kummer(1, rational_lift - lifts[1]))
               for lifts in compatible_lifts)
    sample = kummer(1, rational_lift - compatible_lifts[0][1])
    assert sample.valuation() == 3 and (sample / 27 + 1).valuation() >= 2
    return 0


if __name__ == "__main__":
    contributions = (tame_term(2), wild_term(), tame_term(19))
    assert contributions == (2, 0, 1)
    assert sum(contributions) % 3 == 0
    print("Local contributions at 2, 3, 19 in thirds:", contributions)

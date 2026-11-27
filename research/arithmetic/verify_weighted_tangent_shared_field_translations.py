"""Replay the two profiles, their prime-form cell and an exact same-field pair."""

from math import gcd

from sage.all import EllipticCurve, GF, QQ, ZZ, PolynomialRing, QuadraticField


SOURCE_ROWS = [(1, 0), (1, -3), (1, -2)]
TARGET_ROWS = [(1186923, -1978205), (4226103, -14489496), (892143, -764694)]
C0 = 13209
G = C0**2
K1 = QQ(-4913) / 55566


def homogeneous_quartic(k, u, h):
    w = k * (k + 1)
    return -3 * u**4 + 4 * (2 * k + 1) * u**3 * h - 6 * w * u**2 * h**2 + w**2 * h**4


def forms(rows, u, h):
    return tuple(a * u + b * h for a, b in rows)


def check_parameter_map():
    curve = EllipticCurve(QQ, [0, 576])
    center = curve(24, -120)
    torsion = curve(0, 24)
    image_center = 3 * center
    assert image_center.xy() == (QQ(-2516) / 441, QQ(-182960) / 9261)
    scale = image_center[0] / 24
    assert (-image_center[1] / 24 - 1) / 2 == K1

    for parameter, point in [(0, curve(0)), (2, -torsion), (3, torsion)]:
        translated = point + center
        transformed = 12 * scale**2 * (translated[0] - image_center[0]) / (
            translated[1] + image_center[1]
        )
        assert transformed == QQ(1186923 * parameter - 1978205) / (
            -3333960 * parameter + 13724802
        )

    ring = PolynomialRing(QQ, "u,h")
    u, h = ring.gens()
    u1 = (1186923 * u - 1978205 * h) / G
    h1 = -55566 * (60 * u - 247 * h) / G
    assert forms(TARGET_ROWS, u, h) == tuple(
        G * value for value in (u1, u1 - (K1 + 1) * h1, u1 - K1 * h1)
    )
    assert homogeneous_quartic(K1, u1, h1) == homogeneous_quartic(2, u, h) / G
    polynomial_ring = PolynomialRing(GF(5), "t")
    assert homogeneous_quartic(2, polynomial_ring.gen(), 1).is_irreducible()

    combined_rows = SOURCE_ROWS + TARGET_ROWS
    for index, (a, b) in enumerate(combined_rows):
        for c, d in combined_rows[index + 1 :]:
            assert a * d - b * c != 0
    for row in SOURCE_ROWS:
        assert row[0] * (-247) - row[1] * 60 != 0
    for row in TARGET_ROWS:
        assert row[0] != 0
    assert [gcd(a, b) for a, b in TARGET_ROWS] == [17**2 * 37**2, 3**2 * 7**2 * 37**2, 3**2 * 7**2 * 17**2]


def check_residues():
    table = [
        (2, 1, 0, 0, 0),
        (3, 1, 15, 2, 3),
        (5, 1, 1, 0, 0),
        (7, 18, 1, 2, 3),
        (17, 86, 1, 2, 2),
        (37, 392, 1, 2, 2),
    ]
    for q, u, h, target_v, quartic_v in table:
        assert all(value % q != 0 for value in forms(SOURCE_ROWS, u, h))
        assert all(ZZ(value).valuation(q) == target_v for value in forms(TARGET_ROWS, u, h))
        assert ZZ(homogeneous_quartic(2, u, h)).valuation(q) == quartic_v

    cross_values = [
        (QQ(2), QQ(247) / 60, QQ(-6311) / 3600, QQ(4403) / 1200),
        (K1, QQ(-395641) / 1111320, QQ(278135623) / 400075200, QQ(248858189) / 1200225600),
    ]
    for k, t, b, v in cross_values:
        assert (t**2 - 2 * (k + 1) * t + k * (k + 1)) == b
        assert homogeneous_quartic(k, t, 1) == -3 * v**2
        assert v * (b**2 - v**2) != 0


def reconstruct(primes, profile, d, delta, square_root):
    p1, p2, p3 = primes
    a,b = profile.denominator(), profile.numerator()
    c = -a - b
    x, y, z = abs(c), a, abs(b)
    A, B, C = (y**2 * p1**3, z**2 * p2**3, x**2 * p3**3)
    assert (A - B) % x == (C - B) % y == (A - C) % z == 0
    heron = (4 * A * C - (A + C - B)**2) / QQ((x * y * z)**2)
    assert heron == -d * delta**2
    assert all(ZZ(heron).gcd(p) == 1 for p in primes)
    generators = (
        ((A + B - C) / QQ(y * z) + delta * x * square_root) / 2,
        ((B + C - A) / QQ(x * z) + delta * y * square_root) / 2,
        ((B - A - C) / QQ(x * y) + delta * z * square_root) / 2,
    )
    for element, norm in zip(generators, (p1 * p2, p2 * p3, p1 * p3)):
        assert element.norm() == norm**3
    assert generators[0] * generators[1] == p2**3 * generators[2]

    step = (p1 - p3) / profile
    parameter = p1 / step
    w = profile * (profile + 1)
    sgn = 1 if w > 0 else -1
    U = delta * a * square_root / (sgn * step**3)
    assert U**2 == homogeneous_quartic(profile, parameter, 1)
    beta = step**2 * (parameter**2 - 2 * (profile + 1) * parameter + w - U) / 2
    assert beta**3 == generators[0] * generators[1].conjugate()
    return generators[0], generators[0] * generators[1].conjugate()


def check_example():
    u, h = ZZ(8053580009), ZZ(-8202927966)
    assert gcd(u, h) == 1
    residues = [(2, 1, 0, 1), (3, 1, 15, 4), (5, 1, 1, 1), (7, 18, 1, 4), (17, 86, 1, 3), (37, 392, 1, 3)]
    for q, a, b, exponent in residues:
        assert (u * b - h * a) % q**exponent == 0
        assert gcd(gcd(u, h), q) == 1
    source = forms(SOURCE_ROWS, u, h)
    target = tuple(value // G for value in forms(TARGET_ROWS, u, h))
    assert source == (8053580009, 32662363907, 24459435941)
    assert target == (147789977, 876281423, 77131211)
    assert len(set(source + target)) == 6
    assert all(q.is_prime(proof=True) for q in source + target)
    d = homogeneous_quartic(2, u, h) // G
    assert d == -529788943876447727940522414336123 and d % 4 == 1 and d.is_squarefree()
    field = QuadraticField(d, "s")
    f1, f2 = ZZ(108352475502894), ZZ(14382)
    assert f1 % 23 == f2 % 47 == 0 and f1 % 47 != 0 and f2 % 23 != 0
    first, theta1 = reconstruct(source, QQ(2), d, f1, field.gen())
    second, theta2 = reconstruct(target, K1, d, f2, field.gen())
    assert (first / second).is_nth_power(3) and not first.is_nth_power(3)
    assert QQ(theta1[1]) % 47 == 8 and QQ(theta2[1]) % 23 == 22
    print("Verified the shared field and norm class, and both failures in the common order.")



if __name__ == "__main__":
    check_parameter_map()
    check_residues()
    check_example()

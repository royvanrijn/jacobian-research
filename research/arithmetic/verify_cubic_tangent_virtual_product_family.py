"""Verify the odd inverse profile pair and its actual order-three ideal class."""

from sage.all import EllipticCurve, QQ, QuadraticField


def check_profile_points():
    field = QuadraticField(-31, "a")
    a = field.gen()
    b = 35
    m = b**2
    s = -QQ(1) / 7
    w = m * s**3
    assert w == -QQ(25) / 7
    Y = 48 * a / 49
    assert Y**2 == (w**2 - 1) * (w**2 + 4 * w - 1) == -QQ(71424) / 2401
    x = 2 * m * s**2 / (w - 1)
    z = -2 * m * s / (w + 1)
    y = b * Y / (w - 1)**2
    v = m * Y / (w + 1)**2
    profiles = [
        (b, (x, y), QQ(1) / 3, (-10, -15), QQ(24) / 7, -QQ(135) / 7, QQ(720) / 77),
        (m, (z, v), QQ(16) / 17, (-49, -1176), QQ(561) / 6400, -QQ(240737) / 204800,
         QQ(15421329) / 78438400),
    ]
    for coefficient, point, r, half, scale, k, expected in profiles:
        assert k == -(2 * r + 1) / (r**3 * (r + 2))
        weight, middle = k * (k + 1), 2 * k + 1
        assert coefficient * scale**3 == 4 * weight
        curve = EllipticCurve(field, [0, coefficient**2])
        center = 2 * curve(*half)
        assert scale**2 * center[0] == 4 * weight
        assert scale**3 * center[1] == -4 * weight * middle
        translated = curve(*point) + curve(*half)
        value = 2 * weight * (scale**2 * translated[0] - 4 * weight) / (
            scale**3 * translated[1] - 4 * weight * middle)
        assert value == expected and min(value, value - k, value - k - 1) > 0

    alpha, beta = y + b, v + m
    gamma = Y / (w - 1)**2 + (w + 3) / (w - 1)
    assert gamma == (8 + 3 * a) / 64
    assert gamma.norm() == (QQ(7) / 16)**3
    assert not gamma.is_nth_power(3)
    assert (alpha * beta * gamma).is_nth_power(3)
    assert sorted(gamma.valuation(p) for p in field.primes_above(2)) == [-6, -6]
    assert sorted(gamma.valuation(p) for p in field.primes_above(7)) == [0, 3]
    group = field.class_group(proof=True)
    assert group.invariants() == (3,)
    marked_prime = next(p for p in field.primes_above(7) if gamma.valuation(p) == 3)
    assert group(marked_prime).order() == 3


def check_conductors():
    for coefficient in (35, 1225):
        curve = EllipticCurve(QQ, [0, coefficient**2])
        assert curve.conductor() == 44100
        assert curve.quadratic_twist(-31).conductor() == 44100 * 31**2


if __name__ == "__main__":
    check_profile_points()
    check_conductors()
    print("Verified both positive profiles, the order-three ideal class over Q(sqrt(-31)), and both conductors.")

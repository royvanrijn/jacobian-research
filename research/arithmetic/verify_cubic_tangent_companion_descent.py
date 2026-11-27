"""With Sage's Python, replay the companion points, quadratic ray data and rank bounds."""

from sage.all import EllipticCurve, QQ, ZZ, QuadraticField, pari


def check_fixture(p, h, expected):
    center, step = ZZ(p), ZZ(h)
    assert all(q.is_prime(proof=True) for q in (center - step, center, center + step))
    d = step**4 + 6 * center**2 * step**2 - 3 * center**4
    assert d.is_squarefree()
    m = -d // 3
    field = QuadraticField(d, "s")
    s = field.gen()
    base = EllipticCurve(field, [0, 1])
    torsion = base(-1, 0)
    t = QQ(center) / step
    u = s / step**2
    x = (u + 1 + t**2) / (2 * t**2)
    point = base(x, (x + 1) / t)
    q = 2 * point - torsion
    trace = 2 * center**3 - 3 * center**2 * step - step**3
    alpha = (trace + step * s) / 2
    assert alpha == (point[1] + 1) * (step * (1 - t) / x)**3
    assert not alpha.is_nth_power(3)
    assert d * (q[0] + 1) == 12 * (center * step)**2
    companion = EllipticCurve(QQ, [0, d**3])
    rational_point = companion(QQ(d * q[0]), QQ(d * s * q[1]))
    assert ((rational_point[1] + d * s) / alpha**2).is_nth_power(3)

    variable = pari("x")
    polynomial = variable**2 - variable + (1 - d) // 4
    mirror_polynomial = variable**2 - variable + (1 - m) // 4
    ordinary = polynomial.bnfinit(1)
    mirror = mirror_polynomial.bnfinit(1)
    assert ordinary.bnfcertify() == mirror.bnfcertify() == 1
    ray = ordinary.bnrinit(6, 1)
    prime_over_2 = mirror.idealprimedec(2)[0]
    mirror_generators = tuple(map(int, mirror.bnf_get_cyc()))
    dyadic_class = tuple(map(int, mirror.bnfisprincipal(prime_over_2, 0)))
    ray_generators = tuple(map(int, pari("r->r.cyc")(ray)))
    first_dim = sum(n % 3 == 0 for n in mirror_generators) - any(
        n % 3 == 0 and a % 3 != 0 for n, a in zip(mirror_generators, dyadic_class)
    )
    second_dim = sum(n % 3 == 0 for n in ray_generators)

    norm_ideal = pari.matid(2)
    generator = (pari((trace - step) // 2) + step * variable).Mod(polynomial)
    for ell in (center - step, center):
        factors = ordinary.idealprimedec(ell)
        assert sorted(ordinary.idealval(generator, prime) for prime in factors) == [0, 3]
        for prime in factors:
            if ordinary.idealval(generator, prime) == 3:
                norm_ideal = ordinary.idealmul(norm_ideal, prime)
    edge_class = tuple(map(int, ordinary.bnfisprincipal(norm_ideal, 0)))

    lower, upper = tuple(pari(companion.minimal_model()).ellrank()[:2])
    assert lower == upper
    sha_dimension = second_dim - lower + 1
    actual = (d, lower, first_dim, second_dim, sha_dimension, any(edge_class))
    assert actual == expected
    print((p, h), "field, rank, Selmer dimensions, Sha[3] dimension, nonzero edge class:", actual)


def main():
    base_twist = pari(EllipticCurve(QQ, [0, -27]))
    assert tuple(base_twist.elllocalred(ell)[0] for ell in (2, 3)) == (2, 2)
    cases = [
        (17, 6, (-186843, 4, 1, 3, 0, True)),
        (101, 30, (-256285803, 2, 1, 3, 2, True)),
        (167, 114, (10179717, 1, 1, 2, 2, True)),
        (37, 30, (2580117, 1, 1, 2, 2, False)),
    ]
    for p, h, expected in cases:
        check_fixture(p, h, expected)


if __name__ == "__main__":
    main()

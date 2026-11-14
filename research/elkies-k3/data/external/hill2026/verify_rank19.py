"""Independent exact replay of the paired-cover rank-19 certificate.

Run with Sage's Python from this directory: sage -python verify.py
The only nonstandard dependency is SageMath. All inputs are in data/.
"""
import json
from math import gcd, lcm
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / 'data' / name).read_text())


def model(name):
    data = read(name)
    ring = PolynomialRing(QQ, 'U')
    field = ring.fraction_field()
    a4 = ring(data['jacobian_model']['a4_low_to_high'])
    a6 = ring(data['jacobian_model']['a6_low_to_high'])
    curve = EllipticCurve(field, [0, 0, 0, a4, a6])
    points = []
    for row in data['neighbor_sections']:
        x = field(ring(row['x_numerator'].replace('^', '**'))) / ring(row['x_denominator'].replace('^', '**'))
        y = field(ring(row['y_numerator'].replace('^', '**'))) / ring(row['y_denominator'].replace('^', '**'))
        points.append(curve(x, y))  # Also checks each defining identity.
    return data, ring, field, curve, points


def exact_height(point, chi):
    if not point:
        return QQ(0)
    x = point[0]
    return QQ(max(x.numerator().degree(), x.denominator().degree() + 2 * chi))


def check_gram(data, ring, curve, points, chi):
    delta = ring(curve.discriminant())
    assert ring(curve.a4()).degree() <= 4 * chi
    assert ring(curve.a6()).degree() <= 6 * chi
    assert delta.degree() == 12 * chi
    assert delta.gcd(delta.derivative()).degree() == 0
    heights = [exact_height(point, chi) for point in points]
    gram = matrix(QQ, len(points))
    for i in range(len(points)):
        gram[i, i] = heights[i]
        for j in range(i):
            gram[i, j] = gram[j, i] = (
                exact_height(points[i] + points[j], chi)
                - heights[i] - heights[j]
            ) / 2
    assert gram == matrix(QQ, data['exact_generic_section_basis']['height_matrix'])
    assert gram.is_positive_definite()
    return gram


def compose(value, argument):
    return value.numerator()(argument) / value.denominator()(argument)


def main():
    pair = read('pair.json')
    parent, parent_ring, parent_field, parent_curve, inherited = model('parent-rank17.json')
    assert len(inherited) == 17
    parent_gram = check_gram(parent, parent_ring, parent_curve, inherited, 2)
    assert parent_gram.det() == 1152
    print('Parent: 17 independent sections; exact Gram determinant 1152')

    t_ring = PolynomialRing(QQ, 't')
    first, second = [read(f'cover-{index}.json') for index in (4339, 36773)]
    q1, q2 = [t_ring(d['quadratic_basechange']['branch_coefficients'])
              for d in (first, second)]
    assert all(q.degree() == 2 and q.gcd(q.derivative()).degree() == 0
               for q in (q1, q2))
    assert q1.gcd(q2).degree() == 0

    n = t_ring(pair['parent_numerator'])
    d = t_ring(pair['parent_denominator'])
    w = t_ring(pair['first_root_numerator'])
    quartic = t_ring(pair['quartic'])
    assert q1[2]*n*n + q1[1]*n*d + q1[0]*d*d == w*w
    assert q2[2]*n*n + q2[1]*n*d + q2[0]*d*d == quartic
    assert quartic.degree() == 4
    assert quartic.gcd(quartic.derivative()).degree() == 0
    assert quartic[0].is_square()
    print('Biquadratic cover: exact conic and genus-one quartic identities hold')

    # The quartic w^2 = a*u^4+b*u^3+c*u^2+d*u+e^2 has this Jacobian.
    e = quartic[0].sqrt()
    q_d, c, b, a = [quartic[i] for i in range(1, 5)]
    e0 = EllipticCurve(QQ, [0, c, 0,
        b*q_d - 4*a*e*e,
        e*e*b*b + a*q_d*q_d - 4*a*c*e*e])
    assert list(e0.a_invariants()) == list(map(QQ, pair['auxiliary_a_invariants']))
    # Replay the birational map, not just its stated coefficients.
    function_field = t_ring.fraction_field()
    polynomial_over_field = PolynomialRing(function_field, 'w')
    ww = polynomial_over_field.gen()
    u = t_ring.gen()
    xx = (2*e*(ww+e) + q_d*u)/(u*u)
    yy = ((xx*xx - 4*e*e*a)*u - (xx*q_d + 2*e*e*b))/(2*e)
    relation = yy*yy - (xx**3+c*xx*xx+e0.a4()*xx+e0.a6())
    assert relation % (ww*ww - polynomial_over_field(quartic)) == 0

    p0 = e0(*pair['auxiliary_point'])
    torsion_bound = 0
    for prime, claimed in pair['reduction_orders']:
        assert e0.discriminant().numerator() % prime != 0
        actual = int(e0.change_ring(GF(prime)).cardinality())
        assert actual == claimed
        torsion_bound = gcd(torsion_bound, actual)
    assert torsion_bound == 2 and bool(2*p0)
    print('Parameter curve: rational point P has infinite order; torsion bound 2')

    # The recorded P maps to the seed t=15122/33609.
    x0, y0 = p0.xy()
    u0 = (2*e*y0 + x0*q_d + 2*e*e*b)/(x0*x0 - 4*e*e*a)
    seed = QQ(15122)/33609
    assert u0 == 0 and n(u0)/d(u0) == seed

    for index, lift in ((4339, first), (36773, second)):
        _, ring, field, curve, sections = model(f'cover-{index}.json')
        assert len(sections) == 18
        nn = ring(lift['quadratic_basechange']['parameter_numerator'])
        dd = ring(lift['quadratic_basechange']['parameter_denominator'])
        tt = field(nn)/dd
        assert ring(curve.a4()) == ring(dd**8 * compose(parent_curve.a4(), tt))
        assert ring(curve.a6()) == ring(dd**12 * compose(parent_curve.a6(), tt))
        for old, pulled in zip(inherited, sections[:17]):
            assert pulled[0] == dd**4 * compose(old[0], tt)
            assert pulled[1] == dd**6 * compose(old[1], tt)
        gram = check_gram(lift, ring, curve, sections, 4)
        assert gram[:17, :17] == 2*parent_gram
        assert gram.det() == 452984832

        # Nontrivial deck conjugation: the other root of N(U)-tD(U)=0.
        z = field.gen()
        conjugate_z = (nn[0]-tt*dd[0])/((nn[2]-tt*dd[2])*z)
        assert nn(conjugate_z)/dd(conjugate_z) == tt
        scale = dd(z)/dd(conjugate_z)
        def conjugate(point):
            return curve(compose(point[0], conjugate_z)*scale**4,
                         compose(point[1], conjugate_z)*scale**6)
        assert all(conjugate(point) == point for point in sections[:17])
        difference = sections[17] - conjugate(sections[17])
        assert difference and exact_height(difference, 4) == 12
        print(f'Cover {index}: 18 independent sections; Gram determinant {gram.det()}; '
              'anti-invariant height 12')

    # Let sigma_1 act only on sqrt(q1) and sigma_2 only on sqrt(q2).
    # Applying sigma_i-1 to a relation among the 19 sections kills the
    # parent and the other new section. Each remaining difference is
    # nontorsion by the exact positive height above.
    print('Compositum: 19 independent sections over Q(E0)')

    # Trivial torsion for infinitely many specializations. The auxiliary
    # point P reduces at all three primes; n=1 (mod 1440) fixes its images.
    fiber_a4 = QQ(parent_curve.a4()(seed))
    fiber_a6 = QQ(parent_curve.a6()(seed))
    orders = []
    auxiliary_orders = []
    for prime, expected_fiber, expected_aux in ((31, 42, 30), (37, 46, 32), (59, 71, 72)):
        assert e0.discriminant().numerator() % prime != 0
        assert p0[0].denominator() % prime != 0
        assert p0[1].denominator() % prime != 0
        assert GF(prime)(x0*x0-4*e*e*a) != 0
        assert GF(prime)(d[0]) != 0
        assert seed.denominator() % prime != 0
        assert fiber_a4.denominator() % prime != 0
        assert fiber_a6.denominator() % prime != 0
        reduced = EllipticCurve(GF(prime),
            [0, 0, 0, GF(prime)(fiber_a4), GF(prime)(fiber_a6)])
        assert reduced.discriminant() != 0
        fiber_order = int(reduced.cardinality())
        aux_order = int(e0.change_ring(GF(prime)).cardinality())
        assert (fiber_order, aux_order) == (expected_fiber, expected_aux)
        orders.append(fiber_order)
        auxiliary_orders.append(aux_order)
    assert gcd(*orders) == 1
    assert lcm(*auxiliary_orders) == 1440
    print('Torsion: for n=1 mod 1440, fiber orders mod 31,37,59 are '
          '42,46,71; gcd 1')
    print('CERTIFIED: infinitely many trivial-torsion curves have rank at least 19 '
          '(by specialization)')


if __name__ == '__main__':
    main()

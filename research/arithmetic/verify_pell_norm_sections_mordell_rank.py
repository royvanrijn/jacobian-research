"""Exact trace, conductor, and prime controls for PM1–PM5."""

from sage.all import GF, QQ, ZZ, EllipticCurve, PolynomialRing, matrix, pari

from verify_quartic_pell_two_class_directions import F, H, P, Z, t


rank_ring = PolynomialRing(ZZ, "X")
X = rank_ring.gen()


def rank_characteristic(polynomial, prime, expected_traces, infinity_good):
    assert polynomial.change_ring(GF(prime)).is_squarefree()
    degree = 2*polynomial.degree()+(0 if infinity_good else 2)-4
    half = degree//2
    assert half == len(expected_traces)
    trace_sums = []
    for d in range(1, half+1):
        field = GF(prime**d, "a")
        generator = field.multiplicative_generator()
        sextic_exponent = (prime**d-1)//6
        roots = [generator**(sextic_exponent*j) for j in range(6)]
        counts = dict.fromkeys(roots, 0)
        point_traces = [
            prime**d+1-EllipticCurve(field, [0, generator**j]).cardinality()
            for j in range(6)
        ]
        reduced = polynomial.change_ring(field)
        for value in field:
            twist = reduced(value)
            if twist:
                counts[twist**sextic_exponent] += 1
        if infinity_good:
            counts[reduced.leading_coefficient()**sextic_exponent] += 1
        trace_sums.append(-sum(counts[root]*a for root, a in zip(roots, point_traces)))
    assert trace_sums == expected_traces

    coefficients = [ZZ(1)]
    for j in range(1, half+1):
        numerator = -sum(coefficients[j-i]*trace_sums[i-1] for i in range(1, j+1))
        assert numerator % j == 0
        coefficients.append(numerator//j)
    assert coefficients[-1] != 0  # Forces positive reciprocity sign.
    return (
        sum(coefficients[i]*X**(degree-i) for i in range(half+1))
        + sum(coefficients[i]*prime**(degree-2*i)*X**i for i in range(half))
    )


def verify_quartic_mordell():
    assert P**2 == (-1)**3+F
    assert (P+Z)**2 == H**3+F
    assert (P+Z(-t))**2 == H(-t)**3+F
    polynomial = rank_characteristic(
        F, 7, [31, 177, 1126, 15845, 143491, 1124502, 2706595], False
    )
    expected = (X-7)**6*(X**2+7*X+49)*(
        X**6+4*X**5+42*X**4+686*X**3+2058*X**2+9604*X+117649
    )
    assert polynomial == expected and (polynomial//(X-7)**6)(7) != 0


def verify_prime_mordell():
    B, c, auxiliary = ZZ(27), ZZ(511), ZZ(3449)
    A = 2*(B*c-1)
    assert B*c-1 == 4*auxiliary
    assert auxiliary.is_prime(proof=True) and auxiliary % 3 == 2
    assert c % 27 == 25
    exponent = 3
    parameter = (t**exponent-c*c)/A
    consistent = parameter*(B*B*parameter+2)
    trace = 1+B*B*parameter
    assert trace**2-B*B*consistent == 1
    assert (trace/B)**2 == QQ(1)/9**3+consistent
    assert (B*parameter+c)**2 == t**3+consistent

    q = ZZ(23247793)
    assert q.is_prime(proof=True) and q % 9 == 1
    n = ZZ(parameter(q))
    d = ZZ(consistent(q))
    assert n % 32 == 18 and d.valuation(2) == 3 and (d//8) % 8 == 1
    assert d.valuation(3) == 1 and (d//3) % 9 == 1
    assert trace(q) % 5 == 3
    assert 3 not in {(4*x**3-3*x) % 5 for x in range(5)}
    residue = GF(q)
    constant = QQ(c*(B*c-2))/(B*c-1)
    assert constant.valuation(auxiliary) == -1
    assert 2*residue(B*n+c) == residue(constant)
    assert (residue(constant)**((q-1)//3)).multiplicative_order() == 3

    polynomial = rank_characteristic(consistent, 19, [38, 722, 54872, 260642], True)
    assert polynomial == (X-19)**4*(X**2+19*X+361)**2
    assert (polynomial//(X-19)**4)(19) != 0


def verify_conductor_exponents():
    for constant, p, exponent in [
        (1, 3, 2), (2, 2, 6), (10, 2, 6), (8, 2, 6), (3, 3, 5),
    ]:
        curve = EllipticCurve([0, constant])
        assert int(pari(curve).elllocalred(p)[0]) == exponent


def x_degrees(point):
    coordinate = point[0]
    return coordinate.numerator().degree(), coordinate.denominator().degree()


def verify_generic_lattices():
    function_field = t.parent().fraction_field()
    curve = EllipticCurve(function_field, [0, function_field(F)])
    points = [curve(-1, P), curve(H, P+Z), curve(H(-t), P+Z(-t))]
    pairs = [points[0]+points[1], points[0]+points[2], points[1]+points[2]]
    single_degrees = [x_degrees(3*point) for point in points]
    sum_degrees = [x_degrees(3*point) for point in pairs]
    assert single_degrees == [(24, 16), (24, 20), (24, 20)]
    assert sum_degrees == [(60, 56), (60, 56), (24, 20)]
    heights = [QQ(4+max(d, n-4))/9 for n, d in single_degrees+sum_degrees]
    assert heights == [QQ(8)/3]*3+[QQ(20)/3, QQ(20)/3, QQ(8)/3]
    gram = matrix(QQ, 3)
    for i in range(3):
        gram[i, i] = heights[i]
    for (i, j), combined in zip(((0, 1), (0, 2), (1, 2)), heights[3:]):
        gram[i, j] = gram[j, i] = (combined-heights[i]-heights[j])/2
    assert gram == matrix(QQ, [[8, 2, 2], [2, 8, -4], [2, -4, 8]])/3
    assert gram.det() == QQ(32)/3

    rows = []
    for prime, parameter, root, row in [
        (5, 0, 2, [1, 1, 1]),
        (11, 2, 8, [1, 1, 0]),
        (11, 3, 3, [1, 0, 0]),
    ]:
        field = GF(prime)
        theta = field(root)
        assert theta != 0 and theta**3 == -field(F(parameter))
        coords = [field(-1), field(H(parameter)), field(H(-parameter))]
        values = [coordinate-theta for coordinate in coords]
        assert all(values)
        assert [int(not value.is_square()) for value in values] == row
        rows.append(row)
    assert matrix(GF(2), rows).rank() == 3

    B, c = 27, 511
    n = (t**3-c*c)/(2*(B*c-1))
    d = n*(B*B*n+2)
    curve = EllipticCurve(function_field, [0, function_field(d)])
    first = curve(QQ(1)/9, (1+B*B*n)/B)
    second = curve(t, B*n+c)
    degrees = [x_degrees(point) for point in (first, second, first+second)]
    assert degrees == [(0, 0), (1, 0), (3, 2)]
    assert [2+max(d, n-2) for n, d in degrees] == [2, 2, 4]


if __name__ == "__main__":
    verify_quartic_mordell()
    verify_prime_mordell()
    verify_conductor_exponents()
    verify_generic_lattices()
    print("PM1–PM5 exact conductors, ranks, height lattices, and genuine prime verified")

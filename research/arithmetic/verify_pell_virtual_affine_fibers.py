"""Exact local, divisor-coordinate, and Frobenius controls for VF1–VF5."""

from sage.all import EllipticCurve, GF, QQ, ZZ, Integers, PolynomialRing, RealIntervalField, gcd, matrix, pari


from verify_quartic_pell_two_class_directions import P, t


def verify_laurent_pell_fibers():
    ring = PolynomialRing(QQ, names=("B", "u", "h"))
    B, u, h = ring.gens()
    fraction = ring.fraction_field()
    u = fraction(u)
    n = (u+1/u-2)/(2*B*B)
    y = (u-1/u)/2
    assert y*y == (1+B*B*n)**2-1
    assert 1+B*B*n+y == u
    assert (u+1/h) == u/h*(1/u+h)
    assert (u-h)*(u-1/h) == u*(2*B*B*n+2-h-1/h)
    c = (h+1)/B
    assert (u+h)*(1/u+h)/(B*B) == 2*(B*c-1)*n+c*c


def verify_split_fiber_norm_and_local_residues():
    for prime, q, root, omega, expected in [
        (7, 1, 3, 2, (1, 0)),
        (13, 0, 12, 3, (0, 2)),
    ]:
        field = GF(prime)
        P = field(q)**6
        root, omega = field(root), field(omega)
        assert root*root == P*P+1
        assert omega != 1 and omega**3 == 1
        unit = P+root
        row = 2*(unit-1)
        assert row*2*(P-root-1) == (-2*field(q)**2)**3
        for value, exponent in zip((unit, row), expected):
            assert value != 0 and value**((prime-1)//3) == omega**exponent
    assert matrix(GF(3), [[1, 0], [0, 2]]).determinant() != 0

    ring = PolynomialRing(Integers(9), "x")
    x = ring.gen()
    inert = ring.quotient(x*x-2, "r")
    r = inert.gen()
    assert (1+r)**8 == 1+3*r
    assert (2*r)**8 == 1


def trace_of_even_binomial(prime, exponent, degree):
    field = GF(prime**degree, "a")
    g = gcd(field.order()-1, exponent)
    generator = field.multiplicative_generator()**g
    power = field.one()
    character_sum = 0
    for _ in range((field.order()-1)//g):
        value = power+1
        if value:
            character_sum += 1 if value.is_square() else -1
        power *= generator
    return -2-g*character_sum


def verify_jacobian_orders():
    assert 5+1-trace_of_even_binomial(5, 4, 1) == 4
    traces = [trace_of_even_binomial(5, 12, degree) for degree in range(1, 6)]
    assert traces == [2, -14, -22, 58, 82]
    coefficients = [ZZ(1)]
    for degree in range(1, 6):
        numerator = -sum(coefficients[i]*traces[degree-i-1] for i in range(degree))
        assert numerator % degree == 0
        coefficients.append(numerator//degree)
    coefficients += [5**(degree-5)*coefficients[10-degree] for degree in range(6, 11)]
    assert sum(coefficients) == 2880
    assert ZZ(2880).valuation(3) == 2


def elliptic_quotient(polynomial):
    a, b, c = [polynomial[index] for index in (4, 2, 0)]
    assert polynomial == polynomial(-t)
    return EllipticCurve([0, b*b-4*a*c, 0, -4*a*a, 0])


def verify_fixed_quartic_quotient():
    a, b = P[4], P[2]
    curve = elliptic_quotient(P)
    minimal = EllipticCurve([1, 0, 0, -81381060, -282436946928])
    model_change = curve.isomorphism_to(minimal)
    assert model_change.tuple() == (QQ(8)/1153, QQ(327168)/1329409, QQ(4)/1153, 0)
    assert str(pari(minimal).ellrank()) == "[0, 0, 2, []]"
    for prime, count in [(7, 6), (11, 12)]:
        assert minimal.discriminant() % prime != 0
        assert minimal.change_ring(GF(prime)).cardinality() == count
    points = [
        (-5112, 2556), (-5256, -7812), (-5256, 13068),
        (15624, -1511172), (15624, 1495548),
    ]
    values = set()
    assert (~model_change)(minimal(points[0])) == curve(0, 0)
    for coordinates in points[1:]:
        X, Y = (~model_change)(minimal(coordinates)).xy()
        value = (Y/X-b)/(2*a)
        assert not value.is_square()
        values.add(value)
    assert values == {-QQ(4613)/442368, QQ(593923)/49152}

    conductor = minimal.conductor()
    assert conductor == 1128390 and conductor.is_squarefree()
    assert minimal.root_number() == 1
    enclosures = RealIntervalField(120)
    exponential = (-2*enclosures.pi()/enclosures(conductor).sqrt()).exp()
    cutoff = 1800
    coefficients = minimal.anlist(cutoff)
    partial = 2*sum(enclosures(coefficients[i])/i*exponential**i for i in range(1, cutoff+1))
    error = 4*exponential**(cutoff+1)/(1-exponential)
    assert (partial-error).lower() > 14
    assert (partial+error).upper() < 15


def verify_deformed_quartic_fibers():
    w = QQ(5)/24
    v = 1/(24*(1+2*w*w))
    z = w*v
    s = t-3/(8*v)
    Z = z*(s*s+3)
    H = (s+1)/2+v*(s*s+3)
    quartic = (
        s+3+6*v*(s+1)**2+12*v*v*(s+1)*(s*s+3)
        + 8*v**3*(s*s+3)**2-8*z*z*(s*s+3)
    )/(16*z)
    a, b = quartic[4], quartic[2]
    curve = elliptic_quotient(quartic)
    for t0, root, prime in [
        (QQ(269)/32, QQ(16409)/1280, 7),
        (-QQ(181)/96, QQ(181)/180, 11),
    ]:
        assert H(t0) == 0 and (quartic(t0)+Z(t0))**2 == root*root
        assert root*root == quartic(t0)**2+1
        field = GF(prime)
        assert field(curve.discriminant()) != 0
        reduction = curve.change_ring(field)
        assert reduction.cardinality() == 12
        triples = {3*point for point in reduction.points()}
        for sign in (-1, 1):
            unit = quartic(t0)+sign*root
            X = 2*a*unit
            Y = 2*a*(2*a*t0*t0+b)*unit
            assert reduction(field(X), field(Y)) not in triples


if __name__ == "__main__":
    verify_laurent_pell_fibers()
    verify_split_fiber_norm_and_local_residues()
    verify_jacobian_orders()
    verify_fixed_quartic_quotient()
    verify_deformed_quartic_fibers()
    print("VF1–VF5 exact polynomial and arithmetic controls passed")

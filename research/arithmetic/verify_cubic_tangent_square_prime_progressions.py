"""Exact elliptic and reflection checks for CTQ1–CTQ4."""
from sage.all import GF, QQ, ZZ, EllipticCurve, NumberField, PolynomialRing, pari

ring = PolynomialRing(QQ, "T")
T = ring.gen()
F = -3*T**4 + 6*T**2 + 1
assert F.is_squarefree()

rational = ring.fraction_field()
t = rational(T)
quadratic_polynomials = PolynomialRing(rational, "Z")
Z = quadratic_polynomials.gen()
quadratic = rational.extension(Z**2 - F, "Y")
y = quadratic.gen()
A, B, C = (t-1)**3, 4*t**3, (t+1)**3
alpha1 = (A+B-C+2*y)/4
alpha2 = (B+C-A+2*y)/4
assert alpha1 == alpha2*((y-1-t*t)/(2*t*(t+1)))**3
X = F/(4*t*t) - 1
V = y*(-3*t*t-1/(t*t))/(8*t)
assert V*V == X**3 + 1
assert V + 1 == alpha1*(-(y+3*t*t+2*t+1)/(4*t*t))**3
assert max(X.numerator().degree(), X.denominator().degree()) == 4

E1 = EllipticCurve([0, 0, 0, 0, 1])
assert E1.conductor() == 36
assert E1.quadratic_twist(5).conductor().valuation(2) == 2
assert E1.quadratic_twist(-3).conductor().valuation(3) == 2
assert (T**3 + 4).is_irreducible()

local = GF(19)
t0, y0 = local(5), local(9)
assert y0*y0 == F.change_ring(local)(t0)
alpha = ((t0-1)**3 + 4*t0**3 - (t0+1)**3 + 2*y0)/4
assert alpha == 6 and alpha**6 == 11


def fixture(p, h, expected_D, expected_polynomial):
    p, h = ZZ(p), ZZ(h)
    primes = (p-h, p, p+h)
    assert all(q.is_prime() for q in primes)
    assert p % 3 == 1 and h % 18 == 0
    assert (p % 19, h % 19) == (5, 1)
    a, b, c = (p-h)**3, 4*p**3, (p+h)**3
    radicand = 3*p**4-6*p*p*h*h-h**4
    common = h*h*radicand
    assert 4*common == 4*a*c - (a+c-b)**2
    assert all(common % q for q in primes)
    assert radicand.valuation(3) == 1
    assert (-radicand // 3) % 3 == 2 and (-radicand) % 8 == 5
    D = (-radicand).squarefree_part()
    assert D == expected_D and ZZ(pari(D).quaddisc()) == D
    u1, u2 = (a+b-c)//2, (b+c-a)//2
    assert (u1-2*p**3).valuation(3) >= 3 and (u2-2*p**3).valuation(3) >= 3
    K = NumberField(T*T+common, "radical")
    radical = K.gen()
    first, second = (u1+radical)/2, (u2+radical)/2
    third = (b-a-c)/2 + radical
    assert first.norm() == ((p-h)*p)**3 and second.norm() == (p*(p+h))**3
    assert third.norm() == ((p-h)*(p+h))**3 and first*second == p**3*third
    polynomial = T**3 - 3*(p-h)*p*T - u1
    assert polynomial == expected_polynomial and polynomial.is_irreducible()
    assert ZZ(pari(polynomial).nfdisc()) == -D//3
    cubic = pari(polynomial).nfinit()
    assert sorted((int(q[2]), int(q[3])) for q in cubic.idealprimedec(3)) == [(1, 1)]*3
    for q in primes:
        inert = q != p and GF(q)(2)**((q-1)//3) != 1
        expected = [(1, 3)] if inert else [(1, 1)]*3
        assert sorted((int(v[2]), int(v[3])) for v in cubic.idealprimedec(q)) == expected
    assert EllipticCurve([0, 0, 0, 0, D**3]).conductor() == 4*D**2


fixture(2551, 1008, -86341391658723, T**3-11808579*T-12498602366)
fixture(2437, 2376, 127222490176557, T**3-445971*T+26799849502)

print("Verified cubic norm relations, Mordell map, and signed reflected prime-path fixtures.")

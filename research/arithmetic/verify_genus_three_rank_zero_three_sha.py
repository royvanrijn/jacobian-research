"""Exact HG13–HG14 certificates, SageMath 10.9.

Checks the gluing inputs and an explicit rank-zero Sha lower-bound fixture.
The unrestricted growth and density conclusions use the written normal-order
argument, Cassels' formula and Smith's imported theorem; this script does not
reprove those theorems or compute complete Tate–Shafarevich groups.
"""
from sage.all import (QQ, ZZ, PolynomialRing, EllipticCurve,
                      HyperellipticCurve, GF, gcd, prod, crt)

parameters = PolynomialRing(QQ, "s")
s = parameters.gen()
K = parameters.fraction_field()
h = -3*(s*s+31)/(2*(s*s+49))
b = 2*h+3
a = h*h+12*h+9
B = b**3
Delta = a*a-4*B
u = -(7*h+5)/2
v = (h+1)*s/2
assert b == 54/(s*s+49)
assert h+1 == (5-s*s)/(2*(s*s+49))
assert h-3 == -9*(s*s+43)/(2*(s*s+49))
assert a == -27*(s**4+62*s*s+529)/(4*(s*s+49)**2)
assert Delta == (h+1)*(h-3)**3
assert u*u+v*v == -a/(2*b)

R = PolynomialRing(K, "t")
t = R.gen()
X = u*(1-t*t)-2*v*t
Y = v*(1-t*t)+2*u*t
Z = b*(1+t*t)
qa = -a*a*B/2
qc = a*a*(a*a+4*B)/(16*B)
lam = a*a/(4*B)
delta3 = (1-lam)**2
product_root = Delta*(1-lam)
assert 2*B*(X*X+Y*Y)+a*Z*Z == 0
assert product_root**2 == Delta**2*delta3
assert product_root*(2*a*a/Delta+(1+lam)**2/delta3-1)+2*a*a*(1+lam) == 0
assert -qa*(2/a)*(1/B) == a
assert -qc*(2/a)**2 == -(1+lam)
octic = b*(16*b*b*(X*Y)**2-Delta*(1+t*t)**4)
assert qa*(X**4+Y**4)+qc*Z**4 == a*a/16*octic
assert octic.degree() == 8

# Direct normalized degree-three isogeny E_h -> E'_h, with the target
# translated to its rational two-torsion point. This fixes its direction.
RX = PolynomialRing(K, "x").fraction_field()
x = RX.gen()
source_cubic = x**3+(h*x+h+1)**2
U = x+4*h*(h+1)/x+4*(h+1)**2/x**2-(4*h+3)
y_multiplier = 1-4*h*(h+1)/x**2-8*(h+1)**2/x**3
assert y_multiplier*y_multiplier*source_cubic == U*(U*U+a*U+B)
assert U.numerator().degree() == 3

# The extra prime used to control the companion's root number is always
# available for even s. Its disjointness from all bad E'_h primes is exact.
P = s**4+62*s*s+529
for q, residual in [(s*s-5, 864), (s*s+49, -108), (s*s+43, -288)]:
    assert P % q == residual
for j in range(6):
    if j % 2 == 0:
        assert ZZ(P(j)) % 2 == 1 and ZZ(P(j)) % 3 != 0
print("Generic gluing, rational conic, normalized octic and three-isogeny: PASS")


def curves(ss, d=1):
    hh, aa, bb = h(ss), a(ss), b(ss)
    E = EllipticCurve([0, d*hh*hh, 0,
                      d*d*2*hh*(hh+1), d**3*(hh+1)**2])
    Ep = EllipticCurve([0, d*aa, 0, d*d*bb**3, 0])
    ll = aa*aa/(4*bb**3)
    F = EllipticCurve([0, -1-ll, 0, ll, 0]).quadratic_twist(d)
    return E, Ep, F


def local_ratio(ss, d=1):
    E, Ep, _ = curves(ss, d)
    xx = PolynomialRing(QQ, "x").gen()
    phi = E.isogeny(xx)
    quotient = phi.codomain()
    assert quotient.is_isomorphic(Ep)
    em, epm = E.minimal_model(), quotient.minimal_model()
    psi = E.isomorphism_to(em)*phi.dual()*(~quotient.isomorphism_to(epm))
    sf = psi.scaling_factor()
    primes = (set(ZZ(em.discriminant()).prime_divisors())
              | set(ZZ(epm.discriminant()).prime_divisors()) | {2, 3})
    ratios = {}
    for p in sorted(primes):
        cp = QQ(em.tamagawa_number(p))/epm.tamagawa_number(p)
        if p == 3:
            cp *= QQ(3)**sf.valuation(3)
        ratios[p] = cp
    arch = QQ(1) if d > 0 else QQ(1)/3
    return arch*prod(ratios.values()), ratios


for ss in [0, 1, 2, 3, 4, 5, 8, 11, 12, 20]:
    c, locals_ = local_ratio(ss)
    for p, cp in locals_.items():
        if p < 5:
            assert cp >= QQ(1)/3
            continue
        expected = QQ(1)
        if (ss*ss-5) % p == 0 or (ss*ss+49) % p == 0:
            expected *= 3
        if (ss*ss+43) % p == 0 and GF(p)(-3).is_square():
            expected /= 3
        assert cp == expected, (ss, p, cp, expected)
    W = (sum(p >= 5 for p in ZZ(abs(ss*ss-5)).prime_divisors())
         + sum(p >= 5 for p in ZZ(ss*ss+49).prime_divisors())
         - sum(p >= 5 for p in ZZ(ss*ss+43).prime_divisors()))
    assert c >= QQ(3)**(W-2)
print("Finite local controls agree with the written all-prime table and bound")

# Validate a nontrivial companion root-number control at an extra prime.
for ss in [4, 6, 8]:
    E, Ep, F = curves(ss)
    extra = ZZ(P(ss)).prime_divisors()[0]
    assert extra >= 5 and E.has_good_reduction(extra)
    assert Ep.has_good_reduction(extra) and F.conductor().valuation(extra) == 1
    desired_symbol = Ep.root_number()*F.root_number()
    residue = next(j for j in range(1, extra) if GF(extra)(j).is_square() == (desired_symbol == 1))
    assert ZZ(residue).kronecker(extra) == desired_symbol
    primes = sorted(set(ZZ(E.conductor()*F.conductor()).prime_divisors()) | {2, 3})
    moduli = [ZZ(8)]+[p for p in primes if p != 2]
    residues = [ZZ(1)]+[residue if p == extra else ZZ(1) for p in primes if p != 2]
    representative = crt(residues, moduli)
    assert gcd(representative, prod(moduli)) == 1
print("Extra companion prime and simultaneous root-number residue controls: PASS")

# Concrete fixture independent of Smith's density theorem: two exact elliptic
# rank-zero descents, global Selmer ratio nine, and the prime-to-three gluing.
ss, d = 8, 57
E, Ep, F = curves(ss, d)
assert E.torsion_order() % 3 != 0
assert Ep.torsion_order() % 3 != 0
Ep, F = Ep.minimal_model(), F.minimal_model()
assert Ep.rank(proof=True, use_database=False) == 0
assert F.rank(proof=True, use_database=False) == 0
assert Ep.rank_bound(algorithm="mwrank") == F.rank_bound(algorithm="mwrank") == 0
c, _ = local_ratio(ss, d)
assert c == 9
assert list(Ep.minimal_model().a_invariants()) == [1, 1, 0, -546524405, -4917455709459]
assert list(F.minimal_model().a_invariants()) == [0, -1, 0,
    -32816655601327594443617, -2284648019026739769616608032943135]

RT = PolynomialRing(QQ, "t")
ff = RT([coefficient(ss) for coefficient in octic])
assert ff.is_squarefree()
checked = 0
for p in [5, 17, 31]:
    if any(coefficient.denominator() % p == 0 for coefficient in ff):
        continue
    fp = ff.change_ring(GF(p))*GF(p)(d)
    if not fp.is_squarefree():
        continue
    C = HyperellipticCurve(fp)
    expected = Ep.change_ring(GF(p)).cardinality()**2*F.change_ring(GF(p)).cardinality()
    assert C.frobenius_polynomial()(1) == expected
    checked += 1
assert checked >= 2
print("s=8,d=57: exact elliptic rank bounds (0,0), c(psi_d)=9;")
print("genus-three rank zero and Sha[3] contains (Z/3)^4 by the written proof")
print("HG13–HG14 certificates PASS; no complete Sha computation or independent review")

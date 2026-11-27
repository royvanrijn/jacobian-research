"""Exact quotient, Kummer, conductor and isogeny checks for PQ1–PQ4."""
from sage.all import GF, QQ, ZZ, EllipticCurve, Matrix, PolynomialRing

ring = PolynomialRing(QQ, "T")
T = ring.gen()
F = -3*T**6 + 24*T**4 + 24*T**2 + 4
assert F.is_squarefree()
rational = ring.fraction_field()
t = rational(T)
quadratic_ring = PolynomialRing(rational, "Z")
Z = quadratic_ring.gen()
quadratic = rational.extension(Z**2-F, "Y")
y = quadratic.gen()
alpha1 = (t**3-6*t*t-2+y)/2
alpha2 = (t**3+6*t*t+2+y)/2
gamma = (t**3-6*t*t-2-y)*alpha2/2
eta = alpha1*alpha2
X1, Y1 = -3*t*t-4, -3*y
X2, Y2 = 4*(1/(t*t)-1), 4*y/t**3
for X, Y, b in [(X1, Y1, 14), (X2, Y2, 28)]:
    assert Y**2 == X**3 + (6*X+b)**2
    assert max(X.numerator().degree(), X.denominator().degree()) == 2
assert (Y1+6*X1+14)*gamma == (y+5*t*t+2)**3
assert (Y2+6*X2+28)*eta == -8*(t*t-1)**3


def model(b, D):
    return EllipticCurve([0, 36*D, 0, 12*D*b, D*b*b])


base = [model(b, 1) for b in [14, 28]]
assert [E.conductor() for E in base] == [756, 189]
assert [E.change_ring(GF(13)).trace_of_frobenius() for E in base] == [2, -4]
assert (T**3+48*T*T+336*T+784).is_irreducible()
polynomials = [ring(E.division_polynomial(3).list()) for E in base]
assert polynomials == [3*T*(T**3+48*T*T+336*T+784), 3*T*(T+28)*(T*T+20*T+112)]

for b0, tamagawa, multipliers, ratios, conductor_exponents in [
    (14, [(1, 3), (1, 3), (3, 1)], [0, 1, 0], [3, 9, QQ(1)/3], [2, 3, 1]),
    (28, [(1, 1), (3, 1), (3, 1)], [0, 1, 0], [1, 1, QQ(1)/3], [0, 3, 1]),
]:
    for index, (p, d) in enumerate([(2, 5), (3, -3), (7, 1)]):
        E = model(b0*d, d)
        companion = model((27*b0*d-4*d*6**3)//9, -3*d)
        minimal, target = E.global_minimal_model(), companion.global_minimal_model()
        cE = minimal.local_data(p).tamagawa_number()
        cH = target.local_data(p).tamagawa_number()
        source_scale = E.isomorphism_to(minimal).scaling_factor()
        target_scale = companion.isomorphism_to(target).scaling_factor()
        valuation = (target_scale/source_scale).valuation(p)
        assert (cE, cH) == tamagawa[index] and valuation == multipliers[index]
        assert QQ(cH)/cE*QQ(p)**valuation == ratios[index]
        assert minimal.local_data(p).conductor_valuation() == conductor_exponents[index]



for p, d, expected_ratio in [(2, 5, 1), (3, -3, QQ(1)/3), (7, 1, QQ(1)/3)]:
    recentered = EllipticCurve([0, -48*d, 0, 672*d*d, -2352*d**3])
    expected_target = EllipticCurve([0, 144*d, 0, 96*d*d, 16*d**3])
    assert recentered.is_isomorphic(model(28*d, d))
    isogeny = recentered.isogeny(ring(T), codomain=expected_target)
    assert isogeny.degree() == 3
    assert isogeny.scaling_factor().valuation(p) == 0
    raw_target = expected_target
    source, target = recentered.global_minimal_model(), raw_target.global_minimal_model()
    source_scale = recentered.isomorphism_to(source).scaling_factor()
    target_scale = raw_target.isomorphism_to(target).scaling_factor()
    assert (target_scale/source_scale).valuation(p) == 0
    assert QQ(target.local_data(p).tamagawa_number())/source.local_data(p).tamagawa_number() == expected_ratio


def fixture(p, h, expected_D, residue_rows):
    p, h = ZZ(p), ZZ(h)
    assert h % 126 == 0 and p % 7 != 0
    assert all(q.is_prime() for q in [p-h, p, p+h])
    radicand = -3*p**6+24*p**4*h*h+24*p*p*h**4+4*h**6
    D = radicand.squarefree_part()
    assert D == expected_D and D % 8 == 5 and D.valuation(3) == 1
    assert (D//3) % 3 == 2 and GF(7)(D).is_square()
    A, B, C = (p-h)**3, p**3, (p+h)**3
    traces = (A+B-C, B+C-A)
    rows = []
    for q, radical, zeta, exponents in residue_rows:
        field = GF(q)
        r, zeta = field(radical), field(zeta)
        assert r*r == field(radicand)
        assert zeta**3 == 1 and zeta != 1
        values = [(field(trace)+r)/2 for trace in traces]
        assert tuple(value**((q-1)//3) for value in values) == tuple(zeta**e for e in exponents)
        rows.append(exponents)
    assert Matrix(GF(3), rows).rank() == 2
    assert [model(b*D, D).conductor() for b in [14, 28]] == [84*D**2, 21*D**2]


fixture(383, 126, -367119482556963, [(109, 46, 63, (0, 2)), (199, 99, 106, (2, 0))])
fixture(137, 126, 675712099077, [(37, 17, 26, (2, 0)), (73, 22, 8, (2, 1))])

print("Verified the two elliptic projections, marked cube identities, and local isogeny and conductor data.")

"""Exact HG11–HG12 certificates; SageMath 10.9, no random input.

This is an intermediate two-direction family. It does not certify a third
rational three-torsion class, full torsion [6,6], or independent specialized
ideal classes. Run with sage -python from the repository root.
"""
from sage.all import (QQ, ZZ, PolynomialRing, HyperellipticCurve, GF,
                      QuadraticField, binomial, cyclotomic_polynomial,
                      euler_phi, gcd, prod)

A = PolynomialRing(QQ, "a")
a = A.gen()
R = PolynomialRing(A, "t")
t = R.gen()
q1 = t**3+a*t**2-(4*a+30)*t-1
q2 = t**3+(4*a+27)*t**2+(5*a+27)*t-4
g1 = ((QQ(7)/2*a+QQ(225)/7)*t**4
      +(QQ(241)/14*a+QQ(803)/7)*t**3
      +(QQ(9)/2*a*a+QQ(127)/2*a+QQ(1620)/7)*t**2
      +(6*a+45)*t+1)
g2 = ((-QQ(7)/2*a-QQ(216)/7)*t**4
      +(-QQ(151)/14*a-QQ(569)/7)*t**3
      +(QQ(9)/2*a*a+QQ(71)/2*a+QQ(297)/7)*t**2
      +(-15*a-81)*t+8)
H = QQ(7)/3*(t*t+t+1)
assert g1-g2 == H*(q2-q1)
assert H*(g1+g2) == q1*q1+q1*q2+q2*q2
raw = g1*g1+q1**3
assert raw == g2*g2+q2**3
f, remainder = raw.quo_rem(t*t)
assert not remainder and f.degree() == 7 and f.is_monic()
assert g1[0] == 1 and g2[0] == 8
assert q1[0] == -1 and q2[0] == -4

T = PolynomialRing(QQ, "t")
tt = T.gen()


def specialize(poly, value):
    return T([c(value) for c in poly])


ff = specialize(f, -9)
qq1, qq2, gg1, gg2 = [specialize(h, -9) for h in (q1, q2, g1, g2)]
integral_f = 196*ff
assert integral_f == (196*tt**7-5211*tt**6+41022*tt**5+116149*tt**4
                      -69840*tt**3+13896*tt**2+4480*tt-1008)
assert ff.is_squarefree() and ff.is_irreducible()
assert gcd(gg1, qq1*tt) == gcd(gg2, qq2*tt) == 1

frobenius = {}
for p, expected_count in [(5, 108), (37, 46215)]:
    k = GF(p)
    fp = ff.change_ring(k)
    assert fp.is_squarefree()
    C = HyperellipticCurve(fp)
    F = C.frobenius_polynomial()
    assert F(1) == expected_count
    frobenius[p] = F
    J = C.jacobian()(k)
    marked = []
    for q, g in [(qq1, gg1), (qq2, gg2)]:
        u = q.change_ring(k)
        assert gcd(g.change_ring(k), u*u.parent().gen()) == 1
        vv = (-g.change_ring(k)*u.parent().gen().inverse_mod(u)) % u
        assert (vv*vv-fp) % u == 0
        D = J([u, vv])
        assert 3*D == J(0)
        marked.append(D)
    span = []
    for i in range(3):
        for j in range(3):
            D = i*marked[0]+j*marked[1]
            assert D not in span
            span.append(D)
    print(f"p={p}: good reduction, marked span 9, #J={F(1)}")
assert gcd([F(1) for F in frobenius.values()]) == 9

# A separate elementary enumeration cross-checks the small-prime Frobenius.
counts = []
for n in (1, 2, 3):
    k = GF(5**n, name="b")
    fp = ff.change_ring(k)
    counts.append(1+sum(1 if not fp(x) else (2 if fp(x).is_square() else 0)
                        for x in k))
power_sums = [None]+[5**n+1-counts[n-1] for n in (1, 2, 3)]
cc = [QQ(1)]
for n in (1, 2, 3):
    cc.append(-sum(cc[n-i]*power_sums[i] for i in range(1, n+1))/n)
X = frobenius[5].parent().gen()
assert frobenius[5] == (X**6+cc[1]*X**5+cc[2]*X**4+cc[3]*X**3
                        +5*cc[2]*X**2+25*cc[1]*X+125)
assert frobenius[5] == X**6-X**5+8*X**3-25*X+125

# Exclude every possible root-of-unity quotient of distinct Frobenius roots.
# Such a quotient has degree at most 36, hence order at most 2*36^2:
# phi(m)^2 >= m/2 follows directly from the prime factorization formula.
assert frobenius[5].is_irreducible()
Z = PolynomialRing(QQ, "z")
z = Z.gen()
RX = PolynomialRing(Z, "x")
x = RX.gen()
P = RX(frobenius[5].list())
ratio_resultant = Z(P.resultant(P(x*z)))
residual, remainder = ratio_resultant.quo_rem((z-1)**6)
assert not remainder and residual(1) != 0
orders_checked = []
for m in range(2, 2*36**2+1):
    if euler_phi(m) <= 36:
        assert residual % Z(cyclotomic_polynomial(m)) != 0
        orders_checked.append(m)
print(f"Absolute simplicity: irreducible Frobenius and {len(orders_checked)} "
      "possible nontrivial root-of-unity orders excluded")

# Classical transvectants give a geometric invariant of the binary octic.
B = PolynomialRing(A, names=("X", "Z"))
XX, ZZvar = B.gens()
octic = sum(c*XX**i*ZZvar**(8-i) for i, c in enumerate(f))


def derivative(poly, variable, order):
    for _ in range(order):
        poly = poly.derivative(variable)
    return poly


def transvectant(left, right, order):
    return sum((-1)**i*binomial(order, i)
               *derivative(derivative(left, XX, order-i), ZZvar, i)
               *derivative(derivative(right, XX, i), ZZvar, order-i)
               for i in range(order+1))


J2 = A(transvectant(octic, octic, 8))
quartic_degree_covariant = transvectant(octic, octic, 4)
J4 = A(transvectant(quartic_degree_covariant, quartic_degree_covariant, 8))
invariant = J4/J2**2
assert invariant(-9) == QQ(448591166316421718490)/3947463686715586303
assert invariant(-8) == QQ(182782757428770535)/1009280839235802
assert invariant(-9) != invariant(-8)
print("Nonconstant geometric invariant J4/J2^2: PASS")

# Integral functions retain the same cubic Kummer classes by a rational cube
# scaling. Their common value at t=0 is a cube, enabling simultaneous local
# cube conditions at every exceptional prime.
G1, G2 = 14*gg1, 14*gg2
assert all(c in ZZ for poly in (integral_f, G1, G2, qq1, qq2) for c in poly)
assert G1.resultant(qq1) == -37212507
assert G2.resultant(qq2) == -2596485888
bad = {2, 3, 7, 17, 19, 251, 1163}
assert set(ZZ(37212507*2596485888).prime_divisors()) | {2, 3, 7} == bad
modulus = ZZ(2)**7*3**3*7**4*17*19*251*1163
for p in bad:
    assert modulus.valuation(p) > 3*max(ZZ(14).valuation(p), ZZ(28).valuation(p)) + 2*ZZ(3).valuation(p)
for G, q, cube_root in [(G1, qq1, 14), (G2, qq2, 28)]:
    assert (196*G)**2-(196*tt)**2*integral_f == -(196*q)**3
    assert 196*G(0) == cube_root**3
    # Both integer coefficient differences are divisible by t.
    assert (196*G-cube_root**3) % tt == 0

# Negative control: the congruence restriction matters. At t=1 the second
# function has valuations 4 and 5 at the primes over 2, so is not virtual.
assert integral_f(1) == 36*2769
K = QuadraticField(2769, "s")
bad_value = K(196)*(G2(1)+6*K.gen())
two_valuations = sorted(bad_value.valuation(P) for P in K.primes_above(2))
assert two_valuations == [4, 5]
print("All-prime virtuality anchors, resultants and Hensel bounds: PASS")
print("Negative control t=1: valuations [4,5] above 2, as required")
print("HG11–HG12 PASS: infinite geometrically distinct simple [3,3] family "
      "with rational Weierstrass point; two simultaneous virtual directions "
      "on the stated congruence. No third torsion direction or [6,6] claimed.")

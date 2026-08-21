#!/usr/bin/env sage
"""Verify the discriminant-43 boundary point in the Humbert-8 Kumar chart.

This combines three independent exact fingerprints:

* the norm-43 Gross vector and its complement in the recovered generic
  transcendental lattice;
* a rational point (r,s) of the Humbert-8 plane where two height-5/2
  section loci meet, together with its oriented quadratic coordinate;
* point counts at good primes, which agree with the weight-three CM form of
  discriminant -43.

The finite Frobenius list is a reproducible identification fingerprint, not
by itself a proof of equality of two characteristic-zero K3 models.
"""

from sage.all import *


# The recovered rank-three transcendental lattice and the norm-43 Gross
# vector.  In the exact Gross basis
#
#   (i+j+37*k, 2*j+112*k, 158*k)
#
# the vector has coordinates (-11,5,-1), hence beta=-11*i-j-5*k.
T_generic = matrix(ZZ, (
    (-10, -2, -1),
    (-2, 10, 3),
    (-1, 3, 10),
))
v43 = vector(ZZ, (169, 167, -128))
assert T_generic.det() == -948
assert v43 * T_generic * v43 == -40764
assert gcd([abs(ZZ(value)) for value in T_generic * v43]) == 948
orthogonal = matrix(ZZ, 1, 3, list(v43 * T_generic)).right_kernel_matrix()
T43 = orthogonal * T_generic * orthogonal.transpose()
assert T43 == matrix(ZZ, ((22, 1), (1, 2)))
assert T43.det() == 43

quaternion = QuaternionAlgebra(6)
ii, jj, kk = quaternion.gens()
gross_basis = (ii + jj + 37*kk, 2*jj + 112*kk, 158*kk)
beta43 = -11*gross_basis[0] + 5*gross_basis[1] - gross_basis[2]
assert beta43 == -11*ii - jj - 5*kk
assert beta43.reduced_norm() == 43


# Humbert-8 coordinates.  A height-5/2 section through the nonidentity E7
# component has
#
#   x = T^2*(c0+c*T+W*T^2),
#   y = T^3*(d0+d1*T+d2*T^2+d3*T^3).
#
# Eliminating its coefficients has two nondegenerate components G and F.
# Their rational intersections include the CM-43 point below.
RS = PolynomialRing(QQ, names=("r", "s"))
r, s = RS.gens()
G = -4*r**2 + 16*r*s - 4*r - 1
F = (
    1024*r**4 + 5184*r**3*s + 3584*r**3 + 6561*r**2*s**2
    + 7776*r**2*s + 3840*r**2 - 2592*r*s + 1664*r + 256
)

Rr = PolynomialRing(QQ, "rho")
rho = Rr.gen()
s_on_G = (2*rho + 1)**2 / (16*rho)
intersection_polynomial = Rr((16*rho)**2 * F(r=rho, s=s_on_G))
assert intersection_polynomial == (
    (242*rho + 25) * (722*rho + 1225) * rho**2 * (2*rho + 1)**2
)

r43 = -QQ(1225) / 722
s43 = -QQ(93312) / 442225
assert s43 == s_on_G(rho=r43)
assert G(r=r43, s=s43) == F(r=r43, s=s43) == 0


# The exact Humbert-8 E7+E8 equation.
RT = PolynomialRing(QQ, "T")
T = RT.gen()
A1 = 2*r43*s43**2
A = -(9*r43*s43 + 4*r43**2 + 4*r43 + 1) / 3
B1 = r43*s43**2 * (3*s43 + 8*r43 - 2) / 3
B = -(
    54*r43**2*s43 + 81*r43*s43 - 16*r43**3
    - 24*r43**2 - 12*r43 - 2
) / 27
B2 = r43**2
a4 = RT(A1*T**3 + A*T**4)
a6 = RT(B1*T**5 + B*T**6 + B2*T**7)
Delta = RT(-16 * (4*a4**3 + 27*a6**2))
assert Delta.valuation(T) == 9
assert (8-a4.degree(), 12-a6.degree(), 24-Delta.degree()) == (4, 5, 10)

# The oriented Hilbert modular coordinate has exactly the CM field Q(sqrt(-43)).
branch43 = (
    16*r43*s43**2 + 32*r43**2*s43 - 40*r43*s43 - s43
    + 16*r43**3 + 24*r43**2 + 12*r43 + 2
)
z_scale = QQ(11664) / 6859
assert 2*branch43 == -43*z_scale**2

# In the rational (m,n) chart on the oriented Humbert surface, m is rational
# and the entire CM orientation field is carried by n.  The sign of n follows
# the independent choice of sign for z.
s1_43 = s43 + r43 - QQ(5)/4 - 1/(32*r43)
m43 = 32*r43*s1_43/(16*r43 - 1)
n_scale43 = z_scale/(16*r43 - 1)
assert m43 == -QQ(2468019)/407569
assert n_scale43 == -QQ(1296)/21451
assert r43 == (m43**2 - 1) / (16*(2*(-43*n_scale43**2) - 1))
assert s43 == (
    m43*(16*r43 - 1)/(32*r43) - r43 + QQ(5)/4 + 1/(32*r43)
)

# Kumar's inverse map to the Clebsch--Igusa invariants.
I2 = -4 * (8*r43 + 3*s43 - 2)
I4 = 4 * (4*r43**2 + 9*r43*s43 + 4*r43 + 1)
I6 = -4 * (
    48*r43**3 + 94*r43**2*s43 + 36*r43*s43**2
    + 40*r43**2 - 35*r43*s43 + 4*r43 + 4*s43 - 2
)
I10 = -8*r43**3*s43**2
assert (I2, I4, I6, I10) == (
    QQ(28667544)/442225,
    QQ(4665600)/130321,
    QQ(1855940221440)/2305248169,
    QQ(10666233446400)/6131066257801,
)


# Two exact height-5/2 sections at the intersection G cap F.
c0 = QQ(1194481) / 442225
W = QQ(663613890625) / 34828517376
assert W.is_square() and W.sqrt() == QQ(814625) / 186624
section_data = (
    (
        QQ(684775)/93312,
        (
            QQ(914233879)/294079625,
            QQ(26371835)/1181952,
            QQ(557834834375)/11609505792,
            QQ(540596465650390625)/6499837226778624,
        ),
    ),
    (
        -QQ(1765225)/93312,
        (
            -QQ(1085766121)/294079625,
            QQ(57241835)/1181952,
            -QQ(1437996415625)/11609505792,
            QQ(540596465650390625)/6499837226778624,
        ),
    ),
)
sections = []
for c, d in section_data:
    X = RT(T**2 * (c0 + c*T + W*T**2))
    Y = RT(T**3 * sum(d[index]*T**index for index in range(4)))
    assert Y**2 == X**3 + a4*X + a6
    assert X.valuation(T) == 2 and Y.valuation(T) == 3
    sections.append((X, Y))


def additive_path(X, Y):
    """Follow a section until it reaches a smooth chart over T=0."""
    chart = PolynomialRing(QQ, names=("z", "xx", "yy"))
    z, xx, yy = chart.gens()

    def embed(poly):
        return sum(chart(poly[index])*z**index for index in range(poly.degree()+1))

    surface = yy**2 - xx**3 - embed(a4)*xx - embed(a6)
    section_x, section_y = RT(X), RT(Y)
    center_x = center_y = QQ(0)
    path = []
    for _ in range(10):
        section_x = RT((section_x-center_x)//T)
        section_y = RT((section_y-center_y)//T)
        surface = chart(surface(z, center_x+z*xx, center_y+z*yy)//z**2)
        center_x, center_y = section_x(0), section_y(0)
        center = {z: 0, xx: center_x, yy: center_y}
        gradient = tuple(
            surface.derivative(variable).subs(center) for variable in (z, xx, yy)
        )
        path.append((center_x, center_y, gradient))
        if any(gradient):
            return tuple(path)
    raise RuntimeError("E7 path did not reach a smooth chart")


paths = tuple(additive_path(X, Y) for X, Y in sections)
assert tuple(len(path) for path in paths) == (3, 3)
assert paths[0][:2] == paths[1][:2]
assert paths[0][2][:2] != paths[1][2][:2]

# There are no finite smooth intersections.  At infinity the two sections
# meet transversally once on the smooth identity component.
assert gcd(sections[0][0]-sections[1][0], sections[0][1]-sections[1][1]) == T**3
Sinfinity = PolynomialRing(QQ, "q")
q = Sinfinity.gen()


def reverse(poly, weight):
    return Sinfinity(sum(poly[index]*q**(weight-index) for index in range(poly.degree()+1)))


dx_infinity = reverse(sections[0][0]-sections[1][0], 4)
dy_infinity = reverse(sections[0][1]-sections[1][1], 6)
assert dx_infinity.valuation(q) == dy_infinity.valuation(q) == 1

# E7 contributes 3/2 to each nonidentity section and to their mixed pairing;
# E8 is unimodular.  Both sections are polynomial and hence disjoint from O.
height_self = QQ(4) - QQ(3)/2
height_mixed = QQ(2) - 1 - QQ(3)/2
assert height_self == QQ(5)/2 and height_mixed == -QQ(1)/2


# The skipped D9+E7 two-neighbor recovers the generic height-four Humbert
# section from the second point (T,V)=(0,-s) on its pointed quartic.  In the
# short Kumar model its coordinates are the following compact expressions.
quartic_b = -3*s43*T**2
quartic_c = 3*s43**2*T + (2*r43+1)*T**2
quartic_d = -s43**3 - 2*r43*s43*T + r43*T**2 - 2*s43*T
height4_x_long = quartic_d**2/(4*s43**2) - quartic_c
height4_x = RT(height4_x_long + quartic_c/3)
height4_y = RT(-quartic_b*s43 - quartic_d*height4_x_long/(2*s43))
assert height4_y**2 == height4_x**3 + a4*height4_x + a6
assert (height4_x.degree(), height4_y.degree()) == (4, 6)
assert (height4_x(0), height4_y(0)) == (s43**4/4, s43**6/8)

# In the CM MW basis (P1,P2,P3), where P1,P2 are the two height-5/2
# sections above and P3 is the height-four section, the generic H2 lattice is
# generated by P3 and 4*P1-5*P2+P3.  The latter is the level-79 direction.
height43 = matrix(QQ, (
    (QQ(5)/2, -QQ(1)/2, -1),
    (-QQ(1)/2, QQ(5)/2, 0),
    (-1, 0, 4),
))
height4_coordinates = vector(ZZ, (0, 0, 1))
level79_coordinates = vector(ZZ, (4, -5, 1))
assert height4_coordinates*height43*height4_coordinates == 4
assert level79_coordinates*height43*level79_coordinates == QQ(237)/2
assert height4_coordinates*height43*level79_coordinates == 0

# The exact group law exhibits why a direct H237 equation is so large.  At
# the CM point the level-79 section has a common denominator h with deg(h)=58;
# x and y have denominators h^2 and h^3 and numerator degrees 120 and 180.
KT = RT.fraction_field()
curve43 = EllipticCurve(KT, [0, 0, 0, a4, a6])
point1, point2 = tuple(curve43(point) for point in sections)
point3 = curve43((height4_x, height4_y))
level79_point = 4*point1 - 5*point2 + point3
level79_denominator = level79_point[0].denominator().sqrt()
assert level79_denominator.degree() == 58
assert level79_point[0].denominator() == level79_denominator**2
assert level79_point[1].denominator() == level79_denominator**3
assert (
    level79_point[0].numerator().degree(),
    level79_point[1].numerator().degree(),
) == (120, 180)


# A short exact Frobenius fingerprint.  The raw Weierstrass fibers at zero
# and infinity acquire 7 and 8 exceptional components respectively, hence
# the 15*p correction below.
residual = RT(Delta / T**9)
residual_disc = residual.discriminant()
bad_primes = set(prime_divisors(abs(residual_disc.numerator())))
bad_primes.update(prime_divisors(residual_disc.denominator()))
assert bad_primes == {2, 3, 5, 7, 19, 29, 37, 43, 241, 349}


def cm43_coefficient(p):
    if kronecker(-43, p) == -1:
        return ZZ(0)
    for b in range(0, floor(sqrt(QQ(4*p)/43))+1):
        a2_value = 4*p - 43*b*b
        if a2_value < 0 or not ZZ(a2_value).is_square():
            continue
        aa = ZZ(a2_value).sqrt()
        if (aa-b) % 2 == 0:
            return ZZ((aa*aa - 43*b*b) // 2)
    raise RuntimeError(f"failed to represent split prime {p}")


def k3_point_count(p):
    assert p not in bad_primes
    field = GF(p)

    def reduce(value):
        value = QQ(value)
        return field(value.numerator()) / field(value.denominator())

    a1p, ap, b1p, bp, b2p = map(reduce, (A1, A, B1, B, B2))
    total = ZZ(0)
    for t_value in field:
        a4_value = a1p*t_value**3 + ap*t_value**4
        a6_value = b1p*t_value**5 + bp*t_value**6 + b2p*t_value**7
        character_sum = ZZ(0)
        for x_value in field:
            rhs = x_value**3 + a4_value*x_value + a6_value
            character_sum += 0 if rhs == 0 else (1 if rhs.is_square() else -1)
        total += p + 1 + character_sum
    total += p + 1       # raw cuspidal cubic at infinity
    total += 15*p        # resolve III* and II*
    return total


frobenius = []
for prime in (11, 13, 17, 23, 31, 41, 47, 53, 59, 61, 71, 73, 89, 97):
    count = k3_point_count(prime)
    trace_transcendental = count - 1 - prime**2 - 20*prime
    expected = cm43_coefficient(prime)
    assert trace_transcendental == expected
    frobenius.append((prime, expected))

print(
    f"CM43H8|gross_beta={beta43}|gross_norm=43|k3_vector={tuple(v43)}"
    f"|norm=-40764|div=948|T43={T43}",
    flush=True,
)
print(
    f"CM43H8|r={r43}|s={s43}|z=+-({z_scale})*sqrt(-43)"
    f"|I2={I2}|I4={I4}|I6={I6}|I10={I10}",
    flush=True,
)
print(
    f"CM43H8|oriented_chart|m={m43}"
    f"|n=mp({-n_scale43})*sqrt(-43)|orientation_in_n=yes",
    flush=True,
)
print(
    f"CM43H8|sections=2|height_gram=[5/2,-1/2;-1/2,5/2]"
    f"|frobenius={tuple(frobenius)}",
    flush=True,
)
print(
    "CM43H8|height4_from_D9E7=1|level79=4*P1-5*P2+P3"
    "|level79_height=237/2|level79_denominator_degree=58",
    flush=True,
)
print("CM43H8|status=PASS", flush=True)

#!/usr/bin/env sage
"""Derive the compact E6+D5+A3 ambient chart for the q=80 neighbor.

Place the fibers IV* (E6), I1* (D5), and I4 (A3) at infinity, zero,
and one.  The I1* leading cubic is normalized to (x-1)^2(x+2), so

    A = T^2*(-3 + p*T + q*T^2 + r*T^3),
    B = T^3*( 2 + b1*T + ... + b4*T^4 + e*T^5).

Writing A(1)=-3*d^2 fixes r.  The four I4 discriminant jets at T=1
then determine b1,...,b4 linearly.  This leaves the expected four
parameters d,p,q,e after Weierstrass scaling.
"""

from sage.all import Matrix, PolynomialRing, QQ, vector


parameters = PolynomialRing(QQ, names=("d", "p", "q", "e"))
d, p, q, e = parameters.gens()
K = parameters.fraction_field()
KT = PolynomialRing(K, "T")
T = KT.gen()
KS = PolynomialRing(K, "s")
s = KS.gen()

r = -3 * d**2 + 3 - p - q
A = T**2 * (-3 + p * T + q * T**2 + r * T**3)

# The multiplicative branch at T=1 is the formal square-root branch
# B=2*d^3*(A/(-3*d^2))^(3/2).  Since its correction has zero constant
# term, the cubic binomial truncation gives the exact first four jets.
A_at_one = KS(A(T=1 + s))
assert A_at_one[0] == -3 * d**2
u = (A_at_one + 3 * d**2) / (-3 * d**2)
branch = 2 * d**3 * (1 + QQ(3) / 2 * u + QQ(3) / 8 * u**2 - QQ(1) / 16 * u**3)
branch_jets = vector(K, [branch[j] for j in range(4)])

# Solve the Hermite interpolation problem for the T^4,...,T^7 terms.
# The T^3 and T^8 coefficients have already been fixed to 2 and e.
jet_matrix = Matrix(
    K,
    4,
    4,
    lambda row, column: KS((1 + s) ** (4 + column))[row],
)
fixed_jets = vector(
    K, [KS(2 * (1 + s) ** 3 + e * (1 + s) ** 8)[j] for j in range(4)]
)
b1, b2, b3, b4 = jet_matrix.solve_right(branch_jets - fixed_jets)
B = T**3 * (2 + b1 * T + b2 * T**2 + b3 * T**3 + b4 * T**4 + e * T**5)

raw_discriminant = 4 * A**3 + 27 * B**2
fixed_factor = T**7 * (T - 1) ** 4
residual, remainder = raw_discriminant.quo_rem(fixed_factor)
assert remainder == 0
assert residual.degree() == 5

# Exact generic Kodaira boundaries.  The displayed open coefficients must
# stay nonzero for I1* at zero, I4 at one, and IV* at infinity.
zero_lead = (raw_discriminant // T**7)(T=0)
one_lead = (raw_discriminant // (T - 1) ** 4)(T=1)
infinity_lead = raw_discriminant[16]
assert zero_lead != 0 and one_lead != 0 and infinity_lead != 0

# Local component charts for the marked polynomial section P1.  At I1* write
# x=T*xi and y=T^2*eta.  The first exceptional cubic is
#
#     xi^3-3*xi+2 = (xi-1)^2*(xi+2).
#
# The simple root -2 is the order-two D5 class.  Following the double root
# xi=1 and blowing up once more gives eta^2=p+b1 at T=0; its two signs are
# the inverse spinor classes 1 and 3.  The I1* open coefficient is the same
# nonzero quantity, up to the fixed factor 108.
local = PolynomialRing(K, names=("xi", "eta"))
xi, eta = local.gens()
d5_first_exceptional = xi**3 - 3*xi + 2
assert d5_first_exceptional == (xi-1)**2*(xi+2)
assert zero_lead == 108*(p+b1)
d5_spinor_exceptional = eta**2-(p+b1)

# At IV* use u=1/T and (xbar,ybar)=(u^4*x(1/u),u^6*y(1/u)).  A section with
# deg(x)<=2 and deg(y)=4 has (ord_u xbar,ord_u ybar)>=(2,2); after the two
# standard blowups its outer exceptional equation is eta^2=e.  Its two signs
# are the inverse nonzero E6 classes.  In particular deg(x)=3,deg(y)=5 is not
# this component chart.
e6_outer_exceptional = eta**2-e
assert infinity_lead == 27*e**2

# At the Delta=-24 closure the frame gains exactly A1 and no MW direction
# is lost.  In this chart that is precisely a double root of the residual
# quintic, i.e. Res_T(R,R')=0.  Expanding that large four-parameter
# resultant is unnecessary.  One exact specialization with squarefree R
# proves that the resultant is not identically zero.
specialization = {d: QQ(2), p: QQ(1), q: QQ(2), e: QQ(3)}
QZ = PolynomialRing(QQ, "z")
residual_special = QZ(
    [QQ(coefficient.subs(specialization)) for coefficient in residual]
)
assert residual_special.degree() == 5
assert residual_special.gcd(residual_special.derivative()) == 1

print("Q80NORMAL|A={}".format(A), flush=True)
print(
    "Q80NORMAL|B_coefficients=b1:{};b2:{};b3:{};b4:{};b5:{}".format(
        b1.factor(), b2.factor(), b3.factor(), b4.factor(), e
    ),
    flush=True,
)
print(
    "Q80NORMAL|Delta_fixed=T^7*(T-1)^4|residual_degree=5|"
    "fibers=I1*,I4,IV*|ADE=D5+A3+E6|ambient_dimension=4",
    flush=True,
)
print(
    "Q80NORMAL|Delta=-24_boundary=disc_T(residual_quintic)=0|"
    "enhancement=plus_A1",
    flush=True,
)
print(
    "Q80LOCAL|I1star_first_exceptional={}|double_root=1|simple_root=-2|"
    "spinor_exceptional={}".format(d5_first_exceptional, d5_spinor_exceptional),
    flush=True,
)
print(
    "Q80LOCAL|IVstar_nonzero_profile=ord_xbar>=2,ord_ybar=2|"
    "outer_exceptional={}".format(e6_outer_exceptional),
    flush=True,
)
print("Q80NORMAL|status=PASS", flush=True)

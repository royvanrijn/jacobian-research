#!/usr/bin/env sage
"""Verify the exact discriminant-24 CM anchor in the q=60 chart.

The height-4 section chart is normalized by x_P(0)=1 and by moving the A2
fiber to T=1.  The rational point below was recovered from the nonsingular
solution (16,12,8) modulo 31 by multivariate Hensel lifting and rational
reconstruction.
"""

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, rational_reconstruction, vector


R = PolynomialRing(QQ, "T")
T = R.gen()

c1 = QQ(20) / 9
c2 = QQ(334) / 729
c3 = QQ(68) / 729
d = -QQ(1) / 27

xP = 1 + c1*T + c2*T**2 + c3*T**3 + d**2*T**4
yP = (
    1 + QQ(3)/2*c1*T
    + (QQ(3)/2*c2 + QQ(3)/8*c1**2)*T**2
    + (QQ(3)/2*d*c1 + QQ(3)/4*c3*c2/d
       - QQ(1)/16*c3**3/d**3)*T**3
    + (QQ(3)/2*d*c2 + QQ(3)/8*c3**2/d)*T**4
    + QQ(3)/2*d*c3*T**5 + d**3*T**6
)

a0 = -QQ(8192) / 19683
a1 = -QQ(28672) / 19683
b0 = 0
b1 = -QQ(4980736) / 14348907
b2 = -QQ(8912896) / 14348907
b3 = -QQ(262144) / 14348907

a4 = T**3*(a0 + a1*T)
a6 = T**4*(b0 + b1*T + b2*T**2 + b3*T**3)
discriminant = -16*(4*a4**3 + 27*a6**2)

assert xP == (T**4 + 68*T**3 + 334*T**2 + 1620*T + 729) / 729
assert yP == (-T**6 - 102*T**5 - 2235*T**4 + 188*T**3
              + 49977*T**2 + 65610*T + 19683) / 19683
assert yP**2 == xP**3 + a4*xP + a6

assert a4 == -QQ(4096)/19683*T**3*(7*T + 2)
assert a6 == -QQ(262144)/14348907*T**5*(T**2 + 34*T + 19)
assert discriminant == (
    -QQ(1099511627776)/7625597484987
    * T**9*(T - 1)**3*(T**2 + 71*T + 32)
)

# Kodaira fibers: III* at zero, I3 at one, II* at infinity, and two I1.
assert (a4.valuation(T), a6.valuation(T), discriminant.valuation(T)) == (3, 5, 9)
assert a4(1) != 0 and a6(1) != 0
assert (8-a4.degree(), 12-a6.degree(), 24-discriminant.degree()) == (4, 5, 10)
assert (T**2 + 71*T + 32).discriminant() == 17**3

# P1 meets the nonsingular locus at both reducible finite fibers, hence the
# identity components; with P1.O=0 its Shioda height is exactly four.
assert (xP(0), yP(0)) == (1, 1)
assert yP(1) != 0
height_P1 = 4
root_determinant = 2 * 3  # E7 + A2; E8 is unimodular.
assert root_determinant * height_P1 == 24

# Replay the finite-field seed, Hensel lift, and rational reconstruction rather
# than treating the displayed parameters as unexplained input.
S = PolynomialRing(QQ, names=("C1", "C3", "D"))
C1, C3, DD = S.gens()
K = S.fraction_field()
C2 = (C1**2*DD**2 + 2*C1*C3*DD + 48*DD**3 + 3*C3**2) / (24*DD**2)
KT = PolynomialRing(K, "U")
U = KT.gen()
XP = 1 + C1*U + C2*U**2 + C3*U**3 + DD**2*U**4
YP = (
    1 + QQ(3)/2*C1*U
    + (QQ(3)/2*C2 + QQ(3)/8*C1**2)*U**2
    + (QQ(3)/2*DD*C1 + QQ(3)/4*C3*C2/DD
       - QQ(1)/16*C3**3/DD**3)*U**3
    + (QQ(3)/2*DD*C2 + QQ(3)/8*C3**2/DD)*U**4
    + QQ(3)/2*DD*C3*U**5 + DD**3*U**6
)
section_residual = YP**2 - XP**3
section_coefficients = [section_residual[index] for index in range(13)]
A0 = section_coefficients[3]
A1 = section_coefficients[8] / DD**2
B1 = section_coefficients[5] - A0*C2 - A1*C1
B2 = section_coefficients[6] - A0*C3 - A1*C2
B3 = section_coefficients[7] - A0*DD**2 - A1*C3
residual_discriminant = (
    -16*(4*(U**3*(A0+A1*U))**3
         + 27*(U**5*(B1+B2*U+B3*U**2))**2) // U**9
)

ell = C1*DD - C3
equations = []
for value in (
    residual_discriminant(1),
    residual_discriminant.derivative()(1),
    residual_discriminant.derivative(2)(1),
):
    numerator = S(value.numerator())
    numerator = S(numerator / numerator.content())
    equations.append(S(numerator // ell**6))

prime = 31
seed = (ZZ(16), ZZ(12), ZZ(8))
assert all(ZZ(equation(*seed)) % prime == 0 for equation in equations)
jacobian = matrix(S, [
    [equation.derivative(variable) for variable in (C1, C3, DD)]
    for equation in equations
])
assert ZZ(jacobian(*seed).det()) % prime == 20

values = list(seed)
modulus = ZZ(prime)
for _ in range(24):
    function_values = vector(ZZ, [ZZ(equation(*values)) for equation in equations])
    assert all(value % modulus == 0 for value in function_values)
    rhs = vector(GF(prime), [-ZZ(value // modulus) for value in function_values])
    jacobian_mod_prime = matrix(GF(prime), [
        [ZZ(entry) % prime for entry in row]
        for row in jacobian(*values).rows()
    ])
    digit = jacobian_mod_prime.solve_right(rhs)
    values = [
        (values[index] + modulus*ZZ(digit[index])) % (modulus*prime)
        for index in range(3)
    ]
    modulus *= prime

reconstructed = tuple(rational_reconstruction(value, modulus) for value in values)
assert reconstructed == (c1, c3, d)
assert all(equation(*reconstructed) == 0 for equation in equations)

print(
    "Q60CM24|c1=20/9|c2=334/729|c3=68/729|d=-1/27|"
    "section_height=4|seed_mod31=16,12,8|jacobian_det_mod31=20"
)
print(
    "Q60CM24|a4=-4096/19683*T^3*(7*T+2)|"
    "a6=-262144/14348907*T^5*(T^2+34*T+19)"
)
print(
    "Q60CM24|fibers=II*+III*+I3+2I1|disc=24|"
    "residual_quadratic_discriminant=17^3|status=PASS"
)

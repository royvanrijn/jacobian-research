#!/usr/bin/env sage
"""Verify the two explicit CM24 binary-quartic moves in the Q80 low-q corridor.

The exact field is K=QQ(sqrt(-6)).  Starting from the known CM24 second Q80
child, this script verifies

  * the new q6 marked-chord coordinate and D8+D6+2A1 child;
  * the orbit-424 q4 coordinate using the specialized 2-torsion section and
    the A7+A7 child.

This is a CM24 equation certificate, not a generic-family lift.
"""

from sage.all import PolynomialRing, QuadraticField

K = QuadraticField(-6, "s")
s = K.gen()

# ---------------------------------------------------------------------------
# q6 escape from the CM24 second child.
# ---------------------------------------------------------------------------
VR = PolynomialRing(K, "V")
V = VR.gen()
WR = PolynomialRing(VR, "u")
u = WR.gen()
W = u-K(27)/2

old_A = (
    -27*W**6 + 59049*W**4 + K(13286025)/8*W**3
    + K(129140163)/8*W**2 + K(1162261467)/8*W
    - K(10460353203)/64
)
# Correct the W coefficient to the pinned CM24 model.
old_A = (
    -27*W**6 + 59049*W**4 + K(13286025)/8*W**3
    + K(129140163)/8*W**2 + K(1162261467)/32*W
    - K(10460353203)/64
)
old_B = (
    54*W**9 - 177147*W**7 - K(97253703)/8*W**6
    - K(7360989291)/16*W**5 - K(331244518095)/32*W**4
    - K(4487491524087)/32*W**3 - K(144886352214753)/128*W**2
    - K(1303977169932777)/256*W - K(5147278302366225)/512
)
Qx = (
    -K(8)/27*W**4 + 22*W**3 - K(243)/2*W**2
    + 729*W - K(492075)/8
)
Qy = s*(
    K(16)/243*W**6 - K(22)/3*W**5 + K(333)/2*W**4
    - K(2025)/4*W**3 + K(190269)/4*W**2
    - K(177147)/16*W + K(199290375)/32
)
assert Qy**2 == Qx**3+old_A*Qx+old_B

lam = -K(39)/4*s
zq = V*u**2 + 162*s + lam*u
branch_q6 = (
    (zq**2-3*Qx)**2
    - 4*(2*zq*Qy+3*Qx**2+old_A)
)
assert all(branch_q6[index] == 0 for index in range(4))
quartic_q6, remainder = branch_q6.quo_rem(u**4)
assert remainder == 0 and quartic_q6.degree() == 4

q0, q1, q2, q3, q4 = [VR(quartic_q6[index]) for index in range(5)]
I = 12*q4*q0-3*q3*q1+q2**2
J = 72*q4*q2*q0+9*q3*q2*q1-27*q4*q1**2-27*q3**2*q0-2*q2**3
A6 = -27*I
B6 = -27*J
D6 = -16*(4*A6**3+27*B6**2)

assert A6 == (
    68024448*V**6 - 160849476*s*V**5 - 916676676*V**4
    + 446944068*s*V**3 + 713569986*V**2 - 98644392*s*V
    - 33250608
)
assert B6 == (
    -88159684608*s*V**9 - 1876148288064*V**8
    + 2875406171544*s*V**7 + 15062360165208*V**6
    - 8283115953504*s*V**5 - 17889909720408*V**4
    + 4221255729708*s*V**3 + 3780760503600*V**2
    - 324124340832*s*V - 72965135232
)

f_i2a = V-K(13)/18*s
f_i2b = V-K(5)/9*s
f_i4s = V-K(2)/9*s
f_res = V**2-K(41)/72*s*V-K(481)/864
assert D6.valuation(f_i2a) == 2
assert D6.valuation(f_i2b) == 2
assert D6.valuation(f_i4s) == 10
assert D6.valuation(f_res) == 1
assert A6.valuation(f_i4s) == 2
assert B6.valuation(f_i4s) == 3
assert (8-A6.degree(), 12-B6.degree(), 24-D6.degree()) == (2, 3, 8)

print(
    "Q80LOWQCM24|q6|coordinate="
    "V=(zQ-162*s+(39/4)*s*(W+27/2))/(W+27/2)^2|"
    "fibres=I4*+I2*+2I2+2I1|ADE=D8+D6+2A1|MW=2|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# orbit-424 q4 from the explicit q6 model.
# ---------------------------------------------------------------------------
T = 1944*s*V**3 + 12150*V**2 - 4401*s*V - 3036
assert T**3+A6*T+B6 == 0

UR = PolynomialRing(K, "U")
U = UR.gen()
VUR = PolynomialRing(UR, "V")
VV = VUR.gen()

A6u = VUR(A6(VV))
B6u = VUR(B6(VV))
Tu = VUR(T(VV))
v0 = K(13)/18*s
zT = (VV-v0)*U
branch_q4 = zT**4 - 6*Tu*zT**2 - 3*Tu**2 - 4*A6u
quartic_q4, remainder = branch_q4.quo_rem((VV-v0)**2)
assert remainder == 0 and quartic_q4.degree() == 4

r0, r1, r2, r3, r4 = [UR(quartic_q4[index]) for index in range(5)]
I4 = 12*r4*r0-3*r3*r1+r2**2
J4 = 72*r4*r2*r0+9*r3*r2*r1-27*r4*r1**2-27*r3**2*r0-2*r2**3
A4 = -27*I4
B4 = -27*J4
D4 = -16*(4*A4**3+27*B4**2)

assert A4 == (
    -27*U**8 - 4251528*U**6 - 81616583016*U**4
    - 27113235502176*U**2 - 9882774340543152
)
assert B4 == (
    54*U**12 + 12754584*U**10 + 746946702792*U**8
    + 6181817694496128*U**6 + 5751774666196114464*U**4
    + 1556181178759286886528*U**2 + 378152026438506713426304
)

octic = (
    U**8 + K(314199)/2*U**6 + 1905215985*U**4
    + 509070522546*U**2 + 183014339639688
)
assert D4.valuation(U) == 8
assert D4.valuation(octic) == 1
assert octic.gcd(octic.derivative()) == 1
assert (8-A4.degree(), 12-B4.degree(), 24-D4.degree()) == (0, 0, 8)

print(
    "Q80LOWQCM24|q4-orbit424|"
    "Tx=1944*s*V^3+12150*V^2-4401*s*V-3036|Ty=0|"
    "coordinate=U=Y/((X-Tx)*(V-13*s/18))|"
    "fibres=2I8+8I1|ADE=2A7|MW=4|status=PASS",
    flush=True,
)

print("Q80LOWQCM24|status=PASS_LOWQ_CM24_EQUATIONS", flush=True)

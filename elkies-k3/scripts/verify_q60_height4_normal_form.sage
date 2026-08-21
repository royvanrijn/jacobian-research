#!/usr/bin/env sage
"""Derive a three-parameter q=60 ambient chart from its height-4 section.

We normalize the polynomial x-coordinate to

    x_P = 1 + c1*T + c2*T^2 + c3*T^3 + T^4.

The equation at zero determines the first three y-jets, while the II* fiber
at infinity determines the last four.  The relative leading sign epsilon is
+/-1.  The remaining six coefficients are then read off linearly from
y_P^2-x_P^3.
"""

from sage.all import PolynomialRing, QQ


S = PolynomialRing(QQ, names=("c1", "c2", "c3"))
c1, c2, c3 = S.gens()
R = PolynomialRing(S, "T")
T = R.gen()

x = 1 + c1*T + c2*T**2 + c3*T**3 + T**4


def chart(epsilon):
    y = (
        1
        + QQ(3)/2*c1*T
        + (QQ(3)/2*c2 + QQ(3)/8*c1**2)*T**2
        + epsilon*(QQ(3)/2*c1 + QQ(3)/4*c3*c2
                   - QQ(1)/16*c3**3)*T**3
        + epsilon*(QQ(3)/2*c2 + QQ(3)/8*c3**2)*T**4
        + epsilon*QQ(3)/2*c3*T**5
        + epsilon*T**6
    )
    residual = y**2 - x**3

    # The endpoint jets leave exactly the coefficient window T^3,...,T^8.
    assert all(residual[index] == 0 for index in (0, 1, 2, 9, 10, 11, 12))

    a0 = residual[3]
    a1 = residual[8]
    b0 = residual[4] - a0*c1 - a1
    b1 = residual[5] - a0*c2 - a1*c1
    b2 = residual[6] - a0*c3 - a1*c2
    b3 = residual[7] - a0 - a1*c3

    a4 = T**3*(a0 + a1*T)
    a6 = T**4*(b0 + b1*T + b2*T**2 + b3*T**3)
    assert y**2 == x**3 + a4*x + a6

    # The CM enhancement hyperplanes are already visible in this chart.
    ell = c3 - epsilon*c1
    assert a0 % ell == 0
    assert b0 % ell**2 == 0
    assert b3 % ell**3 == 0
    return y, (a0, a1, b0, b1, b2, b3), ell


for epsilon in (1, -1):
    y, coefficients, ell = chart(epsilon)
    a0, a1, b0, b1, b2, b3 = coefficients
    print(
        "Q60HEIGHT4|epsilon={}|x={}|y={}|a0={}|a1={}|"
        "b0={}|b1={}|b2={}|b3={}".format(
            epsilon, x, y, a0, a1, b0, b1, b2, b3
        )
    )
    print(
        "Q60HEIGHT4|epsilon={}|Delta3_boundary={}=0|"
        "a0_divisibility=1|b0_divisibility=2|b3_divisibility=3".format(
            epsilon, ell
        )
    )

print("Q60HEIGHT4|parameters=3|section_PO=0|status=PASS")

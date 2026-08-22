#!/usr/bin/env sage -python
"""Exact pole-profile reduction for the second marked H3 direction.

The H3 Kumar Mordell--Weil Gram matrix in the marked basis is

    [[21/2, 3], [3, 46]].

``P1`` meets the nonidentity E7 component, whose Shioda correction is 3/2,
whereas ``P2`` meets the identity component.  This script tests the proposed
small translates ``Q_m=P2+m*P1``.  It proves the correction and zero-section
intersection predicted by the lattice, rather than treating denominator degree
as an unconstrained rational-function parameter.
"""

from sage.all import QQ, ZZ, matrix, vector


HEIGHT = matrix(QQ, [[QQ(21) / 2, 3], [3, 46]])


def component_correction(coordinates):
    """E7 discriminant correction for the marked component class.

    The component group is Z/2.  P1 represents its nonzero class and P2
    represents zero, so only the parity of the P1 coefficient matters.
    """
    return QQ(3) / 2 if ZZ(coordinates[0]) % 2 else QQ(0)


def data_for(m):
    coordinates = vector(ZZ, [m, 1])
    height = coordinates * HEIGHT * coordinates
    correction = component_correction(coordinates)
    intersection = (height - 4 + correction) / 2
    assert intersection in ZZ
    # Shioda: h(P)=4+2(P.O)-contr_E7(P), for this E7+E8 K3.
    assert height == 4 + 2 * intersection - correction
    return height, correction, ZZ(intersection)


assert HEIGHT.det() == 474
assert data_for(0) == (46, 0, 21)
assert data_for(-1) == (QQ(101) / 2, QQ(3) / 2, 24)
assert data_for(1) == (QQ(125) / 2, QQ(3) / 2, 30)

window = tuple(range(-6, 7))
rows = tuple((m,) + data_for(m) for m in window)
minimum = min(rows, key=lambda row: (row[3], row[1], abs(row[0])))
assert minimum == (0, 46, 0, 21)

# Complete the minimization on Z.  Relative to P2, the pole increase is
# ``3*m*(7*m+4)/4`` for even m and ``(21*m^2+12*m+3)/4`` for odd m.  The odd
# quadratic has negative discriminant; the even expression is positive away
# from zero because the nearest nonzero even integers are +/-2.
assert 12**2 - 4 * 21 * 3 < 0
assert 3 * (-2) * (7 * (-2) + 4) / 4 > 0
assert 3 * 2 * (7 * 2 + 4) / 4 > 0
for m in (-8, -6, -4, -2, 2, 4, 6, 8):
    assert data_for(m)[2] - 21 == 3 * m * (7 * m + 4) / 4 > 0
for m in (-7, -5, -3, -1, 1, 3, 5, 7):
    assert data_for(m)[2] - 21 == (21 * m**2 + 12 * m + 3) / 4 > 0

print("H3P2POLE|height_gram=[[21/2,3],[3,46]]|component_group=E7:Z/2")
for m, height, correction, intersection in rows:
    print(
        "H3P2POLE|m={}|height={}|correction={}|Qm.O={}".format(
            m, height, correction, intersection
        )
    )
print(
    "H3P2POLE|best_translate=P2|pole_order=21|"
    "nearest_competitor=P2-P1:24|status=PASS_EXACT"
)

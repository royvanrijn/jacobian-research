#!/usr/bin/env sage
"""Reconstruct the CM24 third-q12 Jacobian over GF(73)(V).

The canonical simple infinity branch ``xi=-6`` is used as origin in
``analyze_q80_third_q12_cm24_weierstrass_gf73.sage``.  That exact
Brill--Noether calculation supplies the 56 pinned ``(V,j)`` fibers below.
The first 49 determine a unique bidegree-(24,24) rational function; the final
seven are withheld validation fibers.

The recovered numerator is a cube ``c4^3`` and
``c4^3-1728*Delta`` is a square ``c6^2``.  The resulting short Weierstrass
model has discriminant fibers ``2 I7 + 3 I2 + 4 I1``, exactly the transported
CM24 root system ``2A6+3A1``.  This is an exact finite-characteristic
Jacobian/marking certificate, not yet its characteristic-zero lift.
"""

from sage.all import GF, Matrix, PolynomialRing, vector


finite = GF(73)
polynomial_ring = PolynomialRing(finite, "V")
V = polynomial_ring.gen()

samples = (
    (3, 21), (4, 45), (7, 36), (8, 27), (9, 53), (10, 47),
    (11, 3), (12, 36), (13, 15), (14, 53), (15, 3), (16, 34),
    (17, 26), (18, 45), (19, 15), (20, 50), (21, 47), (22, 34),
    (23, 21), (24, 71), (25, 49), (26, 47), (28, 45), (29, 15),
    (30, 51), (31, 26), (32, 36), (33, 47), (34, 0), (35, 59),
    (37, 0), (38, 19), (39, 36), (40, 38), (41, 3), (42, 11),
    (44, 50), (45, 1), (46, 30), (47, 59), (48, 30), (49, 35),
    (50, 11), (51, 19), (52, 51), (54, 21), (55, 15), (57, 33),
    (58, 20),
    # Withheld after the interpolation kernel is fixed.
    (59, 35), (60, 56), (61, 38), (62, 56), (63, 59), (64, 3),
    (65, 52),
)
training = samples[:49]
withheld = samples[49:]
assert len(training) == 49 and len(withheld) == 7


def interpolation_row(value, j_value):
    value = finite(value)
    j_value = finite(j_value)
    return [value**index for index in range(25)] + [
        -j_value*value**index for index in range(25)
    ]


interpolation = Matrix(
    finite, [interpolation_row(value, j_value) for value, j_value in training]
)
assert interpolation.nrows() == 49 and interpolation.ncols() == 50
assert interpolation.rank() == 49
kernel = interpolation.right_kernel().basis()[0]
kernel /= kernel[0]
expected_kernel = vector(
    finite,
    [
        1, 53, 22, 17, 58, 36, 37, 57, 32, 17, 6, 49, 58, 33, 23,
        30, 50, 51, 58, 33, 2, 25, 71, 1, 64,
        1, 29, 27, 37, 20, 35, 39, 40, 4, 9, 44, 62, 1, 3, 62,
        50, 50, 40, 10, 4, 62, 68, 7, 3, 10,
    ],
)
assert kernel == expected_kernel

j_numerator = polynomial_ring(list(kernel[:25]))
discriminant = polynomial_ring(list(kernel[25:]))
assert all(
    discriminant(value) != 0
    and j_numerator(value) == finite(j_value)*discriminant(value)
    for value, j_value in withheld
)

c4 = 4*(V+36)*(V+39)*(
    V**6+25*V**5+26*V**4+3*V**3+24*V**2+23*V+31
)
c6 = (
    31*V**12+38*V**10+2*V**9+29*V**8+22*V**7+21*V**6
    +64*V**5+23*V**4+16*V**3+27*V**2+53*V+5
)
A = 6*V**8+16*V**7+47*V**6+33*V**5+58*V**4+2*V**3+63*V**2+17*V+23
B = (
    33*V**12+64*V**10+61*V**9+45*V**8+14*V**7+20*V**6
    +54*V**5+8*V**4+50*V**3+57*V**2+47*V+43
)

assert j_numerator == c4**3
assert c4**3-c6**2 == finite(1728)*discriminant
assert -finite(48)*A == c4
assert -finite(864)*B == c6
assert -finite(16)*(4*A**3+27*B**2) == discriminant

expected_discriminant = 10*(V+37)*(V+46)*(V+17)**2*(V+30)**2*(V+68)**2*(V+20)**7*(V+67)**7*(V**2+20*V+67)
assert discriminant == expected_discriminant
assert sum((1, 1, 2, 2, 2, 7, 7, 2)) == 24
root_rank = 2*6+3
assert root_rank == 15
picard_rank = 20
trivial_rank = 2+root_rank
mw_rank = picard_rank-trivial_rank
assert (trivial_rank, mw_rank) == (17, 3)

print(
    "Q80THIRDGF73JACOBIAN|prime=73|origin=xi:-6(simple)|"
    "interpolation=49x50|rank=49|withheld=7|bidegree=24,24|"
    f"A={A}|B={B}",
    flush=True,
)
print(
    "Q80THIRDGF73JACOBIAN|"
    f"Delta={discriminant.factor()}|fibers=2I7+3I2+4I1|"
    "roots=2A6+3A1|rho=20|trivial_rank=17|MW=3|"
    "status=PASS_INTERPOLATED_JACOBIAN",
    flush=True,
)

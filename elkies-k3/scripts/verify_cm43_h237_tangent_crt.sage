#!/usr/bin/env sage
"""Certify the reconstructed CM(-43) canonical-slice cubic.

The companion finite-field verifier ``compute_cm43_h237_tangent.sage``
derives the third-order obstruction from the exact level-79 section after
choosing a canonical right inverse for the second correction.  This script
combines twenty of those normalized slice cubics by CRT, performs rational
reconstruction with a strict uniqueness bound, and factors the result.

The companion script also proves that adding first-order kernel vectors to
the second correction absorbs this whole cubic.  Thus the rational factor
below is a reproducible slice diagnostic, *not* a Humbert-237 tangent.
"""

from sage.all import PolynomialRing, QQ, ZZ, crt, gcd, lcm, prod


MODULAR_CUBICS = (
    (1000003, (1, 886334, 452938, 837495)),
    (1000033, (1, 861594, 637991, 614137)),
    (1000037, (1, 492321, 993513, 477233)),
    (1000039, (1, 823399, 198966, 234386)),
    (1000081, (1, 49221, 704315, 581804)),
    (1000099, (1, 579039, 638745, 718538)),
    (1000117, (1, 653206, 21363, 49389)),
    (1000121, (1, 173529, 279038, 905092)),
    (1000133, (1, 15124, 726307, 106418)),
    (1000151, (1, 408551, 701648, 668129)),
    (1000159, (1, 221492, 237738, 600074)),
    (1000171, (1, 642976, 696747, 125205)),
    (1000183, (1, 584535, 824299, 643656)),
    (1000187, (1, 806390, 43554, 631380)),
    (1000193, (1, 695084, 250624, 959991)),
    (1000199, (1, 831157, 260138, 734807)),
    (1000211, (1, 105002, 694027, 222876)),
    (1000213, (1, 871189, 687964, 323194)),
    (1000231, (1, 981848, 557590, 267435)),
    (1000249, (1, 820735, 430183, 759638)),
)

EXPECTED = (
    QQ(1),
    -QQ(32777647185971477137326047735483125)
    / QQ(5537010211609548283434042242558544),
    -QQ(62256599453976317430685478929919140625)
    / QQ(4252423842516133081677344442284961792),
    QQ(16726681288628536079919155997729849853515625)
    / QQ(446402445291973586382160910173306149076992),
)

primes = tuple(ZZ(prime) for prime, _ in MODULAR_CUBICS)
modulus = prod(primes)
assert modulus.nbits() == 399
reconstructed = []
for column in range(4):
    residue = crt(
        [ZZ(coefficients[column]) for _, coefficients in MODULAR_CUBICS],
        list(primes),
    )
    value = ZZ(residue).rational_reconstruction(modulus)
    # If |a| <= A and 0 < b <= B with 2*A*B < M, the reconstruction a/b
    # modulo M is unique in that rectangle.
    assert 2*abs(value.numerator())*value.denominator() < modulus
    reconstructed.append(value)
assert tuple(reconstructed) == EXPECTED

denominator = lcm(value.denominator() for value in reconstructed)
primitive = [ZZ(denominator*value) for value in reconstructed]
content = gcd(primitive)
primitive = tuple(value//content for value in primitive)
assert primitive == (
    446402445291973586382160910173306149076992,
    -2642585311483808090087377295677498775040000,
    -6535448784280617898603638836147191706250000,
    16726681288628536079919155997729849853515625,
)

R = PolynomialRing(QQ, "z")
z = R.gen()
cubic = sum(reconstructed[index]*z**(3-index) for index in range(4))
rational_slope = QQ(223593125)/30934224
assert cubic(rational_slope) == 0
quadratic = cubic//(z-rational_slope)
assert quadratic.degree() == 2 and quadratic.is_irreducible()

print(
    "CM43H237SLICECRT|primes=20|modulus_bits=399|unique_reconstruction=1",
    flush=True,
)
print(
    "CM43H237SLICECRT|primitive_cubic=" + ",".join(map(str, primitive)),
    flush=True,
)
print(
    f"CM43H237SLICECRT|rational_slice_factor={rational_slope}"
    "|quadratic_slice_factor=irreducible|status=GAUGE_ABSORBED_NOT_A_TANGENT",
    flush=True,
)

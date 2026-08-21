#!/usr/bin/env sage
"""Verify and resolve the q=80 discriminant-24 seed over GF(7).

The companion exhaustive scanner leaves one normalized surface with the two
polynomial section directions that continue to the generic rank-19 branch.
It also finds six one-pole classes with the correct identity/nonidentity
pattern.  This script uses the exact function-field group law and pole counts
to select the signs/classes forced by the transported Mordell--Weil Gram.
"""

from sage.all import EllipticCurve, GF, PolynomialRing


PROTOCOL = "Q80CM7VERIFY"
field = GF(7)
polynomials = PolynomialRing(field, "T")
T = polynomials.gen()
functions = polynomials.fraction_field()


def polynomial(coefficients):
    return polynomials(list(map(field, coefficients)))


A = polynomial((0, 0, -3, 4, 3, 4))
B = polynomial((0, 0, 0, 2, 0, 2, 6, 0, 2))
curve = EllipticCurve(functions, [0, 0, 0, A, B])

P1 = curve(
    polynomial((0, 1, 2, 0, 0)),
    polynomial((0, 0, 2, 1, 4, 0, 0)),
)
P2 = curve(
    polynomial((1, 2, 2, 6, 4)),
    polynomial((1, 3, 3, 2, 3, 3, 6)),
)

raw_candidates = (
    (0, (4, 3, 2, 1, 5, 1, 1), (1, 2, 2, 4, 1, 1, 2, 0, 2, 6)),
    (3, (2, 0, 6, 4, 4, 6, 4), (1, 0, 3, 1, 1, 2, 5, 3, 4, 1)),
    (6, (1, 4, 3, 0, 6, 1, 4), (1, 6, 2, 5, 3, 4, 1, 2, 3, 1)),
    (None, (2, 2, 6, 0, 0, 0, 0), (1, 5, 1, 3, 4, 0, 0, 0, 0, 0)),
    (None, (2, 4, 4, 3, 5, 2, 4), (1, 3, 5, 4, 6, 2, 0, 0, 1, 6)),
    (None, (4, 3, 3, 0, 0, 0, 0), (1, 2, 6, 2, 3, 0, 0, 0, 0, 0)),
)


def candidate_point(record):
    pole, numerator_x, numerator_y = record
    denominator = polynomials.one() if pole is None else T - field(pole)
    return curve(
        functions(polynomial(numerator_x)) / denominator**2,
        functions(polynomial(numerator_y)) / denominator**3,
    )


candidates = tuple(map(candidate_point, raw_candidates))


def zero_intersection(point):
    """Recover P.O from the reduced x-coordinate on this K3 model."""
    if point.is_zero():
        raise ValueError("the zero section has no finite self pole count")
    x_coordinate = functions(point[0])
    numerator = x_coordinate.numerator()
    denominator = x_coordinate.denominator()
    if denominator.degree() % 2:
        raise AssertionError("x denominator is not a square-degree divisor")
    finite = denominator.degree() // 2
    excess = numerator.degree() - denominator.degree() - 4
    if excess < 0:
        excess = 0
    if excess % 2:
        raise AssertionError("the infinity x-pole has odd order")
    return finite + excess // 2


assert zero_intersection(P1) == zero_intersection(P2) == 0
candidate_poles = tuple(zero_intersection(point) for point in candidates)
print(f"{PROTOCOL}|candidate_poles={candidate_poles}", flush=True)
assert candidate_poles == (1, 1, 1, 0, 1, 0)

sign_hits = []
for sign1 in (1, -1):
    for sign2 in (1, -1):
        oriented1 = sign1 * P1
        oriented2 = sign2 * P2
        pole12 = zero_intersection(oriented1 + oriented2)
        if pole12 != 2:
            continue
        for index, point in enumerate(candidates, 1):
            if candidate_poles[index - 1] != 1:
                continue
            for sign3 in (1, -1):
                oriented3 = sign3 * point
                pole13 = zero_intersection(oriented1 + oriented3)
                pole23 = zero_intersection(oriented2 + oriented3)
                print(
                    f"{PROTOCOL}|candidate={index}|signs={sign1},{sign2},{sign3}"
                    f"|sum_poles={pole12},{pole13},{pole23}",
                    flush=True,
                )
                if (pole13, pole23) == (1, 5):
                    sign_hits.append((sign1, sign2, sign3, index))

print(f"{PROTOCOL}|selected={sign_hits}", flush=True)
assert sign_hits

# The two selected-looking third-section representatives are not distinct
# marked deformations.  With the convention used by the formal lift,
# G3=-candidate_3.  The alternative hit has
#
#     G2'=-G2,  G3'=-candidate_5=G3-G2,
#
# an integral basis automorphism that leaves the transported level-79 class
# -3*G1-2*G2+4*G3 unchanged.  Thus candidate 5 cannot repair the failed Q79
# cover by selecting a different CM seed.
G1 = P1
G2 = P2
G3 = -candidates[2]
G2_alternative = -P2
G3_alternative = -candidates[4]
assert G3_alternative == G3-G2
Q79 = -3*G1-2*G2+4*G3
Q79_alternative = -3*G1-2*G2_alternative+4*G3_alternative
assert Q79_alternative == Q79
print(
    f"{PROTOCOL}|candidate3_candidate5_relation=G3alt=G3-G2"
    "|G2alt=-G2|Q79_invariant=1|status=PASS",
    flush=True,
)
print(
    f"{PROTOCOL}|surface=d:3,p:4,q:3,e:2|extra_I2=4"
    f"|polynomial_sections=18|coarse_one_pole=6|status=PASS",
    flush=True,
)

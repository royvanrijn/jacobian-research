#!/usr/bin/env sage
"""Probe the first IV*-denominator saturation of the CM24 q12 module.

The polynomial submodule for ``<1,X,z_Q>`` contains no q12 pencil.  At the
CM24 second child the finite additive fiber is IV* at ``W=-27/2``, which is
``W=-3`` modulo seven.  This bounded experiment permits one common local
denominator there:

    V = z_Q + ((a1*W+a0)*X + b2*W^2+b1*W)/(W+3).

It clears the denominator symbolically once, forms the cubic discriminant,
and then checks all 7^4 coefficient tuples.  The reported branch-support
degree counts distinct finite discriminant factors after undoing the eighth
power introduced by clearing a cubic, plus infinity.  It is a cheap filter
for the saturated local module, not a proof that a survivor has genus one or
the transported ``2A6+3A1`` fibers.
"""

from collections import Counter

from sage.all import GF, PolynomialRing


field = GF(7)
coefficient_ring = PolynomialRing(
    field, names=("a1", "a0", "b2", "b1", "V")
)
a1, a0, b2, b1, V = coefficient_ring.gens()
old_base_ring = PolynomialRing(coefficient_ring, "W")
W = old_base_ring.gen()
x_ring = PolynomialRing(old_base_ring, "X")
X = x_ring.gen()

A = (
    -27*W**6+59049*W**4+field(13286025)/8*W**3
    + field(129140163)/8*W**2+field(1162261467)/32*W
    - field(10460353203)/64
)
B = (
    54*W**9-177147*W**7-field(97253703)/8*W**6
    - field(7360989291)/16*W**5-field(331244518095)/32*W**4
    - field(4487491524087)/32*W**3-field(144886352214753)/128*W**2
    - field(1303977169932777)/256*W-field(5147278302366225)/512
)
Qx = (
    -field(8)/27*W**4+22*W**3-field(243)/2*W**2+729*W
    - field(492075)/8
)
Qy = (
    field(16)/243*W**6-field(22)/3*W**5+field(333)/2*W**4
    - field(2025)/4*W**3+field(190269)/4*W**2
    - field(177147)/16*W+field(199290375)/32
)
assert Qy**2 == Qx**3+A*Qx+B

denominator = W+3
a_numerator = a1*W+a0
b_numerator = b2*W**2+b1*W
cleared_z_numerator = denominator*V-a_numerator*X-b_numerator
cleared_equation = (
    (cleared_z_numerator*(X-Qx)-denominator*Qy)**2
    - denominator**2*(X**3+A*X+B)
)
quotient, remainder = cleared_equation.quo_rem(X-Qx)
assert remainder == 0 and quotient.degree() == 3

# Discriminant of q3*X^3+q2*X^2+q1*X+q0.
q0, q1, q2, q3 = [quotient[index] for index in range(4)]
discriminant = (
    q2**2*q1**2-4*q3*q1**3-4*q2**3*q0-27*q3**2*q0**2
    + 18*q3*q2*q1*q0
)
assert discriminant.degree() >= 23

parameter_ring = PolynomialRing(field, "v")
v = parameter_ring.gen()
special_old_base = PolynomialRing(parameter_ring, "w")
w = special_old_base.gen()
special_denominator = w+3


def valuation(poly, factor):
    result = 0
    while poly % factor == 0:
        poly //= factor
        result += 1
    return result


distribution = Counter()
candidates = []
for a1_value in field:
    for a0_value in field:
        if not a1_value and not a0_value:
            continue
        for b2_value in field:
            for b1_value in field:
                specialization = coefficient_ring.hom(
                    [a1_value, a0_value, b2_value, b1_value, v],
                    parameter_ring,
                )
                specialized = special_old_base(
                    [specialization(value) for value in discriminant.list()]
                )
                denominator_valuation = valuation(
                    specialized, special_denominator
                )
                outside = specialized//special_denominator**denominator_valuation
                outside_squarefree = outside//outside.gcd(outside.derivative())
                infinity_valuation = 8-specialized.degree()
                support_degree = outside_squarefree.degree()
                if denominator_valuation != 8:
                    support_degree += 1
                if infinity_valuation:
                    support_degree += 1
                key = (
                    support_degree,
                    denominator_valuation-8,
                    infinity_valuation,
                )
                distribution[key] += 1
                if support_degree <= 8:
                    candidates.append(
                        (
                            int(a1_value), int(a0_value),
                            int(b2_value), int(b1_value),
                            support_degree, denominator_valuation-8,
                            infinity_valuation,
                        )
                    )

print(
    "Q80THIRDLOCAL|field=GF(7)|denominator=W+3|"
    f"tested={sum(distribution.values())}|"
    f"distribution={tuple(sorted(distribution.items()))}",
    flush=True,
)
print(
    f"Q80THIRDLOCAL|support_degree_le_8={tuple(candidates[:80])}|"
    f"count={len(candidates)}|status=PASS_BOUNDED_EXPERIMENT",
    flush=True,
)

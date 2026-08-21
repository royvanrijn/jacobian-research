#!/usr/bin/env sage
"""Search degree-one sections on the exact Picard-20 q=8 cubic model mod p.

The pointed genus-one equation is ``v^2 = 2*C_W(t)``, where ``C_W`` is the
cubic produced by ``derive_picard20_q8_s_chord.sage``.  This exhausts every
Mobius function t(W) over GF(p), with projective coefficient tuples normalized
once, and reports those for which the right side is a square in GF(p)(W).
"""

import argparse

from sage.all import *
from itertools import product as cartesian_product


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=23)
arguments = parser.parse_args()

field = GF(arguments.p)
K = FunctionField(field, "W")
W = K.gen()
RT = PolynomialRing(K, "t")
t = RT.gen()


def chord_quintic(z):
    return (
        t**5 + (field(21)/50*z-field(323)/200)*t**4
        + (-field(483)/625*z+field(129)/1250)*t**3
        + (field(1323)/62500*z**2+field(11907)/31250*z+field(1)/2)*t**2
        - field(31311)/781250*z**2*t + field(194481)/78125000*z**4
    )


mobius_a = field(3)/17
mobius_b = field(147)/425
completion_denominator = mobius_a*W-mobius_b
z = (t-W)/completion_denominator
completed = RT(completion_denominator**4*chord_quintic(z))
assert completed.valuation(t-field(49)/25) == 2
cubic = completed // (t-field(49)/25)**2
assert cubic.degree() == 3


def projective_tuples():
    for pivot in range(4):
        for tail in cartesian_product(field, repeat=int(3-pivot)):
            values = [field.zero()]*pivot + [field.one()] + list(tail)
            yield tuple(values)


hits = []
tested = 0
for coefficients in projective_tuples():
    numerator_constant, numerator_linear, denominator_constant, denominator_linear = coefficients
    determinant = (
        numerator_linear*denominator_constant
        - numerator_constant*denominator_linear
    )
    if not determinant:
        continue
    tested += 1
    t_value = K(
        (numerator_constant+numerator_linear*W)
        / (denominator_constant+denominator_linear*W)
    )
    value = K(2*cubic(t_value))
    if not value.is_square():
        continue
    square_root = value.sqrt()
    canonical_root = min(square_root, -square_root, key=str)
    hit = (coefficients, t_value, canonical_root)
    hits.append(hit)
    print(
        f"PICARD20Q8SECTION|p={arguments.p}|coefficients="
        + ",".join(map(str, map(int, coefficients)))
        + f"|t={t_value}|v={canonical_root}",
        flush=True,
    )

print(
    f"PICARD20Q8SECTION|p={arguments.p}|mobius_tested={tested}"
    f"|degree_one_hits={len(hits)}|status=PASS",
    flush=True,
)

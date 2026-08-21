#!/usr/bin/env sage
"""Exhaust degree-one S-chord pencils over GF(23), modulo base PGL2.

The ambient functions are ``a(t)+b(t)z`` with ``deg(a),deg(b)<=1``.  A pencil
is a two-plane in the four-dimensional coefficient space, so reduced row
echelon matrices enumerate it exactly once.  Each candidate is converted to
the corresponding hyperelliptic equation; only genus-one squarefree parts
are classified.  This is a bounded finite-field search, not a QQ existence
or nonexistence proof beyond the stated chart.
"""

import argparse

from sage.all import *
from itertools import product as cartesian_product


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=23)
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--limit", type=int, default=0)
parser.add_argument("--progress", type=int, default=10000)
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


def rref_planes():
    width = 4
    for first_pivot in range(width):
        for second_pivot in range(first_pivot+1, width):
            free_positions = []
            for row, pivot in enumerate((first_pivot, second_pivot)):
                for column in range(pivot+1, width):
                    if column not in (first_pivot, second_pivot):
                        free_positions.append((row, column))
            for values in cartesian_product(
                field, repeat=int(len(free_positions))
            ):
                rows = [[field.zero() for _ in range(width)] for _ in range(2)]
                rows[0][first_pivot] = 1
                rows[1][second_pivot] = 1
                for (row, column), value in zip(free_positions, values):
                    rows[row][column] = value
                yield tuple(tuple(row) for row in rows)


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        rank = max(0, ord_delta-1)
        return rank, None if rank == 0 else f"A{rank}"
    if ord_delta == 2:
        return 0, None
    if ord_delta == 3:
        return 1, "A1"
    if ord_delta == 4:
        return 2, "A2"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, "D4"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        rank = ord_delta-2
        return rank, f"D{rank}"
    if ord_delta == 8:
        return 6, "E6"
    if ord_delta == 9:
        return 7, "E7"
    if ord_delta == 10:
        return 8, "E8"
    return -100, f"UNKNOWN({ord_a},{ord_b},{ord_delta})"


def classify(quartic):
    coefficients = [quartic[index] for index in range(5)]
    e, d, c, b, a = coefficients
    invariant_i = 12*a*e-3*b*d+c**2
    invariant_j = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    if 4*invariant_i**3-invariant_j**2 == 0:
        return None
    curve = EllipticCurve(K, [0, 0, 0, -27*invariant_i, -27*invariant_j])
    discriminant = curve.discriminant()
    ade = []
    root_rank = 0
    for factor, delta_order in discriminant.numerator().factor():
        ord_a = curve.a4().numerator().valuation(factor)
        ord_b = curve.a6().numerator().valuation(factor)
        rank, component = kodaira_data(ord_a, ord_b, delta_order)
        if rank < 0:
            return None
        root_rank += factor.degree()*rank
        if component is not None:
            ade.extend([component]*factor.degree())
    infinity = (
        8-curve.a4().numerator().degree(),
        12-curve.a6().numerator().degree(),
        24-discriminant.numerator().degree(),
    )
    infinity_rank, infinity_component = kodaira_data(*infinity)
    if infinity_rank < 0:
        return None
    root_rank += infinity_rank
    if infinity_component is not None:
        ade.append(infinity_component)
    return tuple(sorted(ade)), root_rank, curve, discriminant, infinity


target_ade = tuple(sorted(("A7", "A4", "A3", "A2")))
tested = 0
selected = 0
gcd_survivors = 0
genus_one_hits = 0
target_hits = 0
for index, rows in enumerate(rref_planes()):
    if index < arguments.start:
        continue
    if arguments.limit and selected >= arguments.limit:
        break
    selected += 1
    tested = index+1
    left, right = rows
    left_constant = RT(left[0]+left[1]*t)
    left_z = RT(left[2]+left[3]*t)
    right_constant = RT(right[0]+right[1]*t)
    right_z = RT(right[2]+right[3]*t)
    denominator = W*left_z-right_z
    numerator = right_constant-W*left_constant
    if not denominator:
        continue
    completed = RT(denominator**4*chord_quintic(numerator/denominator))
    if completed.degree() <= 4:
        possible = True
    else:
        common = gcd(completed, completed.derivative())
        possible = common.degree() >= (completed.degree()-4+1)//2
    if not possible:
        continue
    gcd_survivors += 1
    factors = completed.factor()
    odd_part = factors.unit()*prod(
        factor for factor, exponent in factors if exponent % 2
    )
    if odd_part.degree() > 4:
        continue
    while odd_part.degree() < 4:
        # Homogenizing a cubic supplies the fourth branch point at infinity;
        # the invariant formulas accept a zero leading coefficient directly.
        break
    result = classify(odd_part)
    if result is None:
        continue
    genus_one_hits += 1
    ade, root_rank, curve, discriminant, infinity = result
    if ade == target_ade:
        target_hits += 1
        print(
            f"PICARD20Q8GRASS|target={target_hits}|index={index}|rows={rows}"
            f"|completed_degree={completed.degree()}"
            f"|factors={tuple((factor.degree(), exponent) for factor, exponent in factors)}"
            f"|quartic={odd_part}|ADE={'+'.join(ade)}"
            f"|A={curve.a4()}|B={curve.a6()}"
            f"|disc={discriminant.factor()}|infinity={infinity}",
            flush=True,
        )
    elif genus_one_hits <= 20:
        print(
            f"PICARD20Q8GRASS|hit={genus_one_hits}|index={index}|rows={rows}"
            f"|ADE={'+'.join(ade)}|root_rank={root_rank}",
            flush=True,
        )
    if arguments.progress and selected % arguments.progress == 0:
        print(
            f"PICARD20Q8GRASS|progress={selected}|global_index={index}"
            f"|gcd_survivors={gcd_survivors}|genus_one_hits={genus_one_hits}"
            f"|target_hits={target_hits}",
            flush=True,
        )

print(
    f"PICARD20Q8GRASS|p={arguments.p}|start={arguments.start}"
    f"|selected={selected}|global_tested={tested}|gcd_survivors={gcd_survivors}"
    f"|genus_one_hits={genus_one_hits}|target_hits={target_hits}|status=PASS",
    flush=True,
)

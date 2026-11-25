#!/usr/bin/env python3
"""Exact terminal-thirteen fixed-ray certificates for high and private curvature."""

from fractions import Fraction as Q
from math import comb, factorial

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_seven_private_rows import verify_annulus
from verify_multivariate_terminal_six_joint_ray_bands import radial_replay
from verify_multivariate_terminal_six_private_angles import exponential_lower, rational_cap, require


ORDERS = (12, 13)
RAY_DENOMINATOR = 1000000
# Increasing gaps, decreasing gaps, shared loss, ray numerators, radial cells.
CASES = (
    ((), (), Q('63/8'),
     (966020, 891043, 713740, 551444, 284402, 70426, -210314, -446114, -657041, -825038, -938163), 140),
    ((1,), (), Q('6'),
     (966889, 891043, 718210, 552150, 290766, 70426, -204798, -436898, -654030, -824012, -936776), 116),
    ((1, 1, 1, 1), (), Q('0'),
     (967190, 910120, 718780, 578670, 291420, 82110, -204250, -429060, -655710, -821580, -937830), 116),
    ((1, 1), (1, 1), Q('0'),
     (968120, 866950, 726020, 553770, 303190, 27530, -192970, -438510, -641630, -823330, -930570), 56),
    ((1, 1, 1), (1,), Q('0'),
     (967920, 866950, 723900, 553770, 300250, 50940, -195800, -416400, -645880, -825080, -932220), 72),
    ((3,), (), Q('3/8'),
     (966009, 891043, 712130, 566871, 283243, 93306, -210314, -443661, -659225, -820574, -946177), 156),
    ((1, 1, 1), (), Q('9/8'),
     (966140, 891680, 713620, 564510, 284530, 89880, -208470, -444780, -658500, -823330, -943200), 146),
    ((1, 2), (), Q('3/4'),
     (966090, 891680, 713180, 564510, 283540, 93770, -209880, -444780, -658500, -821580, -944690), 154),
    ((1, 1), (), Q('27/8'),
     (966400, 891680, 715770, 568060, 287480, 86000, -207070, -438510, -657100, -823330, -940160), 138),
    ((2,), (), Q('23/8'),
     (966140, 891680, 714050, 568060, 285510, 89880, -208470, -441640, -657100, -823330, -942450), 140),
    ((1,), (1,), Q('13/4'),
     (966720, 891680, 718780, 539320, 293380, 58740, -202840, -438510, -650100, -821580, -933850), 120),
    ((1, 1), (1,), Q('7/8'),
     (966770, 891680, 719210, 539320, 294360, 62630, -200030, -444780, -650100, -825080, -934650), 126),
    ((2,), (1,), Q('3/8'),
     (966560, 891680, 718780, 550180, 294360, 66530, -200030, -444780, -650100, -825080, -935450), 132),
)



def coefficient_box(order, increasing, decreasing, strength):
    caps = [Q(1)]
    for k in range(1, order):
        power = Q(comb(order, k)**order)
        for index, gaps in ((k, increasing), (order-k, decreasing)):
            for gap in gaps:
                power *= Q(factorial(gap*index)**order, factorial(gap*order)**index)
        quarter_argument = strength*k*(order-k)/8
        power /= exponential_lower(quarter_argument)**4
        caps.append(rational_cap(power, order))
    caps.append(Q(1))
    return tuple((cap, cap) if k in (0, order) else (Q(0), cap)
                 for k, cap in enumerate(caps))


def verify_case(increasing, decreasing, strength, numerators, expected_cells):
    cosines = tuple(Q(value, RAY_DENOMINATOR) for value in numerators)
    knots = (Q(1), *cosines, Q(-1))
    require(len(cosines) == ORDERS[-1] - 2, "ray count")
    require(all(left > right for left, right in zip(knots, knots[1:])), "strict ray order")
    radial_cells = anchor_checks = 0
    for order in ORDERS:
        box = coefficient_box(order, increasing, decreasing, strength)
        verify_annulus(box, order)
        for boundary, cosine in enumerate(knots):
            root_count = (boundary + (order == ORDERS[-1])) // 2
            require((-1)**root_count * anchor_value(order, cosine) > 0,
                    ("binomial anchor", order, boundary, cosine))
            anchor_checks += 1
        for cosine in cosines:
            cells = radial_replay(order, box, (cosine, cosine))
            require(cells is not None,
                    ("unresolved ray", increasing, decreasing, strength, order, cosine))
            radial_cells += cells
    require(radial_cells == expected_cells,
            ("radial cell count", increasing, decreasing, radial_cells))
    print(f"{increasing} | {decreasing}, shared loss {strength}: "
          f"{radial_cells} radial cells and {anchor_checks} anchor signs.")
    return radial_cells, anchor_checks


def main():
    radial_cells = anchor_checks = 0
    for case in CASES:
        cells, signs = verify_case(*case)
        radial_cells += cells
        anchor_checks += signs
    print(f"Verified {len(CASES)} cases, {radial_cells} radial cells, "
          f"and {anchor_checks} anchor signs.")


if __name__ == "__main__":
    main()

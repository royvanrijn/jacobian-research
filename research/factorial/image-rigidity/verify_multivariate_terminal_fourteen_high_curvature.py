#!/usr/bin/env python3
"""Exact terminal-fourteen fixed rays for high and private curvature."""

from fractions import Fraction as Q

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_seven_private_rows import verify_annulus
from verify_multivariate_terminal_six_joint_ray_bands import radial_replay
from verify_multivariate_terminal_six_private_angles import require
from verify_multivariate_terminal_thirteen_high_curvature import coefficient_box


ORDERS = (13, 14)
RAY_DENOMINATOR = 1000000
# Increasing gaps, decreasing gaps, shared loss, ray numerators, radial cells.
CASES = (
    ((), (), Q('8'),
     (971097, 909897, 753639, 606991, 374696, 155689, -80717, -300591, -519183, -688954, -848270, -945124), 148),
    ((1,), (), Q('23/4'),
     (971379, 909897, 755762, 605017, 377694, 191386, -76953, -300591, -515860, -688954, -847641, -944594), 140),
    ((1, 1, 1, 1), (), Q('0'),
     (971889, 923880, 757455, 638056, 380689, 193609, -73188, -292761, -515860, -690282, -848270, -944594), 126),
    ((1, 1), (1, 1), Q('0'),
     (972696, 923880, 763759, 641869, 391636, 155689, -63393, -284912, -505006, -692932, -838063, -943526), 66),
    ((1, 1, 1), (1,), Q('0'),
     (972521, 923880, 762085, 641869, 389153, 155689, -64901, -290801, -507518, -683619, -840000, -942989), 78),
    ((3,), (), Q('1/2'),
     (971045, 909897, 752361, 628454, 373195, 193609, -79965, -288839, -519183, -696888, -854495, -948764), 172),
    ((1, 1, 1), (), Q('9/8'),
     (971097, 909897, 753213, 630382, 374696, 193609, -79212, -294721, -518353, -696888, -852640, -948251), 166),
    ((1, 2), (), Q('3/4'),
     (970993, 909897, 752148, 628454, 373195, 193609, -79965, -292761, -519183, -698202, -853878, -949274), 182),
    ((1, 1), (), Q('13/4'),
     (971097, 909897, 753852, 632306, 375196, 191386, -78459, -296679, -517522, -692932, -850775, -946699), 148),
    ((2,), (), Q('23/8'),
     (971045, 909897, 753213, 630382, 375196, 193609, -78459, -294721, -518353, -694253, -852020, -946699), 156),
    ((1,), (1,), Q('3'),
     (971277, 909897, 756398, 606991, 380689, 155689, -73188, -296679, -513362, -690282, -843842, -945124), 152),
    ((1, 1), (1,), Q('5/8'),
     (971379, 909897, 757033, 605017, 382185, 155689, -71682, -300591, -511695, -687624, -843842, -945651), 144),
    ((2,), (1,), Q('1/8'),
     (971199, 909897, 756610, 605017, 381686, 157926, -70928, -302546, -511695, -687624, -844478, -945651), 150),
)


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

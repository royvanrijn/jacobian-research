#!/usr/bin/env python3
"""Exact terminal-eleven ray certificates with shared or private loss."""

from fractions import Fraction as Q

from verify_multivariate_terminal_nine_high_curvature import (
    anchor_value, coefficient_bounds as shared_coefficient_bounds,
)
from verify_multivariate_terminal_seven_private_rows import verify_annulus
from verify_multivariate_terminal_six_joint_ray_bands import radial_replay
from verify_multivariate_terminal_six_private_angles import coefficient_bounds, require


ORDERS = (10, 11)
RAY_DENOMINATOR = 100000
# Increasing gaps, decreasing gaps, shared loss, ray numerators, radial cells.
CASES = (
    ((), (), Q(8),
     (95195, 86554, 60463, 37608, 4516, -24236, -52163, -74496, -91132), 98),
    ((1,), (), Q(6),
     (95284, 86554, 60811, 38066, 5243, -22962, -51814, -73733, -91024), 102),
    ((1, 1, 1, 1), (), Q(0),
     (95311, 86554, 60880, 38066, 5243, -22536, -51814, -73733, -91132), 88),
    ((1, 1), (1, 1), Q(0),
     (95451, 86554, 61778, 34372, 6552, -20401, -50764, -73989, -90248), 62),
    ((1, 1, 1), (1,), Q(0),
     (95425, 86554, 61502, 38066, 6262, -20829, -50939, -73989, -90473), 68),
    ((3,), (), Q(1, 4),
     (95124, 85710, 59974, 40345, 4370, -21683, -52337, -73989, -92078), 118),
    ((1, 1, 1), (), Q(9, 8),
     (95160, 85994, 60184, 37608, 4516, -22962, -52163, -74243, -91769), 114),
    ((1, 2), (), Q(5, 8),
     (95124, 83038, 60044, 37608, 4370, -22536, -52337, -73989, -91976), 116),
    ((1, 1), (), Q(27, 8),
     (95195, 86275, 60463, 38066, 4807, -23387, -51989, -73989, -91347), 112),
    ((2,), (), Q(11, 4),
     (95133, 85994, 60184, 37608, 4516, -22536, -52163, -73989, -91769), 122),
    ((1,), (1,), Q(13, 4),
     (95213, 86554, 60811, 37608, 5243, -24236, -51640, -74748, -90805), 102),
    ((1, 1), (1,), Q(7, 8),
     (95204, 83038, 60880, 37148, 5534, -24236, -51465, -74496, -90805), 110),
    ((2,), (1,), Q(3, 8),
     (95195, 83038, 60811, 37148, 5534, -24660, -51290, -74243, -90915), 114),
)


def coefficient_box(order, increasing, decreasing, strength):
    if increasing or decreasing:
        caps = coefficient_bounds(order, increasing, strength, decreasing)
    else:
        caps = shared_coefficient_bounds(order, strength)
    return [(cap, cap) if k in (0, order) else (Q(0), cap)
            for k, cap in enumerate(caps)]


def verify_case(increasing, decreasing, strength, numerators, expected_cells):
    cosines = tuple(Q(k, RAY_DENOMINATOR) for k in numerators)
    knots = (Q(1), *cosines, Q(-1))
    require(len(cosines) == 9, "ray count")
    require(all(a > b for a, b in zip(knots, knots[1:])), "strict ray order")
    radial_cells = anchor_checks = 0
    for order in ORDERS:
        box = coefficient_box(order, increasing, decreasing, strength)
        verify_annulus(box, order)
        for index, cosine in enumerate(knots):
            root_count = (index + (order == 11)) // 2
            require((-1)**root_count * anchor_value(order, cosine) > 0,
                    ("binomial anchor", order, index, cosine))
            anchor_checks += 1
        for cosine in cosines:
            cells = radial_replay(order, box, (cosine, cosine))
            require(cells is not None, ("unresolved ray", increasing, decreasing, strength, order, cosine))
            radial_cells += cells
    require(radial_cells == expected_cells, ("radial cell count", increasing, decreasing, radial_cells))
    print(f"{increasing} | {decreasing}, shared loss {strength}: "
          f"{radial_cells} radial cells and {anchor_checks} anchor signs.")
    return radial_cells, anchor_checks


def main():
    radial_cells = anchor_checks = 0
    for case in CASES:
        cells, signs = verify_case(*case)
        radial_cells += cells
        anchor_checks += signs
    print(f"Verified {len(CASES)} cases, {radial_cells} radial cells, and {anchor_checks} anchor signs.")


if __name__ == "__main__":
    main()

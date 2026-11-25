#!/usr/bin/env python3
"""Exact terminal-twelve fixed-ray certificates for high and private curvature."""

from fractions import Fraction as Q

from verify_multivariate_terminal_eleven_high_curvature import coefficient_box
from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_seven_private_rows import verify_annulus
from verify_multivariate_terminal_six_joint_ray_bands import radial_replay
from verify_multivariate_terminal_six_private_angles import require


ORDERS = (11, 12)
RAY_DENOMINATOR = 100000
# Increasing gaps, decreasing gaps, shared loss, ray numerators, radial cells.
CASES = (
    ((), (), Q(8),
     (96004, 86517, 66635, 45239, 17827, -7594, -35742, -58814, -79217, -92201), 114),
    ((1,), (), Q(6),
     (96078, 86517, 66960, 45239, 18305, -5946, -35265, -58814, -78949, -92078), 114),
    ((1, 1, 1, 1), (), Q(0),
     (96105, 86517, 67014, 45239, 18424, -5122, -35265, -58223, -79217, -92201), 90),
    ((1, 1), (1, 1), Q(0),
     (96225, 86517, 67875, 45239, 19735, -4297, -34149, -59108, -77864, -91830), 50),
    ((1, 1, 1), (1,), Q(0),
     (96199, 86517, 67660, 45239, 19378, -5946, -34468, -56733, -78274, -91830), 58),
    ((3,), (), Q(3, 8),
     (95970, 86517, 66363, 44414, 17468, -5122, -35900, -58814, -80010, -92801), 124),
    ((1, 1, 1), (), Q(9, 8),
     (95983, 86517, 66417, 44827, 17708, -5534, -35742, -59694, -79747, -92683), 126),
    ((1, 2), (), Q(3, 4),
     (95977, 86517, 66417, 44414, 17588, -5534, -35742, -59402, -79879, -92801), 126),
    ((1, 1), (), Q(27, 8),
     (96017, 86517, 66635, 44827, 17946, -5946, -35583, -59402, -79483, -92444), 120),
    ((2,), (), Q(23, 8),
     (95983, 86517, 66526, 44827, 17827, -5122, -35583, -59402, -79615, -92564), 122),
    ((1,), (1,), Q(13, 4),
     (96038, 86517, 67014, 45239, 18543, -8417, -34946, -58519, -78545, -92078), 100),
    ((1, 1), (1,), Q(7, 8),
     (96045, 86517, 67068, 45239, 18663, -8006, -34946, -58814, -78680, -92201), 108),
    ((2,), (1,), Q(3, 8),
     (96024, 86517, 67014, 45239, 18663, -7594, -34946, -58814, -78680, -92201), 120),
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

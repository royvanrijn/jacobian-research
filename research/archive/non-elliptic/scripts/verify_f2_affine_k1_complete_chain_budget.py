#!/usr/bin/env python3
"""Verify the complete-chain point budget for the F2 k=1 E8 packet."""

from __future__ import annotations

from fractions import Fraction


# (degree, action multiplicity, fixed sheets, q_1, q_2, q_4)
SIMPLE_INERTIA_ROWS = (
    (6, 1, 2, 2, 0, 0),
    (10, 1, 2, 4, 0, 0),
    (12, 1, 4, 0, 2, 0),
    (15, 1, 3, 6, 0, 0),
    (20, 1, 4, 0, 4, 0),
    (24, 2, 4, 0, 1, 2),
    (30, 1, 4, 1, 6, 0),
    (30, 1, 2, 14, 0, 0),
    (40, 2, 4, 0, 1, 4),
    (60, 1, 4, 0, 14, 0),
    (120, 1, 4, 0, 1, 14),
)


def global_budget(case: str, degree: int, smooth_blowups: int) -> Fraction:
    """Return ch_2 of the global logarithmic boundary complex."""

    source_square = {
        "squarefree": -6 - smooth_blowups,
        "double": -11 - smooth_blowups,
    }[case]
    target_square = -5
    return Fraction(
        source_square - degree * target_square + 2 * (degree - 1), 2
    )


def determinant_square(
    case: str,
    degree: int,
    ramified_sheets: int,
    contact: int,
    smooth_blowups: int,
) -> int:
    """Recover R_log^2 from L_X=f^*L_Y+R_log."""

    source_square = {
        "squarefree": -6 - smooth_blowups,
        "double": -11 - smooth_blowups,
    }[case]
    target_square = -5
    target_pairing = ramified_sheets * (contact - 8)
    return source_square - degree * target_square - 2 * target_pairing


def complete_chain_point_budget(
    case: str,
    degree: int,
    fixed_sheets: int,
    contact: int,
    smooth_blowups: int,
) -> Fraction:
    """Subtract the complete Cartier cycle and its conormal kernel degree."""

    ramified_sheets = degree - fixed_sheets
    square = determinant_square(
        case, degree, ramified_sheets, contact, smooth_blowups
    )
    kernel_degree = ramified_sheets * (contact - 7)
    return global_budget(case, degree, smooth_blowups) - (
        kernel_degree + Fraction(square, 2)
    )


def abstract_point_budget(
    source_square: int,
    target_square: int,
    degree: int,
    fixed_sheets: int,
    target_intersection: int,
) -> Fraction:
    """Evaluate the degree-independent complete-chain identity."""

    moved_sheets = degree - fixed_sheets
    global_ch2 = Fraction(
        source_square - degree * target_square + 2 * (degree - 1), 2
    )
    determinant_square_value = (
        source_square
        - degree * target_square
        - 2 * moved_sheets * target_intersection
    )
    kernel_degree = moved_sheets * (target_intersection + 1)
    return global_ch2 - (
        kernel_degree + Fraction(determinant_square_value, 2)
    )


def main() -> None:
    for source_square in range(-20, 8):
        for target_square in range(-12, 8):
            for degree in range(2, 18):
                for fixed in range(1, degree):
                    for target_intersection in range(-10, 5):
                        assert abstract_point_budget(
                            source_square,
                            target_square,
                            degree,
                            fixed,
                            target_intersection,
                        ) == fixed - 1

    # The cancellation is independent of the carrier contact, common-model
    # smooth blowups, and squarefree/double terminal row.
    for case in ("squarefree", "double"):
        for degree, _, fixed, q1, q2, q4 in SIMPLE_INERTIA_ROWS:
            residue_weight = q1 + 2 * q2 + 4 * q4
            assert degree == fixed + 2 * residue_weight
            for contact in range(9):
                for smooth_blowups in range(8):
                    point_budget = complete_chain_point_budget(
                        case, degree, fixed, contact, smooth_blowups
                    )
                    assert point_budget == fixed - 1

    expected_deficits = (
        (6, 1, 2, 2, 3),
        (10, 1, 2, 4, 7),
        (12, 1, 4, 4, 5),
        (15, 1, 3, 6, 10),
        (20, 1, 4, 8, 13),
        (24, 2, 4, 10, 17),
        (30, 1, 4, 13, 23),
        (30, 1, 2, 14, 27),
        (40, 2, 4, 18, 33),
        (60, 1, 4, 28, 53),
        (120, 1, 4, 58, 113),
    )
    actual_deficits = []
    for degree, multiplicity, fixed, q1, q2, q4 in SIMPLE_INERTIA_ROWS:
        residue_weight = q1 + 2 * q2 + 4 * q4
        cusp_lower = 2 * residue_weight
        deficit = cusp_lower - (fixed - 1)
        assert deficit == degree - 2 * fixed + 1
        assert deficit > 0
        actual_deficits.append(
            (degree, multiplicity, fixed, residue_weight, deficit)
        )
    assert tuple(actual_deficits) == expected_deficits

    print(
        "PASS: after the complete Cartier determinant cycle and conormal "
        "kernel are subtracted, every squarefree/double F2 k=1 one-curve "
        "model has point budget u-1; all 13 simple-inertia E8 actions have "
        "cusp lower 2R>u-1 and require a negative noncyclic correction of "
        "exact minimum 3,7,5,10,13,17,23,27,33,53,113 by atlas row"
    )


if __name__ == "__main__":
    main()

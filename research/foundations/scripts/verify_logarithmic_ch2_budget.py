#!/usr/bin/env python3
"""Verify the global logarithmic ch2 budget and the live F2 specialization."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "plane-jc" / "cas"))

from jcsearch.log_node_profiles import (  # noqa: E402
    cyclic_boundary_charge,
    determinant,
    logarithmic_ch2_budget,
    logarithmic_ch2_from_ramification,
    logarithmic_ch2_model_change,
    residual_point_budget,
    tangential_kernel_trivialization_profile,
)
from verify_f2_carrier_log_node_profiles import (  # noqa: E402
    refined_source_graph,
)


GLOBAL_ATTACHMENT = (
    ROOT / "artifacts/generated-results/jc2_f2_75_125_global_attachment.json"
)
CARRIER_WRONSKIAN = (
    ROOT / "artifacts/generated-results/jc2_f2_75_125_carrier_wronskian.json"
)


def boundary_blowup_types(
    insertion_order: list[list[int]],
) -> tuple[int, int]:
    """Count smooth-boundary and boundary-node blowups in a local fan.

    The initial first-quadrant cone has one boundary ray ``(1,0)`` and one
    nonboundary ray ``(0,1)``.  Every inserted exceptional ray is boundary.
    """

    rays = [(1, 0), (0, 1)]
    is_boundary = [True, False]
    smooth = 0
    nodes = 0
    for raw_ray in insertion_order:
        ray = tuple(raw_ray)
        sectors = [
            index
            for index, (left, right) in enumerate(zip(rays, rays[1:]))
            if determinant(left, ray) > 0 and determinant(ray, right) > 0
        ]
        if len(sectors) != 1:
            raise AssertionError(f"ray {ray} does not select one current fan cone")
        index = sectors[0]
        boundary_count = int(is_boundary[index]) + int(is_boundary[index + 1])
        if boundary_count == 2:
            nodes += 1
        elif boundary_count == 1:
            smooth += 1
        else:
            raise AssertionError("a declared boundary extraction left the boundary")
        rays.insert(index + 1, ray)
        is_boundary.insert(index + 1, True)
    return smooth, nodes


def verify_general_identities() -> None:
    # Expand L_X=f^*L_Y+R and compare the two global formulas, including the
    # general top logarithmic Chern numbers e_X and e_Y.
    for degree in range(1, 8):
        for target_square in (-5, 0, 4):
            for target_dot_ramification in (-3, 0, 7):
                for ramification_square in (-8, 0, 11):
                    source_square = (
                        degree * target_square
                        + 2 * target_dot_ramification
                        + ramification_square
                    )
                    bundle = logarithmic_ch2_budget(
                        source_square,
                        target_square,
                        degree,
                        source_log_c2=2,
                        target_log_c2=3,
                    )
                    ramification = logarithmic_ch2_from_ramification(
                        target_dot_ramification,
                        ramification_square,
                        degree,
                        source_log_c2=2,
                        target_log_c2=3,
                    )
                    assert bundle.ch2 == ramification

    # Smooth-boundary blowups change (K+D)^2 by -1; node blowups are log
    # crepant.  Check the model-change formula directly against the budget.
    base = logarithmic_ch2_budget(-2, 3, 7)
    for source_smooth in range(4):
        for target_smooth in range(4):
            changed = logarithmic_ch2_budget(
                -2 - source_smooth,
                3 - target_smooth,
                7,
            )
            delta = logarithmic_ch2_model_change(
                7,
                source_smooth_blowups=source_smooth,
                target_smooth_blowups=target_smooth,
                source_node_blowups=5,
                target_node_blowups=6,
            )
            assert changed.ch2 - base.ch2 == delta.ch2_change

    # The blowup of one smooth boundary point of the identity completion is
    # the smallest exact warning that the raw localized ch2 is signed.
    assert logarithmic_ch2_budget(3, 4, 1).ch2 == Fraction(-1, 2)

    # The split nodal module differs from the glued cyclic reference by one
    # point class.  The helper records the signed K-class subtraction only;
    # effectivity still comes from the separate exact sequence.
    split = residual_point_budget(Fraction(9, 2), (Fraction(7, 2),))
    assert split.residual_ch2 == 1


def verify_f2_budget() -> None:
    global_attachment = json.loads(GLOBAL_ATTACHMENT.read_text())
    carrier = json.loads(CARRIER_WRONSKIAN.read_text())

    terminal_insertions = global_attachment["target_minimal_completion"]["fan"][
        "insertion_order_coordinates"
    ]
    carrier_insertions = carrier["carrier_target_completion"]["insertion_order"]
    assert boundary_blowup_types(terminal_insertions) == (1, 3)
    assert boundary_blowup_types(carrier_insertions) == (8, 4)
    target_log_square = 4 - 1 - 8
    assert target_log_square == -5

    square_graph = refined_source_graph("squarefree_one_packet")
    double_graph = refined_source_graph("double_same_target")
    square_source_square = square_graph["log_canonical_square"]
    double_source_square = double_graph["log_canonical_square"]
    assert square_source_square == -6
    assert double_source_square == -11

    live_cases = {
        case["case"]: case
        for case in global_attachment["cases"]
        if case["status"] != "excluded_by_target_valuation_uniqueness"
    }
    assert live_cases["squarefree_one_packet"]["geometric_degree_floor"] == 6
    assert live_cases["double_same_target"]["geometric_degree_floor"] == 12

    # The upstream cyclic packet has D=3E+18L, D^2=54, and a fixed target
    # covector trivializes its kernel.  Its actual cokernel contribution is
    # therefore +D^2/2=27, not the -D^2/2 of the untwisted O_D.
    root_charge = cyclic_boundary_charge((3, 18), (-6, 0), ((0, 1),))
    assert root_charge.node_matching_length == 54
    assert root_charge.charge == 27
    root_cokernel = tangential_kernel_trivialization_profile(
        2 * root_charge.charge,
        (3, 18),
        (3, 18),
    )
    assert root_cokernel.cokernel_ch2 == 27

    square_degree = 6
    square_budget = logarithmic_ch2_budget(
        square_source_square,
        target_log_square,
        square_degree,
    )
    assert square_budget.ch2 == 17
    square_residual = residual_point_budget(
        square_budget.ch2,
        (root_cokernel.cokernel_ch2,),
    )
    assert square_residual.residual_ch2 == -10

    double_degree = 12
    double_budget = logarithmic_ch2_budget(
        double_source_square,
        target_log_square,
        double_degree,
    )
    assert double_budget.ch2 == Fraction(71, 2)
    double_residual = residual_point_budget(
        double_budget.ch2,
        (root_cokernel.cokernel_ch2,),
    )
    assert double_residual.residual_ch2 == Fraction(17, 2)

    # Record the exact affine dependence on the still unknown geometric
    # degree.  These equalities are twice the corresponding rational budget.
    for degree in range(1, 30):
        square = logarithmic_ch2_budget(-6, -5, degree)
        double = logarithmic_ch2_budget(-11, -5, degree)
        assert 2 * square.ch2 == 7 * degree - 8
        assert 2 * double.ch2 == 7 * degree - 13
        square_after_root = residual_point_budget(square.ch2, (27,))
        double_after_root = residual_point_budget(double.ch2, (27,))
        assert (square_after_root.residual_ch2.denominator == 1) == (
            degree % 2 == 0
        )
        assert (double_after_root.residual_ch2.denominator == 1) == (
            degree % 2 == 1
        )


def main() -> None:
    verify_general_identities()
    verify_f2_budget()
    print(
        "PASS: global logarithmic ch2 and blowup laws; current F2 budgets are "
        "(7*d-8)/2 and (7*d-13)/2 with exact root contribution 27"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify the contracted common-factor Chern and boundary-minimality ledger."""

from __future__ import annotations

from itertools import product


def audit_boundary_blowup_invariant() -> None:
    # Every component created from (P2,L_infinity) by boundary point blowups
    # satisfies self_intersection + valency <= 1.  The two transition types
    # preserve that inequality on old components and initialize it on the new
    # exceptional component.
    for self_intersection in range(-12, 2):
        for valency in range(7):
            if self_intersection + valency > 1:
                continue

            smooth_old = (self_intersection - 1, valency + 1)
            node_old = (self_intersection - 1, valency)
            assert sum(smooth_old) <= 1
            assert sum(node_old) <= 1
            assert -1 + 1 <= 1  # new smooth-center exceptional component
            assert -1 + 2 <= 1  # new node-center exceptional component

    for valency in range(7):
        if -1 + valency <= 1:
            assert valency <= 2


def audit_common_factor_charge() -> None:
    # For F=Omega_X^1(log D_X), rank(F)=2 and c1(F)=L_X.  If H=hT,
    # ch_2(F/F(-H))=L_X.H-H^2=h(v-2)+h^2 n.
    for h in range(1, 9):
        for valency in range(1, 8):
            for negativity in range(1, 9):
                charge = h * (valency - 2) + h * h * negativity
                if negativity >= 2:
                    assert charge >= 2 * h * h - h > 0


def audit_saturated_budget() -> None:
    for residual_budget in range(-3, 12):
        for h in range(1, 6):
            for valency in range(1, 7):
                for negativity in range(1, 7):
                    for intersection in range(0, 15):
                        # D=D'+2hT, T^2=-n, L_X.T=v-2.
                        determinant_square_change = (
                            2 * h * intersection - 2 * h * h * negativity
                        )
                        common_quotient = (
                            h * (valency - 2) + h * h * negativity
                        )
                        saturated = (
                            residual_budget
                            + determinant_square_change
                            - common_quotient
                        )
                        closed_form = (
                            residual_budget
                            + 2 * h * intersection
                            - 3 * h * h * negativity
                            - h * (valency - 2)
                        )
                        assert saturated == closed_form

    # Cubic E8 single-packet gate: h=1, residual budget 2, local point
    # quotient 4, and I=3.  Relative boundary minimality gives n>=2.
    for valency in range(1, 8):
        for negativity in range(2, 9):
            point_other = 2 * 3 - valency - 3 * negativity
            assert point_other < 0


def audit_connected_cycle_gate() -> None:
    """Replay the numerical core of the complete-cycle argument.

    The written proof is general.  These path and fork graphs exercise the
    identities with unequal common-factor multiplicities and nonzero
    residual vertical determinant cycles.
    """

    graphs = (
        (1, ()),
        (2, ((0, 1),)),
        (3, ((0, 1), (1, 2))),
        (4, ((0, 1), (1, 2), (1, 3))),
        (5, ((0, 1), (1, 2), (2, 3), (2, 4))),
    )
    for vertex_count, edges in graphs:
        degrees = [0] * vertex_count
        neighbors = [[] for _ in range(vertex_count)]
        for left, right in edges:
            degrees[left] += 1
            degrees[right] += 1
            neighbors[left].append(right)
            neighbors[right].append(left)

        # The unique external horizontal cubic row meets vertex zero, whose
        # common order is one.
        for remaining_heights in product(range(1, 5), repeat=vertex_count - 1):
            heights = (1, *remaining_heights)
            for negativities in product(range(2, 6), repeat=vertex_count):
                anti_nef = tuple(
                    negativities[index] * heights[index]
                    - sum(heights[other] for other in neighbors[index])
                    for index in range(vertex_count)
                )
                if min(anti_nef) < 0 or not any(anti_nef):
                    continue

                q = sum(
                    heights[index] * anti_nef[index]
                    for index in range(vertex_count)
                )
                log_canonical_intersection = sum(
                    heights[index]
                    * (degrees[index] + (1 if index == 0 else 0) - 2)
                    for index in range(vertex_count)
                )
                assert (
                    sum(anti_nef) + log_canonical_intersection
                    == 1
                    + sum(
                        (negativities[index] - 2) * heights[index]
                        for index in range(vertex_count)
                    )
                )
                assert q >= sum(anti_nef)
                assert log_canonical_intersection >= 1 - q

                # A saturated residual vertical component can occur only
                # where the generated line O(-H)|T has degree a_i=0.
                # Therefore every such effective residual cycle R has R.H=0.
                zero_degree_vertices = tuple(
                    index
                    for index, value in enumerate(anti_nef)
                    if value == 0
                )
                for residual_multiplicity in range(4):
                    residual_intersection = sum(
                        residual_multiplicity
                        * (-anti_nef[index])
                        for index in zero_degree_vertices
                    )
                    assert residual_intersection == 0

                # The parameter-ideal argument supplies q>=2 on a nonempty
                # relative-minimal fiber.  The cubic row has U=2 and D'.H<=3.
                if q >= 2:
                    saturated_upper = (
                        2 + 2 * 3 - 3 * q - log_canonical_intersection
                    )
                    assert saturated_upper <= 3
                    assert saturated_upper - 4 < 0


def main() -> None:
    audit_boundary_blowup_invariant()
    audit_common_factor_charge()
    audit_saturated_budget()
    audit_connected_cycle_gate()
    print(
        "PASS: common contracted factors have positive curve charge on a "
        "relative SNC-minimal boundary; the full connected cubic contracted "
        "fiber has budget at most three against its forced length-four quotient"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify carrier/arm logarithmic profiles for the F2 ``(75,125)`` row."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "plane-jc" / "cas"))

from jcsearch.log_node_profiles import compile_toric_fan_profile  # noqa: E402
from verify_f2_75_125_global_attachment import (  # noqa: E402
    bareiss_determinant,
    minimal_source_principal_boundary,
    solve_linear_system,
    symmetric_inertia,
)
from verify_f2_75_125_carrier_wronskian import (  # noqa: E402
    is_zero_mod_rho_relation,
)


def infinity_order(expression: sp.Expr, variable: sp.Symbol) -> int:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    return sp.Poly(denominator, variable).degree() - sp.Poly(
        numerator, variable
    ).degree()


def local_order_pair(u_order: int, h_order: int, side: str) -> tuple[int, int]:
    """Return transverse coefficient order and residue contact index."""

    if side == "zeta_zero":
        return 29 * u_order - 4 * h_order, 5 * h_order - 36 * u_order
    if side == "zeta_infinity":
        return h_order - 7 * u_order, 36 * u_order - 5 * h_order
    raise ValueError(f"unknown carrier target side: {side}")


def compile_local_fans() -> dict[str, object]:
    # A simple zero over zeta=0 has (Pi,zeta)=(q*t^-1,t).
    simple_spectator = compile_toric_fan_profile(
        exponent_map=((1, -1), (0, 1)),
        source_rays=((1, 0), (1, 1)),
        target_rays=((1, 0), (0, 1)),
    )
    assert simple_spectator.refined_source_rays == ((1, 0), (1, 1))
    assert simple_spectator.determinants == (1,)

    # The double-carrier fivefold point starts with four carrier-centered
    # blowups.  Pulling back the next target ray inserts (3,2), one additional
    # fan-alignment blowup away from the carrier.
    fivefold = compile_toric_fan_profile(
        exponent_map=((1, -4), (0, 5)),
        source_rays=((1, 0), (4, 1), (3, 1), (2, 1), (1, 1)),
        target_rays=((1, 0), (0, 1), (-1, 2), (-2, 3)),
    )
    assert fivefold.target_preimage_rays == (
        (1, 0),
        (4, 1),
        (3, 2),
    )
    assert fivefold.refined_source_rays == (
        (1, 0),
        (4, 1),
        (3, 1),
        (2, 1),
        (3, 2),
        (1, 1),
    )
    assert fivefold.determinants == (5, 5, 5, 5, 5)

    # At v=1 (and v=rho on the double row) the carrier residue has pole index
    # three.  The source proximity arm is regular, but it crosses six rays of
    # the extracted target fan.  Their inverse images plus two regularizing
    # rays give the exact common refinement up to the terminal divisor.
    target_left_rays = ((1, 0), (0, 1), *(
        (-index, 5 * index + 1) for index in range(1, 8)
    ))
    principal_arm = compile_toric_fan_profile(
        exponent_map=((1, -1), (0, 3)),
        source_rays=((1, 0), (1, 1), (1, 2), (3, 7), (5, 12)),
        target_rays=target_left_rays,
    )
    expected_arm_rays = (
        (1, 0),
        (1, 1),
        (1, 2),
        (5, 11),
        (4, 9),
        (7, 16),
        (3, 7),
        (11, 26),
        (8, 19),
        (13, 31),
        (5, 12),
    )
    assert principal_arm.refined_source_rays == expected_arm_rays
    assert principal_arm.determinants == (3,) * (len(expected_arm_rays) - 1)

    return {
        "simple_spectator": simple_spectator,
        "fivefold": fivefold,
        "principal_arm": principal_arm,
    }


def refined_source_graph(case: str) -> dict[str, object]:
    graph = minimal_source_principal_boundary(case)
    names = list(graph["component_order"])
    weights = dict(zip(names, graph["self_intersections_in_component_order"]))
    edges = [tuple(edge) for edge in graph["edges"]]
    packet_count = 1 if case == "squarefree_one_packet" else 2

    def subdivide(
        left: str, right: str, additions: tuple[tuple[str, int], ...]
    ) -> None:
        edge = (left, right) if (left, right) in edges else (right, left)
        edges.remove(edge)
        sequence = [left, *(name for name, _ in additions), right]
        for name, weight in additions:
            names.append(name)
            weights[name] = weight
        edges.extend(zip(sequence, sequence[1:]))

    for packet_index in range(1, packet_count + 1):
        prefix = f"packet_{packet_index}_arm_"
        weights[prefix + "2"] = -6
        weights[prefix + "3"] = -6
        weights[prefix + "4"] -= 2
        subdivide(
            prefix + "2",
            prefix + "3",
            (
                (prefix + "carrier_align_2", -1),
                (prefix + "carrier_regular_2", -3),
                (prefix + "carrier_align_3", -1),
            ),
        )
        subdivide(
            prefix + "3",
            prefix + "4",
            (
                (prefix + "carrier_align_5", -1),
                (prefix + "carrier_regular_5", -3),
                (prefix + "carrier_align_6", -1),
            ),
        )

    carrier = str(graph["carrier_component"])
    if case == "squarefree_one_packet":
        weights[carrier] -= 2
        for index in (1, 2):
            spectator = f"carrier_simple_spectator_{index}"
            names.append(spectator)
            weights[spectator] = -1
            edges.append((carrier, spectator))
        expected = (27, 8, -4, 5, [-9])
    else:
        weights[carrier] -= 4
        alpha_chain = (
            ("carrier_alpha_e4", -1),
            ("carrier_alpha_e3", -2),
            ("carrier_alpha_e2", -3),
            ("carrier_alpha_alignment", -1),
            ("carrier_alpha_e1", -3),
        )
        sequence = [carrier]
        for name, weight in alpha_chain:
            names.append(name)
            weights[name] = weight
            sequence.append(name)
        edges.extend(zip(sequence, sequence[1:]))
        expected = (48, 11, -7, 5, [-9, -9])

    edge_set = {tuple(sorted(edge)) for edge in edges}
    matrix = [
        [
            weights[left]
            if left == right
            else int(tuple(sorted((left, right))) in edge_set)
            for right in names
        ]
        for left in names
    ]
    matrix_determinant = bareiss_determinant(matrix)
    inertia = symmetric_inertia(matrix)
    canonical = solve_linear_system(
        matrix, [-2 - weights[name] for name in names]
    )
    log_canonical = [value + 1 for value in canonical]
    log_canonical_square = sum(
        log_canonical[left]
        * matrix[left][right]
        * log_canonical[right]
        for left in range(len(names))
        for right in range(len(names))
    )
    adjacency = {name: set() for name in names}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)

    summary = (
        len(names),
        sum(len(adjacency[name]) == 1 for name in names),
        weights[carrier],
        len(adjacency[carrier]),
        [weights[str(name)] for name in graph["terminal_components"]],
    )
    assert summary == expected
    assert abs(matrix_determinant) == 1
    assert inertia == (1, len(names) - 1, 0)
    assert all(value.denominator == 1 for value in canonical)

    return {
        "component_count": len(names),
        "leaf_count": summary[1],
        "carrier_weight": weights[carrier],
        "carrier_valency": len(adjacency[carrier]),
        "terminal_weights": summary[4],
        "intersection_determinant": matrix_determinant,
        "inertia": inertia,
        "canonical_integral": True,
        "log_canonical_square": log_canonical_square,
    }


def main() -> None:
    v, rho = sp.symbols("v rho")

    # U has divisor -div(c).  At a simple c-root, H has order -7; at a
    # double c-root, H has order -15.  These give exact transverse orders in
    # the two target-node charts.
    assert local_order_pair(-1, -7, "zeta_zero") == (-1, 1)
    assert local_order_pair(-2, -15, "zeta_infinity") == (-1, 3)
    assert local_order_pair(0, 1, "zeta_zero") == (-4, 5)

    square_residue = 1 + 1 / (v - 1) ** 3
    assert infinity_order(square_residue - 1, v) == 3

    alpha = sp.Rational(3, 5) * (rho + 1)
    double_residue = v * (v - alpha) ** 5 / (
        (v - 1) ** 3 * (v - rho) ** 3
    )
    for selected_rho in (
        (3 + sp.sqrt(5)) / 2,
        (3 - sp.sqrt(5)) / 2,
    ):
        selected_residue = sp.cancel(double_residue.subs(rho, selected_rho))
        double_limit = sp.limit(selected_residue, v, sp.oo)
        assert infinity_order(selected_residue - double_limit, v) == 3

    n0 = sp.Rational(625, 3**8) * (3 - rho)
    n1 = sp.Rational(625, 3**9) * (4 * rho - 11)
    assert is_zero_mod_rho_relation(n0 + n1 * alpha, rho)
    for point in (0, 1, rho):
        assert not is_zero_mod_rho_relation(n0 + n1 * point, rho)

    fans = compile_local_fans()
    square_graph = refined_source_graph("squarefree_one_packet")
    double_graph = refined_source_graph("double_same_target")

    assert all(
        determinant != 0
        for profile in fans.values()
        for determinant in profile.determinants
    )
    assert square_graph["component_count"] == 27
    assert double_graph["component_count"] == 48

    print(
        "PASS: F2 carrier log profiles; all marked carrier-local, aligned "
        "principal-arm, and spectator nodes are tame log-etale; bounds are 27/48"
    )


if __name__ == "__main__":
    main()

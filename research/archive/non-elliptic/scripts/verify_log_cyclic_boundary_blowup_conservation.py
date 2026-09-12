#!/usr/bin/env python3
"""Verify cyclic SNC matching/Chern conservation under boundary blowups."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "plane-jc" / "cas"))

from jcsearch.log_node_profiles import cyclic_boundary_charge  # noqa: E402
from verify_f2_75_125_global_attachment import (  # noqa: E402
    source_proximity_fan_audit,
)


def symbolic_blowup_identities() -> None:
    a, b, source_left, source_right = sp.symbols(
        "a b source_left source_right"
    )
    old_node = a**2 * source_left + b**2 * source_right + 2 * a * b
    exceptional = a + b
    new_node = (
        a**2 * (source_left - 1)
        + b**2 * (source_right - 1)
        - exceptional**2
        + 2 * a * exceptional
        + 2 * b * exceptional
    )
    assert sp.expand(new_node - old_node) == 0

    multiplicity, source_self = sp.symbols("multiplicity source_self")
    old_smooth = multiplicity**2 * source_self
    new_smooth = (
        multiplicity**2 * (source_self - 1)
        - multiplicity**2
        + 2 * multiplicity**2
    )
    assert sp.expand(new_smooth - old_smooth) == 0


def blow_up_node(
    multiplicities: tuple[int, ...],
    weights: tuple[int, ...],
    edges: tuple[tuple[int, int], ...],
    selected_edge: tuple[int, int],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[tuple[int, int], ...]]:
    left, right = tuple(sorted(selected_edge))
    edge_set = {tuple(sorted(edge)) for edge in edges}
    if (left, right) not in edge_set:
        raise ValueError("the selected components do not meet")
    exceptional = len(multiplicities)
    new_multiplicities = (*multiplicities, multiplicities[left] + multiplicities[right])
    new_weights = list(weights)
    new_weights[left] -= 1
    new_weights[right] -= 1
    new_weights.append(-1)
    edge_set.remove((left, right))
    edge_set.update({tuple(sorted((left, exceptional))), tuple(sorted((right, exceptional)))})
    return new_multiplicities, tuple(new_weights), tuple(sorted(edge_set))


def exact_graph_regression() -> None:
    multiplicities = (2, 3, 5, 7)
    weights = (-2, -3, -4, -1)
    edges = ((0, 1), (1, 2), (2, 3))
    initial = cyclic_boundary_charge(multiplicities, weights, edges)
    for selected in ((1, 2), (1, 4), (2, 4)):
        multiplicities, weights, edges = blow_up_node(
            multiplicities, weights, edges, selected
        )
        transformed = cyclic_boundary_charge(multiplicities, weights, edges)
        assert transformed.doubled_charge == initial.doubled_charge


def f2_extraction_charge() -> None:
    source = source_proximity_fan_audit()["standard_P2_carrier_fan"]
    weights = tuple(source["boundary_self_intersections"])
    # carrier_chain_5 is the first exceptional (-1,0); carrier_chain_6 is the
    # strict line at infinity (-1,-1).
    extraction_root_weights = (weights[5], weights[6])
    assert extraction_root_weights == (-6, 0)

    root = cyclic_boundary_charge((3, 18), extraction_root_weights, ((0, 1),))
    assert root.node_matching_length == 54
    assert root.doubled_charge == 54
    assert root.charge == 27
    assert root.cyclic_cokernel_ch2 == -27

    blown_data = blow_up_node((3, 18), extraction_root_weights, ((0, 1),), (0, 1))
    blown = cyclic_boundary_charge(*blown_data)
    assert blown.multiplicities == (3, 18, 21)
    assert blown.self_intersections == (-7, -1, -1)
    assert blown.node_matching_length == 3 * 21 + 18 * 21 == 441
    assert blown.charge == root.charge == 27
    assert blown.cyclic_cokernel_ch2 == -27


def main() -> None:
    symbolic_blowup_identities()
    exact_graph_regression()
    f2_extraction_charge()
    print(
        "PASS: cyclic SNC matching plus component self-intersection is "
        "blowup-conserved; the F2 extraction-root corrected charge is 27"
    )


if __name__ == "__main__":
    main()

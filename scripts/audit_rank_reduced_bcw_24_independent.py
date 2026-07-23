#!/usr/bin/env python3
"""Dependency-free sparse audit of the 24D rank-17 cubic witness."""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path

from audit_index_reduced_bcw_22_independent import (
    decode_h,
    matrix_product,
    rational_rank,
)
from audit_shared_bcw_33_independent import derivative, evaluate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "rank_reduced_bcw_24_counterexample.json"
)


def main() -> None:
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "rank-reduced-bcw-sparse-cubic-homogeneous-map-v1"
    dimension = int(stored["dimension"])
    assert dimension == 24
    h = decode_h(stored)
    assert len(h) == dimension
    assert all(sum(exponents) == 3 for poly in h for exponents in poly)

    points = [[Q(value) for value in point] for point in stored["collision_points"]]
    common_image = [Q(value) for value in stored["common_image"]]
    assert len({tuple(point) for point in points}) == 3
    for point in points:
        assert [
            coordinate + evaluate(poly, point)
            for coordinate, poly in zip(point, h)
        ] == common_image

    jacobian = [
        [derivative(h[row], column) for column in range(dimension)]
        for row in range(dimension)
    ]
    point = [Q(index * index + 3 * index + 5) for index in range(dimension)]
    assert rational_rank(
        [[evaluate(entry, point) for entry in row] for row in jacobian]
    ) == 17

    power = jacobian
    nonzero_counts: list[int] = []
    term_counts: list[int] = []
    for _ in range(1, 19):
        nonzero_counts.append(sum(bool(entry) for row in power for entry in row))
        term_counts.append(sum(len(entry) for row in power for entry in row))
        if len(nonzero_counts) < 18:
            power = matrix_product(power, jacobian)
    assert nonzero_counts[-2:] == [8, 0]
    assert term_counts[-1] == 0

    statistics = stored["statistics"]
    assert statistics["generic_rank_JH_over_QQ_x"] == 17
    assert statistics["nilpotency_index_JH"] == 18
    assert statistics["specialized_jordan_type"] == [18, 1, 1, 1, 1, 1, 1]
    certificate = statistics["exact_certificate"]
    assert certificate["independent_generic_kernel_columns"] == 7
    assert certificate["specialized_rank_lower_bound"] == 17
    assert certificate["nonzero_entries_JH_power_17"] == 8
    assert certificate["nonzero_entries_JH_power_18"] == 0

    print("PASS (stdlib) rank-reduced BCW 24: parsed a cubic-homogeneous 24D correction")
    print("PASS (stdlib) rank-reduced BCW 24: verified three exact collision points")
    print("PASS (stdlib) rank-reduced BCW 24: exact specialized Jacobian rank is 17")
    print("PASS (stdlib) rank-reduced BCW 24: (JH)^17 has 8 nonzero entries and (JH)^18=0")
    print("PASS (stdlib) rank-reduced BCW 24: direct nilpotence certifies the Keller property")


if __name__ == "__main__":
    main()

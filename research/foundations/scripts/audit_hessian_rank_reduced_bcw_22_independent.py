#!/usr/bin/env python3
"""Independent final-artifact audit of the rank-37 quartic HN source."""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path

from audit_fixed_rank_hessian_witness import (
    cotangent_hessian,
    deterministic_point,
    exact_singular_certificate,
    nilpotency_profile,
    specialization_profile,
)
from audit_index_reduced_bcw_22_independent import (
    matrix_product,
    rational_rank,
)
from audit_shared_bcw_33_independent import dense, derivative, evaluate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_reduced_bcw_22_counterexample.json"
)


def decode_h(stored: dict[str, object]) -> list[dict[tuple[int, ...], Q]]:
    dimension = int(stored["dimension"])
    return [
        {
            dense(term["monomial"], dimension): Q(term["coefficient"])
            for term in component
        }
        for component in stored["H"]
    ]


def main() -> None:
    stored = json.loads(ARTIFACT.read_text())
    assert (
        stored["format"]
        == "hessian-rank-reduced-bcw-sparse-cubic-homogeneous-map-v1"
    )
    dimension = int(stored["dimension"])
    assert dimension == 22
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

    jacobian_small = [
        [derivative(h[row], column) for column in range(dimension)]
        for row in range(dimension)
    ]
    point = [Q(index * index + 3 * index + 5) for index in range(dimension)]
    evaluated = [
        [evaluate(entry, point) for entry in row] for row in jacobian_small
    ]
    assert rational_rank(evaluated) == 18

    power = jacobian_small
    nonzero_counts = []
    for exponent in range(1, 19):
        nonzero_counts.append(sum(bool(entry) for row in power for entry in row))
        if exponent < 18:
            power = matrix_product(power, jacobian_small)
    assert nonzero_counts[-2:] == [7, 0]

    jacobian, upper_left, hessian = cotangent_hessian(h)
    profiles = [
        specialization_profile(jacobian, upper_left, hessian, 20_260_730 + offset)
        for offset in range(3)
    ]
    assert profiles == [(18, 16, 37)] * 3
    exact_point = deterministic_point(44, 1_000_003, 20_260_730)
    syzygies, specialized_rank, kernel_check = exact_singular_certificate(
        hessian,
        exact_point,
        expected_kernel_rank=7,
    )
    assert syzygies == 12
    assert specialized_rank == 37
    assert kernel_check
    transformed_power_ranks = nilpotency_profile(
        jacobian, upper_left, 20_260_730
    )
    assert transformed_power_ranks == [
        37,
        34,
        *range(32, -1, -1),
    ]

    statistics = stored["statistics"]
    lift = stored["homogeneous_quartic_HN_lift"]
    assert lift["field"] == "QQ(I)"
    assert lift["dimension"] == 44
    assert lift["degree"] == 4
    assert lift["generic_hessian_rank"] == 37
    assert statistics["generic_rank_JH_over_QQ_x"] == 18
    assert statistics["nilpotency_index_JH"] == 18
    assert statistics["generic_cotangent_hessian_rank_over_QQ_xy"] == 37
    assert (
        statistics["cotangent_hessian_exact_certificate"][
            "independent_generic_kernel_columns"
        ]
        == 7
    )
    assert (
        statistics["transformed_hessian_sampled_power_ranks_mod_1000033"]
        == transformed_power_ranks
    )
    assert statistics["transformed_hessian_sampled_nilpotency_index"] == 35

    print("PASS (independent) rank-37 HN: parsed a cubic-homogeneous 22D source")
    print("PASS (independent) rank-37 HN: verified three exact collision points")
    print("PASS (independent) rank-37 HN: (JH)^17!=0 and (JH)^18=0")
    print("PASS (independent) rank-37 HN: three new specializations have rank 37")
    print("PASS (independent) rank-37 HN: 7 exact syzygies certify generic corank 7")
    print("PASS (independent) rank-37 HN: sampled transformed index is 35")


if __name__ == "__main__":
    main()

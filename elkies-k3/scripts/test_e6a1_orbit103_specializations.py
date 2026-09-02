#!/usr/bin/env python3
"""Small exact regressions for the orbit-103 specialization pipeline."""

from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).with_name("search_e6a1_orbit103_specializations.py")
PROBE_ARTIFACT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-e6a1-orbit103-specialization-rank-probes-v1.json"
)
SPEC = importlib.util.spec_from_file_location("orbit103_search", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SEARCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SEARCH)


def direct_affine_coefficients(k: int, r: int, prime: int) -> tuple[int, int]:
    a4 = (
        k**4 * r**4
        + 16 * k * (k**2 + 12) * r**3
        + 8 * k**2 * (k**2 + 152) * r**2
        + 64 * k * (31 * k**2 - 12) * r
        + 16 * k**2 * (61 * k**2 - 48)
    )
    b6 = (
        k**6 * r**6
        + 24 * k**3 * (k**2 + 12) * r**5
        + 12 * (k**6 + 160 * k**4 + 192 * k**2 + 1152) * r**4
        + k * (3072 * k**4 + 35072 * k**2 + 27648) * r**3
        + (1488 * k**6 + 91776 * k**4 + 4608 * k**2 - 55296) * r**2
        + k * (85632 * k**4 - 13824 * k**2 - 110592) * r
        + k**2 * (26560 * k**4 - 4608 * k**2 - 55296)
    )
    return (
        -27 * (r**2 - 4) ** 2 * a4 % prime,
        54 * (r**2 - 4) ** 3 * b6 % prime,
    )


def on_curve(model: list[str], point: list[str]) -> bool:
    a1, a2, a3, a4, a6 = map(Fraction, model)
    x_coordinate, y_coordinate = map(Fraction, point)
    return (
        y_coordinate**2 + a1 * x_coordinate * y_coordinate + a3 * y_coordinate
        == x_coordinate**3 + a2 * x_coordinate**2 + a4 * x_coordinate + a6
    )


class Orbit103SpecializationTest(unittest.TestCase):
    def test_bihomogeneous_affine_evaluation(self) -> None:
        for prime, k, r in ((19, 3, 5), (23, 7, 11), (31, 2, 17)):
            self.assertEqual(
                SEARCH.orbit103_coefficients_mod((k, 1), (r, 1), prime),
                direct_affine_coefficients(k, r, prime),
            )

    def test_simultaneous_sign_symmetry(self) -> None:
        for prime in (19, 29, 41):
            self.assertEqual(
                SEARCH.orbit103_coefficients_mod((5, 7), (11, 13), prime),
                SEARCH.orbit103_coefficients_mod((-5, 7), (-11, 13), prime),
            )

    def test_exact_known_points(self) -> None:
        record = SEARCH.exact_fibre_record((1, 1), (-83, 223))
        self.assertTrue(record["nonsingular"])
        self.assertEqual(len(record["known_points_Q_plus_Q_minus"]), 2)

    def test_pinned_probe_lower_bounds(self) -> None:
        document = json.loads(PROBE_ARTIFACT.read_text())
        counts: dict[int, int] = {}
        for record in document["records"]:
            rank = record["certified_rank_lower_bound"]
            counts[rank] = counts.get(rank, 0) + 1
            points = (
                record["baseline_points"]
                + record["finite_quotient_independent_new_points"]
            )
            self.assertEqual(len(points), rank)
            self.assertTrue(
                record["combined_finite_quotient_certificate"]["certified_independent"]
            )
            self.assertTrue(
                all(on_curve(record["short_integral_model"], point) for point in points)
            )
        self.assertEqual(counts, {2: 33, 3: 251, 4: 366, 5: 240, 6: 103, 7: 7})
        self.assertEqual(document["best_certified_rank_lower_bound"], 7)
        self.assertEqual(document["fibres_of_certified_rank_at_least_8"], 0)
        self.assertEqual(len(document["errors"]), 0)
        self.assertEqual(
            sum(
                record["pari_ellrank"]["backend"].startswith("system_gp_fallback")
                for record in document["records"]
            ),
            617,
        )


if __name__ == "__main__":
    unittest.main()

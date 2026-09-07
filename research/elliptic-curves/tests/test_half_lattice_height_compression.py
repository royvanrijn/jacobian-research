from __future__ import annotations

from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_height_compression_analysis_v1.json.gz"
)
SCRIPT = CAS / "analyze_half_lattice_height_compression.sage"
sys.path.insert(0, str(CAS))

from alternate_quartic_covers import alternate_cover, short_add  # noqa: E402


def load_artifact() -> dict:
    with gzip.open(ARTIFACT, "rt") as handle:
        return json.load(handle)


class HalfLatticeHeightCompressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = load_artifact()

    def test_artifact_is_bound_to_analyzer(self) -> None:
        key = "elliptic-curves/cas/analyze_half_lattice_height_compression.sage"
        self.assertEqual(
            self.data["input_hashes"][key], hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(
            self.data["status"],
            "PASS_EXACT_IDENTITIES_AND_BOUNDED_RETROSPECTIVE_MECHANISM_AUDIT",
        )

    def test_complete_detailed_chart_census(self) -> None:
        census = self.data["chart_census"]
        self.assertEqual(census["chart_count"], 3865)
        counts = {row["id"]: row["chart_count"] for row in census["datasets"]}
        self.assertEqual(
            counts,
            {
                "rank28-selected-union": 64,
                "curve385-first-quotient-iteration": 301,
                "curve398-initial-deepest12": 12,
                "curve398-first-quotient-iteration": 372,
                "curve385-primary-natural-weight-1": 394,
                "curve385-primary-natural-weight-2": 2722,
            },
        )
        self.assertTrue(
            all(
                row["raw_binary_quartic_invariants"]["constant_across_every_chart"]
                and row["reduced_discriminant_quartic_invariants"][
                    "constant_across_every_chart"
                ]
                for row in census["datasets"]
            )
        )

    def test_rank28_coordinate_visibility_replays_sources_and_gains(self) -> None:
        audit = self.data["rank28_posthoc_target_oracle"]
        self.assertEqual(
            audit["pairwise_visibility_confusion"],
            {
                "true_positive": 40,
                "false_positive": 0,
                "false_negative": 0,
                "true_negative": 2520,
            },
        )
        self.assertTrue(
            audit["reduced_coordinate_visibility_matches_recorded_sources_exactly"]
        )
        self.assertTrue(
            audit["coordinate_visibility_reconstructs_every_prefix_quotient_gain"]
        )
        self.assertEqual(audit["reconstructed_final_quotient_rank"], 11)
        gap = audit["centered_height_strict_separation"]
        self.assertGreater(
            float(gap["minimum_over_nonproductive_charts"]),
            float(gap["maximum_over_productive_charts"]),
        )

    def test_no_target_free_scalar_is_stable(self) -> None:
        conclusion = self.data["predictor_conclusion"]
        self.assertEqual(conclusion["positive_chart_order_count"], 11)
        self.assertEqual(
            conclusion["stable_target_free_scalar_predictors_at_auc_at_least_0_7"],
            [],
        )

    def test_compact_positive_control_prefix_gains(self) -> None:
        panel = self.data["compact_control_predictor_panel"]
        self.assertEqual(panel["case_count"], 7)
        self.assertEqual(panel["chart_count"], 394)
        recovered = {
            row["id"]: row["blind_recovered_quotient_dimension"]
            for row in panel["cases"]
        }
        self.assertEqual(
            recovered,
            {
                "r17-control-a": 8,
                "r17-control-b": 9,
                "r17-control-c": 10,
                "r17-control-d": 9,
                "curve12-2024-rank29": 10,
                "curve356-rank29": 12,
                "curve385-rank29": 4,
            },
        )

    def test_exact_fiber_coordinate_identity(self) -> None:
        blind = json.loads(
            (
                ROOT
                / "artifacts/generated-results/elliptic-curves"
                / "half_lattice_fake_descent_rank28_blind_v1.json"
            ).read_text()
        )
        model = tuple(Fraction(value) for value in blind["fibre"]["short_model"])
        record = blind["search_records"][0]
        base = (
            Fraction(record["base_point"]["x"]),
            Fraction(record["base_point"]["y"]),
        )
        target = (
            Fraction(record["finite_curve_points"][0]["x"]),
            Fraction(record["finite_curve_points"][0]["y"]),
        )
        parameter = Fraction(record["finite_raw_points"][0]["x"])
        companion = short_add(model, base, (target[0], -target[1]))
        self.assertIsNotNone(companion)
        self.assertEqual(
            parameter**2,
            target[0] + companion[0] + base[0],
        )
        cover = alternate_cover(model, base)
        self.assertEqual(
            cover.cover_point_to_curve(
                (
                    parameter,
                    Fraction(record["finite_raw_points"][0]["y"]),
                )
            ),
            target,
        )


if __name__ == "__main__":
    unittest.main()

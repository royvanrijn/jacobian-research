#!/usr/bin/env python3
"""Focused checks for the mapped-height-one-million pair-product charts."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_fermigier_published_pair_fiber_products_h1000000_charts.py"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_published_pair_fiber_products_h1000000_charts.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "c24323961ab2877e57e1efbf61cd73efafb2aa18f05bb3e3206a5328942b2e35"
)
EXPECTED_ARTIFACT_SHA256 = (
    "615f32c9c9c9ed2b9cde8f7ef7ebed60aed75b9d80104f6d36330b4031f62b25"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierPairFiberProductH1000000ChartsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.rows = [
            row
            for chart in cls.data["chart_searches"]
            for row in chart["pair_searches"]
        ]

    def test_pinned_files_and_stable_exact_dependencies(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        source = self.data["source"]
        self.assertEqual(
            source["published_preimage_sha256"],
            "6224da9ce4db3150a197a2cf1d9bc6c1a7d0cc6f01245b3f834945f76775ab15",
        )
        self.assertEqual(
            source["H50000_exact_pair_result_sha256"],
            "dea8b716c5aec56817a172afd6e894e7748aaddc482a2d29c0a3360abe55bf4b",
        )
        prior = self.data["prior_decontamination"]
        self.assertEqual(prior["base_prior_parameter_count"], 590)
        self.assertEqual(
            prior["base_prior_parameter_sha256"],
            "64c09a13b427938a44251a91f74a116f7f9e685aed07c6159550e7ec3ea51291",
        )
        self.assertEqual(prior["terminal_prior_parameter_count"], 593)
        self.assertEqual(
            prior["terminal_prior_parameter_sha256"],
            "a4d06e4662d2e30c1a0f8873f91d8d348dae10f2abaffce88dcc0f480cfeede0",
        )

    def test_chart_union_is_exact_bounded_and_not_overclaimed(self) -> None:
        scope = self.data["scope"]
        self.assertEqual(scope["chart_height"], 50_000)
        self.assertEqual(scope["mapped_projective_height_ceiling"], 1_000_000)
        self.assertFalse(scope["full_H1000000_exhaustive"])
        charts = self.data["chart_population"]
        self.assertEqual(len(charts), 24)
        self.assertEqual(
            Counter(chart["kind"] for chart in charts),
            Counter({"translate": 12, "reciprocal-shift": 12}),
        )
        self.assertEqual(len({chart["chart_id"] for chart in charts}), 24)
        self.assertEqual(
            max(chart["mapped_projective_height_upper_bound"] for chart in charts),
            1_000_000,
        )
        for chart in charts:
            a_value, b_value, c_value, d_value = chart["matrix"]
            self.assertEqual(abs(a_value * d_value - b_value * c_value), 1)
            self.assertLessEqual(
                chart["mapped_projective_height_upper_bound"], 1_000_000
            )

    def test_all_5280_searches_completed_once(self) -> None:
        self.assertEqual(len(self.data["chart_searches"]), 24)
        self.assertEqual(len(self.rows), 5_280)
        self.assertTrue(
            all(row["search"]["search"]["status"] == "completed" for row in self.rows)
        )
        self.assertTrue(
            all(not row["search"]["search"]["retried"] for row in self.rows)
        )
        self.assertTrue(
            all(row["search"]["search"]["height_bound"] == 50_000 for row in self.rows)
        )
        outcome = self.data["outcome"]
        self.assertEqual(outcome["declared_search_call_count"], 5_280)
        self.assertEqual(outcome["search_calls_attempted"], 5_280)
        self.assertEqual(outcome["search_calls_completed"], 5_280)
        self.assertEqual(outcome["search_calls_timed_out_or_errored"], 0)
        self.assertEqual(self.data["execution"]["charts_completed"], 24)

    def test_every_noncalibration_incidence_fails_individual_square_gate(self) -> None:
        classifications = Counter()
        product_only = []
        for row in self.rows:
            for incidence in row["search"]["incidences"]:
                classifications[incidence["classification"]] += 1
                if incidence["classification"] == "record-fiber-excluded":
                    self.assertEqual(incidence["canonical_parameter_t"], "39508/39")
                    self.assertTrue(incidence["left_factor_is_square"])
                    self.assertTrue(incidence["right_factor_is_square"])
                    self.assertEqual(len(incidence["exact_forced_quartic_points"]), 2)
                    self.assertTrue(
                        all(
                            point["exact_membership_checked"]
                            for point in incidence["exact_forced_quartic_points"]
                        )
                    )
                elif incidence["classification"] == "product-square-only":
                    product_only.append(incidence)
                    self.assertFalse(incidence["left_factor_is_square"])
                    self.assertFalse(incidence["right_factor_is_square"])
                else:
                    self.fail(f"unexpected incidence classification: {incidence}")
        self.assertEqual(
            classifications,
            Counter({"record-fiber-excluded": 2_640, "product-square-only": 38}),
        )
        self.assertEqual(
            {abs(Fraction(incidence["signed_parameter_t"])) for incidence in product_only},
            {Fraction(48_363, 26), Fraction(23_317, 6), Fraction(42_058, 25)},
        )
        self.assertEqual(
            self.data["outcome"]["incidence_classification_counts"],
            {"product-square-only": 38, "record-fiber-excluded": 2_640},
        )

    def test_negative_target_frontier_is_exactly_scoped(self) -> None:
        self.assertEqual(self.data["candidates"], [])
        outcome = self.data["outcome"]
        self.assertEqual(outcome["genuinely_new_double_forced_fibers"], 0)
        self.assertEqual(outcome["completed_conductors"], 0)
        self.assertEqual(outcome["subtarget_conductors"], 0)
        self.assertEqual(outcome["rank_triage_count"], 0)
        self.assertIsNone(outcome["maximum_stable_numerical_rank"])
        self.assertEqual(
            outcome["exact_chart_result_sha256"],
            "e8a4421d6d9b66667b6fbcdb25eff179f719d638e87bde7771460e8bb4b5cdaf",
        )
        self.assertEqual(outcome["maximum_mapped_projective_height_observed"], 48_363)
        self.assertFalse(self.data["target"]["hit"])
        self.assertTrue(self.data["parameters"]["no_retries"])


if __name__ == "__main__":
    unittest.main()

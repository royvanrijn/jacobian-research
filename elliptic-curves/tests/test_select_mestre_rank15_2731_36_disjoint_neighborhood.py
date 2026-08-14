#!/usr/bin/env python3
"""Focused checks for the frozen disjoint T=2731/36 neighborhood."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/select_mestre_rank15_2731_36_disjoint_neighborhood.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results"
    / "elliptic_mestre_rank15_2731_36_disjoint_neighborhood_selection.json"
)
EXPECTED_SCRIPT_SHA256 = "31e1e73818855679239ea0b336ce8a9838d252b0c69339a3eadee870933e1995"
EXPECTED_ARTIFACT_SHA256 = "512c0fdaa47ccc23794c104c9ad4e96af5b286b4bf46815cafb871a44e931c63"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRank15DisjointNeighborhoodSelectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_result(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["provenance"]["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(
            self.data["result_sha256"],
            "fb7dd7cef970e4e7d84c5120e37326d7c3258baa3067605e46adf64a96ba0a1a",
        )
        self.assertEqual(self.data["provenance"]["external_conductor_point_height_rank_calls"], 0)

    def test_exact_disjoint_population_and_generator_manifest(self) -> None:
        scope = self.data["scope"]
        self.assertEqual(scope["roots"], [0, 7, 93, 154, 161, 191])
        self.assertEqual(scope["center"], "2731/36")
        self.assertTrue(scope["center_excluded"])
        self.assertEqual(scope["reduced_denominator_interval"], [257, 2048])
        self.assertEqual(scope["prior_common_box_reduced_denominator_maximum"], 256)
        self.assertTrue(scope["exactly_disjoint_from_prior_common_box"])
        generator = self.data["generator"]
        self.assertEqual(generator["unique_reduced_parameter_count"], 57_755)
        self.assertEqual(generator["raw_attempt_count"], 95_584)
        self.assertEqual(
            generator["canonical_parameter_source_manifest_sha256"],
            "c3c9dcbb3b5be8f11ff933ff5f96b0ad5c8cc51aaa4c0d3affcfa4cf6b2e2bf2",
        )
        self.assertEqual(
            generator["parameters_by_source_category_after_dedup"],
            {
                "ordinary-annulus": 52_607,
                "gauss-farey": 2_789,
                "discriminant-power": 0,
                "local-trace": 3_504,
            },
        )
        self.assertTrue(
            all(not roots for roots in generator["simple_hensel_roots_mod_p2"].values())
        )

    def test_fresh_trace_feature_and_rank_blind_selection_gates(self) -> None:
        trace = self.data["local_trace_screen"]
        self.assertEqual(trace["discovery_primes"], [919, 929, 937, 941, 947, 953, 967, 971])
        self.assertEqual(trace["held_primes"], [977, 983, 991, 997, 1009, 1013, 1019, 1021])
        self.assertTrue(trace["bands_disjoint"])
        self.assertTrue(trace["fresh_relative_to_prior_multifamily_bands_through_911"])
        self.assertEqual(trace["complete_population_scored_on_discovery_band"], 57_755)
        self.assertEqual(trace["discovery_survivors_retained"], 2_048)
        self.assertEqual(
            trace["discovery_table_digest"],
            "962b74b72b6e01a70f3e0ee27081f7c2516ab98575ff5fc6d3f21680d93daf8d",
        )
        self.assertEqual(
            trace["held_table_digest"],
            "9a8ee6ba995a4498105316161578d7e19feb00b069bb4c57f35d26f24db2461a",
        )
        features = self.data["exact_discriminant_feature_screen"]
        self.assertEqual(features["admissible_feature_pool_count"], 2_048)
        self.assertEqual(features["exact_singular_rejections"], 0)
        self.assertEqual(
            features["exact_feature_population_sha256"],
            "8291e275848a225661eb952cb3186e49d3ecb3c0042e3e1a2d68e9946c07203e",
        )
        selection = self.data["conductor_selection"]
        self.assertTrue(selection["discovery_population_closed_before_held_scores"])
        self.assertFalse(selection["selection_uses_conductor_or_point_or_rank_data"])
        self.assertEqual(selection["selected_population"], 46)
        self.assertEqual(len(self.data["selected_records"]), 46)
        self.assertEqual(
            selection["selected_population_sha256"],
            "42fa8afe79fe912e3e24f79f33c3ae78922273d85b31fd8335ffe8cb95f72744",
        )


if __name__ == "__main__":
    unittest.main()

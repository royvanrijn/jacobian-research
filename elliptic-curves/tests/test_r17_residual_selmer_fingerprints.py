from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-residual-selmer-fingerprints-v1.json"
)


class R17ResidualSelmerFingerprintTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.document = json.loads(ARTIFACT.read_text())

    def test_claim_boundary_distinguishes_known_localizations_from_selmer(self) -> None:
        summary = self.document["summary"]
        self.assertFalse(summary["complete_two_selmer_groups_computed"])
        self.assertFalse(summary["full_selmer_leave_one_place_out_matrices_computed"])
        self.assertTrue(summary["known_residual_localization_matrices_computed"])

    def test_control_dimensions_and_known_localization_rank(self) -> None:
        expected = {351: 8, 356: 12, 376: 5, 377: 6, 385: 12, 12: 12}
        rows = {int(row["curve_id"]): row for row in self.document["fingerprints"]}
        self.assertEqual(set(rows), set(expected))
        for curve_id, dimension in expected.items():
            row = rows[curve_id]
            self.assertEqual(row["certified_known_residual_dimension"], dimension)
            intersections = row["localization_intersections"]
            self.assertEqual(intersections["full_stacked_localization_rank"], dimension)
            self.assertEqual(
                intersections["full_simultaneous_localization_kernel_dimension"], 0
            )
            self.assertTrue(
                all(item["rank_drop"] == 0 for item in intersections["delete_one_place"])
            )

    def test_two_adic_features_separate_plus_twelve_from_plus_five(self) -> None:
        comparison = {
            int(row["prime"]): row
            for row in self.document["plus12_vs_plus5"]["placewise_comparisons"]
        }
        names = {
            item["feature"]
            for item in comparison[2][
                "scalar_features_separating_both_plus12_samples_from_plus5"
            ]
        }
        self.assertEqual(
            names,
            {
                "kodaira_symbol",
                "minimal_discriminant_valuation",
                "tamagawa_two_part",
                "known_residual_localization_kernel_dimension",
            },
        )

    def test_place_block_support_spectrum_explains_zero_delete_one_ranks(self) -> None:
        expected_minimum_weights = {
            351: (4, 4),
            356: (5, 6),
            376: (7, 7),
            377: (5, 5),
            385: (6, 7),
            12: (4, 4),
        }
        rows = {int(row["curve_id"]): row for row in self.document["fingerprints"]}
        for curve_id, (minimum_weight, selected_minimum_weight) in (
            expected_minimum_weights.items()
        ):
            code = rows[curve_id]["known_residual_place_block_support_code"]
            self.assertEqual(code["minimum_nonzero_place_block_weight"], minimum_weight)
            self.assertEqual(
                code["every_deletion_of_at_most_this_many_places_preserves_injectivity"],
                minimum_weight - 1,
            )
            self.assertEqual(
                sum(code["block_weight_enumerator"].values()),
                (1 << rows[curve_id]["certified_known_residual_dimension"]) - 1,
            )
            self.assertTrue(code["minimum_support_words"])
            self.assertTrue(
                all(
                    len(word["place_block_support"]) == minimum_weight
                    for word in code["minimum_support_words"]
                )
            )
            selected = rows[curve_id][
                "selected_comparison_block_place_support_code"
            ]
            self.assertEqual(
                selected["minimum_nonzero_place_block_weight"],
                selected_minimum_weight,
            )
            self.assertEqual(
                selected["source_dimension"],
                len(rows[curve_id]["selected_comparison_block"]),
            )

    def test_record_ten_blocks_are_place_matroid_indecomposable(self) -> None:
        rows = {int(row["curve_id"]): row for row in self.document["fingerprints"]}
        for curve_id in (356, 385):
            block = rows[curve_id]["selected_comparison_block_dual_components"]
            self.assertEqual(block["source_dimension"], 10)
            self.assertEqual(block["global_localization_dual_rank"], 10)
            self.assertEqual(block["component_count"], 1)
            self.assertTrue(block["place_block_matroid_indecomposable"])
            self.assertIn("not a Cassels-pairing", block["claim_boundary"])

    def test_crt_classes_are_hypotheses_not_inherited_certificates(self) -> None:
        block = self.document["crt_search_prototypes"]
        self.assertIn("does not yet prove", block["claim_boundary"])
        self.assertEqual(
            [row["prototype_curve_id"] for row in block["prototypes"]], [356, 385]
        )
        self.assertTrue(
            all(
                row["status"]
                == "EXACT_CRT_CLASS_HEURISTIC_FINGERPRINT_PRESERVATION"
                for row in block["prototypes"]
            )
        )


if __name__ == "__main__":
    unittest.main()

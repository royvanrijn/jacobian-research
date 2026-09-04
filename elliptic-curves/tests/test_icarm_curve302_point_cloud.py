import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve302_point_cloud_v1.json"
)


class Curve302PointCloudTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        cls.curves = {item["label"]: item for item in cls.payload["curves"]}

    def test_mod2_and_mod3_codes_separate_every_submitted_direction(self) -> None:
        expected = {
            "curve302": (31, 17, 14),
            "curve273": (30, 17, 13),
        }
        for label, (full_rank, first_rank, quotient_rank) in expected.items():
            for code in self.curves[label]["finite_reduction_kummer_codes"].values():
                self.assertEqual(code["combined_image_dimension"], full_rank)
                self.assertEqual(code["combined_relation_dimension"], 0)
                self.assertEqual(code["subset_dimensions"]["first_17"], first_rank)
                self.assertEqual(
                    code["subset_dimensions"]["remaining_dimension_mod_first_17"],
                    quotient_rank,
                )

    def test_exact_coordinate_pattern_search_has_no_hits(self) -> None:
        for curve in self.curves.values():
            patterns = curve["coordinate_patterns"]
            for records in patterns["exact_low_complexity_relations"].values():
                self.assertFalse(records)
            self.assertFalse(patterns["repeated_pair_sums"])
            self.assertFalse(patterns["repeated_oriented_pair_differences"])
            self.assertFalse(patterns["repeated_pair_products"])

    def test_held_out_interpolation_and_deformation_fail_closed(self) -> None:
        for curve in self.curves.values():
            for group in curve["held_out_low_degree_interpolation"]["groups"]:
                self.assertEqual(group["best_held_out_exact_match_count"], {"x": 0, "y": 0})
                self.assertEqual(group["best_joint_xy_held_out_exact_match_count"], 0)
            deformation = curve["fixed_x_quadratic_deformations"]
            self.assertEqual(deformation["maximum_preserved_subset_size"], 2)
            self.assertFalse(deformation["candidates_with_at_least_one_held_out_hit"])

    def test_curve302_denominator_clusters_cross_the_first17_boundary(self) -> None:
        patterns = self.curves["curve302"]["coordinate_patterns"]
        self.assertEqual(
            patterns["repeated_denominator_root_clusters"]["3"], [16, 21, 25]
        )
        self.assertIn(20, patterns["integral_x_point_indices"])
        self.assertIn(23, patterns["integral_x_point_indices"])
        bad_large_clusters = {
            item["prime"]: len(item["point_indices"])
            for item in patterns["small_prime_numerator_divisibility_clusters"]
            if not item["good_reduction_on_short_model"]
        }
        self.assertEqual(bad_large_clusters, {5: 25, 11: 23, 23: 15})

    def test_mestre_quartic_control_has_one_mod2_coset(self) -> None:
        control = self.payload["mestre_two_cover_coset_calibration"]
        self.assertTrue(control["all_thirteen_transport_to_one_mod2_coset"])
        self.assertTrue(control["twelve_visible_points_sum_to_zero_in_public_basis"])
        self.assertEqual(sum(control["common_public_basis_parity_vector"]), 5)
        self.assertEqual(
            self.curves["curve302"]["finite_reduction_kummer_codes"]["mod_2"][
                "combined_relation_dimension"
            ],
            0,
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import unittest

from verify_mestre_diameter235_eight_companion_component import verify_component


class Diameter235EightCompanionComponentTest(unittest.TestCase):
    def test_exact_component(self) -> None:
        result = verify_component()
        certificate = result["residual_identity_certificate"]
        self.assertTrue(certificate["all_recursive_residuals_vanish"])
        self.assertGreater(
            certificate["admissible_exact_sample_count"],
            certificate["cleared_numerator_degree_bound"],
        )
        escape = result["finite_reduction_visible_subgroup_escape_at_seed"]["records"]
        self.assertEqual(
            [(row["T"], row["visible_mod3_rank"], row["visible_plus_eight_companions_mod3_rank"]) for row in escape],
            [("1", 10, 10), ("2", 10, 11), ("3", 10, 11), ("-1", 10, 10)],
        )
        generic_escape = result["generic_visible_subgroup_escape"]
        self.assertTrue(generic_escape["characteristic_zero_fibre_is_smooth"])
        self.assertEqual(generic_escape["visible_mod3_rank"], 10)
        self.assertEqual(
            generic_escape["visible_plus_each_selected_companion_mod3_rank"], [11, 11]
        )
        self.assertEqual(generic_escape["generic_rank_lower_bound"], 11)
        self.assertTrue(
            all(
                item["independent_subset_indices_one_based"] == list(range(1, 11)) + [13]
                for item in generic_escape["individual_exact_rank_certificates"]
            )
        )
        self.assertTrue(
            result["pair_abscissa_collision_at_seed"][
                "common_triangular_orientation_gives_hyperelliptic_conjugates"
            ]
        )


if __name__ == "__main__":
    unittest.main()

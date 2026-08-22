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


if __name__ == "__main__":
    unittest.main()

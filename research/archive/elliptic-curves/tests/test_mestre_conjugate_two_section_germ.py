from __future__ import annotations

import unittest

from probe_mestre_conjugate_two_section_germ import replay


class ConjugateTwoSectionGermTest(unittest.TestCase):
    def test_declared_local_lift(self) -> None:
        result = replay(bivariate_order=3, slice_order=6)
        self.assertEqual(result["exact_tangent_rank_over_Q"], 6)
        self.assertEqual(result["modular_tangent_ranks"], {"17": 6, "23": 6, "29": 6})
        self.assertTrue(
            result["bivariate_germ"]["all_E0_2_coefficients_through_total_order_vanish"]
        )
        self.assertTrue(
            result["intercept_slice"]["all_E0_2_coefficients_through_order_vanish"]
        )
        self.assertEqual(
            result["root_motion_slice"]["recognized_common_intercept"],
            "233/(113-7t)=233/(240-7r3)",
        )
        self.assertEqual(result["comparison_seed"]["exact_tangent_rank_over_Q"], 6)
        self.assertTrue(
            result["comparison_seed"]["all_E0_2_coefficients_through_total_order_vanish"]
        )


if __name__ == "__main__":
    unittest.main()

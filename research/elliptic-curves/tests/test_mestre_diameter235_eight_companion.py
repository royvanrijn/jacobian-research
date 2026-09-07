from __future__ import annotations

import unittest

from probe_mestre_diameter235_eight_companion import replay


class Diameter235EightCompanionTest(unittest.TestCase):
    def test_local_continuation_and_quotient_escape(self) -> None:
        result = replay(order=12, precision=4)
        self.assertEqual(result["exact_jacobian_rank"], 7)
        self.assertEqual(result["all_reconstructed_section_count"], 8)
        self.assertFalse(result["expanded_two_section_residual_materialized"])
        self.assertEqual(
            [
                (item["T"], item["visible_mod3_rank"], item["visible_plus_eight_companions_mod3_rank"])
                for item in result["finite_reduction_quotient_ranks"]
            ],
            [("1", 10, 10), ("2", 10, 11), ("3", 10, 11), ("-1", 10, 10)],
        )
        for record in result["finite_reduction_quotient_ranks"]:
            expected = 11 if record["T"] in {"2", "3"} else 10
            self.assertEqual(
                record["visible_plus_individual_companion_mod3_ranks"], [expected] * 8
            )


if __name__ == "__main__":
    unittest.main()

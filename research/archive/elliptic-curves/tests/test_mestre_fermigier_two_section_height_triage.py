from __future__ import annotations

import unittest

from screen_mestre_fermigier_two_section_height_triage import replay


class FermigierTwoSectionHeightTriageTest(unittest.TestCase):
    def test_small_exact_panel_is_stable(self) -> None:
        result = replay(height=2, parameters=(1,))
        self.assertEqual(result["smooth_specialization_count"], 2)
        self.assertEqual(result["best_numerical_height_rank"], 9)
        self.assertEqual(
            result["best_specializations"],
            [{"u": "-1/2", "T": "1"}, {"u": "1/2", "T": "1"}],
        )


if __name__ == "__main__":
    unittest.main()

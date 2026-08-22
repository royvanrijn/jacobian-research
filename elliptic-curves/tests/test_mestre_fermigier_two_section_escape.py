import unittest

from screen_mestre_fermigier_two_section_escape import replay


class MestreFermigierTwoSectionEscapeTest(unittest.TestCase):
    def test_second_line_has_no_bounded_quotient_escape(self):
        result = replay()
        self.assertEqual(len(result["records"]), 40)
        self.assertEqual(result["mod2_escape_count"], 0)
        self.assertEqual(result["mod3_escape_count"], 0)
        self.assertTrue(
            all(
                record["mod2_augmented_rank"] <= record["mod2_baseline_rank"]
                and record["mod3_augmented_rank"] <= record["mod3_baseline_rank"]
                for record in result["records"]
            )
        )


if __name__ == "__main__":
    unittest.main()

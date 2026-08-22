import unittest

from screen_mestre_transverse_conic_component_independence import screen


class MestreTransverseConicComponentIndependenceTest(unittest.TestCase):
    def test_small_exact_screen(self):
        result = screen(root_height=2, parameter_height=2, prime_bound=31)
        self.assertGreater(result["admissible_specialization_count"], 0)
        self.assertIn("strict_quotient_rank_gains", result)


if __name__ == "__main__":
    unittest.main()

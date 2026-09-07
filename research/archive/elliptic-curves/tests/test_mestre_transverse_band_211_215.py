import unittest

from audit_mestre_transverse_band_211_215 import replay


class MestreTransverseBand211215Test(unittest.TestCase):
    def test_rank_seven_tangents_match_the_known_component(self):
        result = replay()
        self.assertEqual(result["rank_seven_base_point_count"], 2)
        self.assertEqual(result["rank_seven_pair_count"], 6)
        self.assertTrue(
            all(
                pair["matches_known_component_tangent"]
                for audit in result["audits"]
                for pair in audit["pairs"]
            )
        )


if __name__ == "__main__":
    unittest.main()

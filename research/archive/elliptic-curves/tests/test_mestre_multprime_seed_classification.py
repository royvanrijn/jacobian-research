import unittest

from audit_mestre_multprime_seed_classification import replay


class MestreMultprimeSeedClassificationTest(unittest.TestCase):
    def test_all_rank_seven_base_points_are_classified(self):
        result = replay()
        self.assertEqual(result["rank_seven_base_point_count"], 8)
        self.assertEqual(result["rank_seven_pair_count"], 38)
        classifications = {tuple(item["roots"]): item["base_classification"] for item in result["classifications"]}
        self.assertEqual(
            classifications[(0, 8, 58, 77, 85, 102)],
            "Fermigier two-parameter base family at u=-3, v=-8/3",
        )
        self.assertEqual(
            sum(
                kind == "previous rational two-section component base curve"
                for kind in classifications.values()
            ),
            7,
        )


if __name__ == "__main__":
    unittest.main()

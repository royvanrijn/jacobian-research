import unittest

from audit_mestre_transverse_conic_component_relations import replay


class MestreTransverseConicComponentRelationsTest(unittest.TestCase):
    def test_seed_is_rank_neutral(self):
        result = replay()
        self.assertEqual(result["visible_mod3_rank"], 9)
        self.assertEqual(result["augmented_mod3_rank"], 9)
        self.assertEqual(
            result["exact_affine_relations"],
            ["P1=-V1-V4-V5-V7-V10", "P2=-V1-V4-V5"],
        )


if __name__ == "__main__":
    unittest.main()

import unittest

from verify_mestre_transverse_component_generic_relations import replay


class MestreTransverseComponentGenericRelationsTest(unittest.TestCase):
    def test_both_affine_sections_are_generically_visible(self):
        result = replay()
        self.assertEqual(result["base_field"], "Q(r,T)")
        self.assertEqual(
            result["relations"],
            ["P1=-V(0,-)+V(0,+)-V(1,+)", "P2=V(0,-)+V(1,+)+V(r3,-)"],
        )


if __name__ == "__main__":
    unittest.main()

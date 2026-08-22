import unittest

from verify_mestre_transverse_conic_component_generic_relations import replay


class MestreTransverseConicComponentGenericRelationsTest(unittest.TestCase):
    def test_both_affine_sections_are_generically_visible(self):
        result = replay(include_second=True)
        self.assertEqual(result["base_field"], "Q(s,T)")
        self.assertEqual(
            result["relation"],
            "P1=V(0,-)+V(1,+)+V(r3,-)+V(r4,-)+V(r5,+)",
        )
        self.assertEqual(result["second_relation"], "P2=-V(0,-)-V(1,+)-V(r3,-)")


if __name__ == "__main__":
    unittest.main()

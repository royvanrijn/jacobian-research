from __future__ import annotations

import unittest

from audit_mestre_conjugate_two_section_seed_relations import replay


class ConjugateTwoSectionSeedRelationsTest(unittest.TestCase):
    def test_seed_relations(self) -> None:
        result = replay()
        self.assertEqual(len(result["records"]), 2)
        for record in result["records"]:
            self.assertEqual(record["visible_mod3_rank"], 9)
            self.assertEqual(record["augmented_mod3_rank"], 9)


if __name__ == "__main__":
    unittest.main()

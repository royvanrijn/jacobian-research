from __future__ import annotations

import unittest

from audit_mestre_dsquare_all_companion_relations import replay


class MestreDSquareAllCompanionRelationsTest(unittest.TestCase):
    def test_generic_relations(self) -> None:
        result = replay()
        self.assertEqual(
            result["status"], "exact generic D-square companion relations verified"
        )
        self.assertEqual(len(result["remaining_companion_relations"]), 5)
        self.assertFalse(result["expanded_two_section_residual_materialized"])


if __name__ == "__main__":
    unittest.main()

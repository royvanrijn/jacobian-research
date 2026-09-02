from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_relative_2selmer_checkpointed as checkpointed  # noqa: E402


class CheckpointedRelativeSelmerTests(unittest.TestCase):
    def test_workers_preserve_blind_order(self) -> None:
        self.assertIn("bnfinit", checkpointed.BNF_WORKER)
        self.assertIn("bnfcertify", checkpointed.BNF_WORKER)
        self.assertIn("writebin", checkpointed.BNF_WORKER)
        self.assertIn("polredbest", checkpointed.BNF_WORKER)
        self.assertIn("curve_theta_in_field", checkpointed.BNF_WORKER)
        self.assertIn("ell2selmer_basis_gen", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("elllocalimage_mapped", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("nfeltsign", checkpointed.SIMON_GP_FUNCTION)
        self.assertNotIn('payload["exceptional_points"]', checkpointed.SELMER_WORKER)
        self.assertNotIn('payload["generic_points"]', checkpointed.SELMER_WORKER)
        self.assertIn("cover_for", checkpointed.QUOTIENT_COVER_WORKER)
        self.assertNotIn('payload["exceptional_points"]', checkpointed.QUOTIENT_COVER_WORKER)

    def test_standard_extension(self) -> None:
        initial = [[1, 1, 0, 0], [0, 1, 1, 0]]
        extension = checkpointed.extend_standard(initial, 4)
        self.assertEqual(len(extension), 2)
        self.assertEqual(checkpointed.f2_rank([*initial, *extension]), 4)
        self.assertTrue(all(sum(row) == 1 for row in extension))

    def test_control_classification(self) -> None:
        class Case:
            exceptional_points = ((0, 0),) * 2

        selmer = {"two_selmer_dimension": 20}
        generic_rows = [
            [1 if column == row else 0 for column in range(20)]
            for row in range(17)
        ]
        exceptional_rows = [
            [1 if column == row else 0 for column in range(20)]
            for row in (17, 18)
        ]
        result = checkpointed.combine_results(
            Case(),
            selmer,
            {"point_selmer_rows": generic_rows},
            {"point_selmer_rows": exceptional_rows},
            {
                "quotient_basis": [
                    [1 if column == row else 0 for column in range(20)]
                    for row in (17, 18, 19)
                ],
                "enumeration_complete": True,
                "classes": [
                    {
                        "quotient_class_integer": mask,
                        "quotient_bits": [(mask >> index) & 1 for index in range(3)],
                        "blind_search": {"status": "no_point_within_bound"},
                    }
                    for mask in range(1, 8)
                ],
            },
        )
        self.assertEqual(result["quotient_dimension"], 3)
        self.assertEqual(result["exceptional_quotient_rank"], 2)
        self.assertEqual(result["classes_not_realized_by_known_exceptional_subgroup"], 4)
        labeled = {
            row["quotient_class_integer"]: row
            for row in result["quotient_cover_classification"]
        }
        self.assertTrue(labeled[3]["known_exceptional_subgroup_realizes_class"])
        self.assertFalse(labeled[4]["known_exceptional_subgroup_realizes_class"])


if __name__ == "__main__":
    unittest.main()

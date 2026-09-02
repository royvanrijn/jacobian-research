from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_relative_2selmer_open as open_runner  # noqa: E402


class ElkiesRelative2SelmerOpenTests(unittest.TestCase):
    def test_worker_is_blind_and_uses_certified_pari_cover_basis(self) -> None:
        source = open_runner.PARI_WORKER
        self.assertIn("ellrankinit", source)
        self.assertIn("bnfcertify", source)
        self.assertIn("ell2cover", source)
        self.assertIn("hyperellratpoints", source)
        self.assertNotIn('payload["exceptional_points"]', source)
        self.assertNotIn('payload["generic_points"]', source)

    def test_binary_linear_solver(self) -> None:
        vectors = ([1, 0, 1], [0, 1, 1])
        self.assertEqual(
            open_runner.f2_linear_combination_coefficients(vectors, [1, 1, 0]),
            [1, 1],
        )
        self.assertIsNone(
            open_runner.f2_linear_combination_coefficients(vectors, [0, 0, 1])
        )
        self.assertEqual(open_runner.f2_rank(vectors), 2)

    def test_standard_quotient_basis_extension(self) -> None:
        generic = ([1, 1, 0, 0], [0, 1, 1, 0])
        extension = open_runner.extend_rows_to_standard_basis(generic, 4)
        self.assertEqual(len(extension), 2)
        self.assertEqual(open_runner.f2_rank([*generic, *extension]), 4)
        self.assertTrue(all(sum(row) == 1 for row in extension))

    def test_control_quotient_classification(self) -> None:
        class Case:
            exceptional_points = ((0, 0),) * 4
            certified_rank_lower_bound = 21

        worker = {"total_two_selmer_dimension": 23, "two_torsion_dimension": 0}
        result = open_runner.quotient_classification(Case(), worker)
        self.assertEqual(result["relative_quotient_dimension"], 6)
        self.assertEqual(result["known_exceptional_quotient_dimension"], 4)
        self.assertEqual(result["unexplained_quotient_dimension"], 2)
        self.assertEqual(result["classes_not_realized_by_known_exceptional_subgroup"], 48)


if __name__ == "__main__":
    unittest.main()

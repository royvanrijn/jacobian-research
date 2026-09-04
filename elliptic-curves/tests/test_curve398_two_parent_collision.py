from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve398_two_parent_collision_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def multiply(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [
        [sum(a * b for a, b in zip(row, column)) for column in zip(*right)]
        for row in left
    ]


class Curve398TwoParentCollisionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.document = json.loads(ARTIFACT.read_text())

    def test_second_parent_compilation(self) -> None:
        second = self.document["second_parent"]
        self.assertEqual(second["priority_rank"], 63669)
        self.assertEqual(second["orbit_hex"], "0x06119")
        self.assertEqual(second["fibre_configuration"], "I2 at infinity + 22 I1")
        self.assertEqual(second["generic_mw16"]["complete_old_degree_one_section_count"], 180)
        self.assertEqual(second["generic_mw16"]["height_gram_determinant"], "474")
        self.assertTrue(second["generic_mw16"]["saturated"])
        self.assertTrue(second["parameter_recovery"]["isomorphic_to_curve398_over_Q"])

    def test_exact_collision_summary(self) -> None:
        collision = self.document["collision"]
        self.assertEqual(collision["rank_G1"], 16)
        self.assertEqual(collision["rank_G2"], 16)
        self.assertEqual(collision["rank_intersection"], 16)
        self.assertEqual(collision["rank_sum"], 16)
        self.assertTrue(collision["integral_subgroups_equal"])
        self.assertEqual(collision["smith_diagonal_nonzero"], [1] * 16)
        self.assertEqual(collision["quotient_free_rank"], 14)
        self.assertEqual(collision["quotient_torsion_invariant_factors_nontrivial"], [])
        self.assertEqual(collision["smith_index_in_public_M30"], "infinite")

    def test_unimodular_transition_replay(self) -> None:
        first = self.document["first_parent"]["matrix_16_by_30_rows"]
        second = self.document["second_parent"]["public_rank30_embedding"]["matrix_16_by_30_rows"]
        collision = self.document["collision"]
        second_in_first = collision["G2_basis_rows_in_G1_basis"]
        first_in_second = collision["G1_basis_rows_in_G2_basis"]
        self.assertEqual(multiply(second_in_first, first), second)
        self.assertEqual(multiply(first_in_second, second), first)
        self.assertEqual(multiply(second_in_first, first_in_second), [
            [int(row == column) for column in range(16)] for row in range(16)
        ])
        self.assertEqual(collision["basis_transition_determinants"], [1, 1])

    def test_hash_chain(self) -> None:
        for path, expected in self.document["inputs"].items():
            self.assertEqual(digest(ROOT / path), expected)
        generation = self.document["generation"]
        self.assertEqual(digest(ROOT / generation["script"]), generation["script_sha256"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path
import sys
import unittest

from sympy import Matrix, eye


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from audit_kumar_kuwata_f6_galois import (  # noqa: E402
    PRIMARY_SOURCE_BUNDLE_SHA256,
    PRIMARY_SOURCE_FILE_SHA256,
)


ARTIFACT = (
    Path(__file__).resolve().parents[2]
    / "artifacts"
    / "generated-results"
    / "elliptic_kumar_kuwata_f6_galois.json"
)
SCRIPT = CAS_DIRECTORY / "audit_kumar_kuwata_f6_galois.py"


class KumarKuwataF6GaloisAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_source_and_scope_are_pinned(self) -> None:
        source = self.report["source"]
        self.assertEqual(source["bundle_sha256"], PRIMARY_SOURCE_BUNDLE_SHA256)
        self.assertEqual(
            source["ancillary_file_sha256"], PRIMARY_SOURCE_FILE_SHA256
        )
        self.assertEqual(self.report["geometric_rank"], 18)
        self.assertEqual(self.report["exact_method"]["matching_fibres"], [2, 3, 5, 7])
        self.assertEqual(self.report["exact_method"]["p_minimal_vector_count"], 84)
        self.assertEqual(self.report["exact_method"]["q_minimal_vector_count"], 240)

    def test_reproduction_metadata_is_pinned(self) -> None:
        self.assertEqual(
            self.report["script_sha256"], sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(
            self.report["software"],
            {
                "python": "3.14.6",
                "sympy": "1.14.0",
                "pari_gp": "[2, 17, 4]",
            },
        )

    def test_action_matrices_satisfy_the_galois_relations(self) -> None:
        actions = self.report["generator_actions"]

        def full(name: str) -> Matrix:
            p = Matrix(actions[name]["p_matrix_rows"])
            q = Matrix(actions[name]["q_matrix_rows"])
            return Matrix.diag(p, q)

        conjugation = full("complex_conjugation")
        rotation = full("quartic_rotation")
        beta_flip = full("sqrt2_flip")
        identity = eye(18)
        self.assertEqual(conjugation**2, identity)
        self.assertEqual(rotation**4, identity)
        self.assertEqual(beta_flip**2, identity)
        self.assertEqual(conjugation * rotation * conjugation, rotation**-1)
        self.assertEqual(beta_flip * conjugation, conjugation * beta_flip)
        self.assertEqual(beta_flip * rotation, rotation * beta_flip)

    def test_fixed_and_quadratic_character_ranks(self) -> None:
        eigenspaces = self.report["quadratic_character_eigenspaces"]
        base = next(row for row in eigenspaces if row["base_curve_over_Q_t"])
        twists = [row for row in eigenspaces if not row["base_curve_over_Q_t"]]
        self.assertEqual(base["full_rank"], 5)
        self.assertEqual(base["p_subspace_rank"], 2)
        self.assertEqual(base["q_subspace_rank"], 3)
        self.assertEqual(sorted((row["full_rank"] for row in twists), reverse=True), [3, 2, 2, 1, 1, 0, 0])
        ranks_by_squareclass = {
            row["quadratic_character_squareclass"]: row["full_rank"]
            for row in twists
        }
        self.assertEqual(
            ranks_by_squareclass,
            {-1: 2, 2: 3, -2: 0, 3: 1, -3: 2, 6: 0, -6: 1},
        )
        conclusion = self.report["conclusion"]
        self.assertEqual(conclusion["rank_over_Q_t"], 5)
        self.assertEqual(conclusion["largest_quadratic_twist_rank_over_Q_t"], 3)
        self.assertFalse(conclusion["competitive_for_specialization_search"])
        self.assertFalse(conclusion["specialization_search_performed"])
        self.assertIn("finite-fibre", conclusion["mathematical_status"])


if __name__ == "__main__":
    unittest.main()

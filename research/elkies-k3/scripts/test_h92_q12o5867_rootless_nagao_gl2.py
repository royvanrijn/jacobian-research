#!/usr/bin/env python3
"""Tests for the q12/orbit5867 PGL2 Nagao chart wrapper."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "scan_h92_q12o5867_rootless_nagao_gl2.py"
SPEC = importlib.util.spec_from_file_location("q12o5867_gl2", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
GL2 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = GL2
SPEC.loader.exec_module(GL2)


class ProjectiveGL2ChartTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = GL2.load_family_model()

    def test_matrix_normalization_and_projective_transport(self) -> None:
        matrix = GL2.normalize_matrix((-6, -10, 2, 2))
        self.assertEqual(matrix, (3, 5, -1, -1))
        # u=(3a+5b)/(-a-b), evaluated at v=2/3.
        self.assertEqual(GL2.map_projective_pair(2, 3, matrix), (-21, 5, 21))

    def test_singular_matrix_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonzero determinant"):
            GL2.normalize_matrix((1, 2, 2, 4))

    def test_transformed_table_is_a_projective_permutation(self) -> None:
        blocks, rejected = GL2.build_residue_tables(self.model, ((19,),))
        self.assertEqual(rejected, ())
        matrix = (7, -11, 3, 5)
        transformed = GL2.transformed_blocks(blocks, matrix)[0][19]
        original = blocks[0][19]
        self.assertEqual(len(transformed), 20)
        # Compare identities rather than LocalSymbol.projective_index, which
        # intentionally retains its original-table provenance.
        self.assertEqual(
            sorted(symbol.projective_index for symbol in transformed),
            sorted(symbol.projective_index for symbol in original),
        )

    def test_noninvertible_discovery_prime_is_rejected(self) -> None:
        blocks, _ = GL2.build_residue_tables(self.model, ((19,),))
        with self.assertRaisesRegex(ValueError, "p=19"):
            GL2.transformed_blocks(blocks, (19, 0, 0, 1))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Tests for bounded norm-one cover direction enumeration."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/generate_q12o5867_norm_one_cover_candidates.py"
SPEC = importlib.util.spec_from_file_location("q12_norm_one", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Q12O5867NormOneCandidateTests(unittest.TestCase):
    def test_projective_direction_normalization(self) -> None:
        self.assertIsNone(MODULE.normalize_direction((0, 0, 0)))
        self.assertEqual(MODULE.normalize_direction((-2, 4, -6)), (1, -2, 3))
        self.assertEqual(MODULE.normalize_direction((0, -3, 6)), (0, 1, -2))
        self.assertEqual(MODULE.normalize_direction((0, 0, -5)), (0, 0, 1))


if __name__ == "__main__":
    unittest.main()

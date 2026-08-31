#!/usr/bin/env python3
"""Tests for the q12o5867 normalized slope-quartic search helpers."""

from __future__ import annotations

from fractions import Fraction
import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/scripts/search_q12o5867_section_slope_slices.py"
SPEC = importlib.util.spec_from_file_location("q12o5867_slope_slices", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SLOPES = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SLOPES
SPEC.loader.exec_module(SLOPES)


class Q12O5867SlopeSliceTests(unittest.TestCase):
    def test_ratpoints_abscissae_recover_both_exact_ordinate_signs(self) -> None:
        points = SLOPES.quartic_points_from_abscissae(
            (Fraction(1), Fraction(2), Fraction(1)),
            (Fraction(0), Fraction(2)),
        )
        self.assertEqual(
            points,
            (
                (Fraction(0), Fraction(1)),
                (Fraction(0), Fraction(-1)),
                (Fraction(2), Fraction(3)),
                (Fraction(2), Fraction(-3)),
            ),
        )

    def test_ratpoints_abscissa_is_rechecked_exactly(self) -> None:
        with self.assertRaises(AssertionError):
            SLOPES.quartic_points_from_abscissae(
                (Fraction(1), Fraction(1)), (Fraction(1),)
            )

    def test_exact_signed_pair_relations_are_recognized(self) -> None:
        model = tuple(map(Fraction, (0, 0, 0, 0, -2)))
        point = (Fraction(3), Fraction(5))
        doubled = SLOPES.add_rational_points(model, point, point)
        self.assertIsNotNone(doubled)
        lookup = SLOPES.signed_pair_relation_lookup(model, (point,))
        relation = lookup[SLOPES.canonical_short_point(doubled)]
        self.assertEqual(relation["left_section_index"], 0)
        self.assertEqual(relation["right_section_index"], 0)


if __name__ == "__main__":
    unittest.main()

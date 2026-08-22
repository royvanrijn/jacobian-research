#!/usr/bin/env python3
"""Focused replay for the bounded transverse-seed census."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest
from fractions import Fraction
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
SCRIPT = CAS / "screen_mestre_two_section_transverse_seeds.py"


def load() -> object:
    sys.path.insert(0, str(CAS))
    spec = importlib.util.spec_from_file_location("mestre_transverse_seed_screen", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class MestreTwoSectionSeedScreenTest(unittest.TestCase):
    def test_all_prime_union_does_not_drop_a_later_smooth_reduction(self) -> None:
        screen = load()
        calls = []

        def at_prime(moduli, prime, precision):
            calls.append(prime)
            if prime == 7:
                return ((Fraction(0), Fraction(1)),)
            if prime == 11:
                return ((Fraction(2), Fraction(3)),)
            raise ValueError("bad reduction")

        with patch.object(screen, "sections_at_prime", side_effect=at_prime):
            sections = screen.sections_across_primes(
                (Fraction(1),) * 4, (7, 11, 13), 8
            )
        self.assertEqual(
            sections,
            ((Fraction(0), Fraction(1)), (Fraction(2), Fraction(3))),
        )
        self.assertEqual(calls, [7, 11, 13])

    def test_rank_seven_seed_at_known_census_offset(self) -> None:
        screen = load()
        result = screen.screen(
            screen.DEFAULT_INPUT, screen.DEFAULT_PRIMES, 8, start=45, count=1
        )
        self.assertEqual(result["candidate_count"], 167)
        self.assertEqual(result["records_with_reconstructed_sections"], 1)
        self.assertEqual(result["exact_transverse_pair_count"], 4)
        record = result["records"][0]
        self.assertEqual(record["roots"], [0, 7, 79, 81, 128, 137])
        self.assertEqual(len(record["reconstructed_sections"]), 5)
        self.assertEqual(len(record["transverse_rank_seven_pairs"]), 4)


if __name__ == "__main__":
    unittest.main()

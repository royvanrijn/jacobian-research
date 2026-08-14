#!/usr/bin/env python3
"""Focused exact tests for the new family-2 companion section."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))

from verify_mestre_02595143168205_discriminants import (  # noqa: E402
    replay as replay_discriminants,
)
from verify_mestre_02595143168205_rank13_section import (  # noqa: E402
    replay as replay_rank13,
)


class Mestre02595143168205Rank13Test(unittest.TestCase):
    def test_nonvisible_section_certifies_generic_rank_13(self) -> None:
        result = replay_rank13()
        self.assertEqual(result["family_roots"], [0, 25, 95, 143, 168, 205])
        self.assertEqual(result["T"], "337/394")
        self.assertEqual(result["new_section_original_x"], "(-2375+37*T)/23")
        self.assertEqual(
            result["generic_companion_identity"]["identity_verified_over"],
            "Q[T]",
        )
        self.assertEqual(result["visible_plus_infinity_rank_lower_bound"], 12)
        self.assertEqual(result["one_companion_rank_lower_bound"], 13)
        self.assertEqual(result["all_six_companions_rank_lower_bound"], 13)
        self.assertEqual(
            result["one_companion_pivots"],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14],
        )

    def test_discriminant_cores_are_irreducible_degrees_20_and_40(self) -> None:
        result = replay_discriminants()
        self.assertEqual(result["base_change"], "T=(39146-u^2)/(2u)")
        self.assertEqual(result["direct"]["factor_degrees_over_Q"], [20])
        self.assertEqual(result["pullback"]["factor_degrees_over_Q"], [40])
        self.assertTrue(result["direct"]["squarefree_over_Q"])
        self.assertTrue(result["pullback"]["squarefree_over_Q"])
        self.assertEqual(
            result["direct"]["coefficient_sha256"],
            "fc36f00ad71a6b30126402aae310cdd2c9d35553e9f22910334c2ba4b9a05590",
        )
        self.assertEqual(
            result["pullback"]["coefficient_sha256"],
            "876a5e46a21c20cf531eb63469b55fe2cecf58d4fd5fdedfedacb3950a0e3a41",
        )


if __name__ == "__main__":
    unittest.main()

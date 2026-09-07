#!/usr/bin/env python3
"""Focused exact tests for the four split-infinity Mestre families."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))

from search_mestre_dsquare_four import FAMILIES, family_geometry  # noqa: E402
from verify_mestre_dsquare_four_u197 import PINNED_ARTIFACT, replay  # noqa: E402


ARTIFACT_ROOT = ROOT / "artifacts/local/elliptic-curves/mestre-dsquare-four-v1"


class MestreDSquareFourTest(unittest.TestCase):
    def test_family_geometry_and_split_infinity_baselines(self) -> None:
        expected = (
            ((0, 7, 225, 232, 235, 265), "84878/3", 3, "a4bcf120aa4c6f7bf7cf0121cd3d6403413fc822d3e8c2fd7296fc69e06782b4"),
            ((0, 9, 213, 247, 256, 291), "196250/3", 3, "de45534f1e3bdc15f8be29429a9cbdb7c1d645d8ec81446cf1c6973a415dde41"),
            ((0, 25, 95, 143, 168, 205), "39146", 1, "c0adcd48da63f3e949dbcd9df1d70c9c8f2b7fbb2510acb917defb2cb2c05a91"),
            ((0, 43, 128, 197, 231, 289), "55950", 3, "f02dc4c28fa73e0a9543286ad03560fdff2f10f3f0f100deb5a8be0c527d87f5"),
        )
        for family, (roots, constant, multiplier, digest) in zip(FAMILIES, expected):
            self.assertEqual(family.roots, roots)
            self.assertEqual(str(family.base_constant), constant)
            self.assertEqual(family.leading_square_multiplier, multiplier)
            geometry = family_geometry(family)
            self.assertEqual(geometry["quartic_condition"], "0")
            self.assertEqual(geometry["primitive_discriminant_degree_in_T"], 20)
            self.assertTrue(geometry["primitive_discriminant_even"])
            self.assertEqual(geometry["u1_visible_plus_infinity_point_count"], 13)
            self.assertEqual(geometry["u1_exact_mod3_rank_lower_bound"], 12)
            self.assertEqual(geometry["u1_point_sha256"], digest)

    def test_u197_promotion_replays_exactly(self) -> None:
        result = replay(PINNED_ARTIFACT, verify_pari=False, pari_timeout=0)
        self.assertEqual(result["T"], "337/394")
        self.assertEqual(result["generic_rank_lower_bound_after_base_change"], 13)
        self.assertEqual(result["generic_companion_identity_verified_over"], "Q[T]")
        self.assertEqual(
            result["primitive_discriminant_degrees"],
            {"direct_T": 20, "base_changed_u": 40},
        )
        self.assertEqual(result["exact_rank_lower_bound"], 17)
        self.assertEqual(result["root_number"], -1)
        self.assertEqual(
            result["conductor"],
            "2462086522751621334987931952469307556796057284118717977320345864383117775914",
        )

    @unittest.skipUnless(
        (ARTIFACT_ROOT / "point-certificates/f2_u197_1.json").exists(),
        "the discovery outputs are intentionally local",
    )
    def test_u197_discovery_outputs_match_pinned_subset(self) -> None:
        result = replay(
            PINNED_ARTIFACT,
            verify_pari=False,
            pari_timeout=0,
            discovery_root=ARTIFACT_ROOT,
        )
        self.assertTrue(result["discovery_raw_audited"])

if __name__ == "__main__":
    unittest.main()

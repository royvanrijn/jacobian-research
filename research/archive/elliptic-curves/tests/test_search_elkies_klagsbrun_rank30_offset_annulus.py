#!/usr/bin/env python3
"""Regression tests for the public-center rank-30 offset annulus."""

from __future__ import annotations

import hashlib
import json
from math import gcd
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from search_elkies_klagsbrun_rank30_offset_annulus import (  # noqa: E402
    EXPECTED_PRIOR_BASE_SHA256,
    EXPECTED_SOURCE_ENGINE_SHA256,
    PRIOR_BASE_ARTIFACT,
    SOURCE_ENGINE,
    annulus_mask,
    exact_center_separation,
    positive_coprime_annulus_count,
)


SCRIPT = CAS / "search_elkies_klagsbrun_rank30_offset_annulus.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_offset_annulus.json"
)
EXPECTED_PROCESSED = 80_474_692_618
EXPECTED_SURVIVORS = 269_139
EXPECTED_SURVIVOR_SHA256 = "8bf8fbe73b7fdf201a43e67c03b0857e12c506df73f8b0d918b2330ebd94da49"


class Rank30OffsetAnnulusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_annulus_mask_and_coprime_count(self) -> None:
        minimum, maximum = 3, 9
        mask = annulus_mask(minimum, maximum)
        values = [
            index - maximum
            for index in range(2 * maximum + 1)
            if (mask >> index) & 1
        ]
        self.assertEqual(values, list(range(-maximum, -minimum + 1)) + list(range(minimum, maximum + 1)))
        for denominator in range(1, 20):
            expected = sum(gcd(value, denominator) == 1 for value in range(minimum, maximum + 1))
            self.assertEqual(
                positive_coprime_annulus_count(minimum, maximum, denominator), expected
            )

    def test_pinned_sources_and_exact_nonoverlap(self) -> None:
        self.assertEqual(hashlib.sha256(SOURCE_ENGINE.read_bytes()).hexdigest(), EXPECTED_SOURCE_ENGINE_SHA256)
        self.assertEqual(hashlib.sha256(PRIOR_BASE_ARTIFACT.read_bytes()).hexdigest(), EXPECTED_PRIOR_BASE_SHA256)
        separation = exact_center_separation(3163, 65536)
        self.assertTrue(separation["pairwise_center_boxes_disjoint"])
        self.assertTrue(self.data["nonoverlap"]["disjoint_from_prior_base_by_offset_magnitude"])
        self.assertTrue(self.data["nonoverlap"]["disjoint_from_prior_extension_by_denominator"])

    def test_complete_exact_negative_search(self) -> None:
        result = self.data["search_result"]
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertEqual(result["processed_primitive_candidate_count"], EXPECTED_PROCESSED)
        self.assertEqual(result["declared_primitive_candidate_count"], EXPECTED_PROCESSED)
        self.assertEqual(result["modular_survivor_count_after_primitivity"], EXPECTED_SURVIVORS)
        self.assertEqual(result["modular_survivor_manifest_sha256"], EXPECTED_SURVIVOR_SHA256)
        self.assertEqual(result["exact_nonsquare_count_after_sieve"], EXPECTED_SURVIVORS)
        self.assertEqual(result["exact_square_abscissa_count"], 0)
        self.assertFalse(result["rank30_target_hit"])

    def test_artifact_script_hash(self) -> None:
        self.assertFalse(self.data["target_hit"])
        self.assertEqual(
            self.data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()

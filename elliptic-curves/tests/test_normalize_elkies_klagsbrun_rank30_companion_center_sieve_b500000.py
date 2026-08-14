#!/usr/bin/env python3
"""Regression tests for the extended companion-center rank-30 sieve."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from normalize_elkies_klagsbrun_rank30_companion_center_sieve_b500000 import (  # noqa: E402
    EXPECTED_CENTER_COUNT,
    EXPECTED_CENTER_SHA256,
    EXPECTED_ENGINE_SHA256,
    EXPECTED_PROCESSED_COUNT,
    EXPECTED_RAW_SHA256,
    EXPECTED_SURVIVOR_COUNT,
    EXPECTED_SURVIVOR_SHA256,
    validate,
)


SCRIPT = CAS / "normalize_elkies_klagsbrun_rank30_companion_center_sieve_b500000.py"
ENGINE = CAS / "search_elkies_klagsbrun_rank30_companion_center_sieve.py"
GENERATED = ROOT / "artifacts" / "generated-results"
RAW = GENERATED / "elliptic_elkies_klagsbrun_rank30_companion_center_sieve_b500000_raw.json"
ARTIFACT = GENERATED / "elliptic_elkies_klagsbrun_rank30_companion_center_sieve_b500000.json"


class Rank30CompanionCenterB500000Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.raw = json.loads(RAW.read_text(encoding="utf-8"))
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_pinned_inputs_and_validator(self) -> None:
        self.assertEqual(hashlib.sha256(ENGINE.read_bytes()).hexdigest(), EXPECTED_ENGINE_SHA256)
        self.assertEqual(hashlib.sha256(RAW.read_bytes()).hexdigest(), EXPECTED_RAW_SHA256)
        validate(self.raw)

    def test_center_manifest_and_complete_region(self) -> None:
        centers = self.data["center_manifest"]
        self.assertEqual(centers["count"], EXPECTED_CENTER_COUNT)
        self.assertEqual(centers["sha256"], EXPECTED_CENTER_SHA256)
        self.assertEqual(self.data["parameters"]["denominator_interval"], [50001, 500000])
        result = self.data["search_result"]
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertEqual(result["processed_primitive_candidate_count"], EXPECTED_PROCESSED_COUNT)

    def test_every_survivor_is_an_exact_nonsquare(self) -> None:
        result = self.data["search_result"]
        self.assertEqual(result["modular_survivor_count_after_primitivity"], EXPECTED_SURVIVOR_COUNT)
        self.assertEqual(result["modular_survivor_manifest_sha256"], EXPECTED_SURVIVOR_SHA256)
        self.assertEqual(result["exact_nonsquare_count_after_sieve"], EXPECTED_SURVIVOR_COUNT)
        self.assertEqual(result["exact_square_abscissa_count"], 0)
        self.assertFalse(result["rank30_target_hit"])

    def test_normalized_provenance(self) -> None:
        reproduction = self.data["reproduction"]
        self.assertEqual(reproduction["raw_execution_artifact_sha256"], EXPECTED_RAW_SHA256)
        self.assertEqual(
            reproduction["normalizer_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertFalse(reproduction["metadata_normalization_reran_bounded_search"])
        self.assertIn("--denominator-max 500000", reproduction["bounded_execution_command"])


if __name__ == "__main__":
    unittest.main()

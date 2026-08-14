#!/usr/bin/env python3
"""Checks for the conditional Delta=11/5 closure of the new rank-15 fiber."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest

import mpmath


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "explicit_formula_mestre_02136217261290_t2_delta22.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_mestre_02136217261290_t2_explicit_formula_delta22.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "27965afe7551bb35ec231428b3cf4355445fbcd4893bd75a7b0de190e000934f"
)
EXPECTED_ARTIFACT_SHA256 = (
    "6092f3b547e53275cd16842a595488ba97230f299fafae2305d75f92239b070a"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "explicit_formula_mestre_02136217261290_t2_delta22", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Mestre02136217261290Rank15ExplicitFormulaTests(unittest.TestCase):
    def test_pinned_inputs_and_support(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(MODULE.prime_limit(Fraction(11, 5)), 1_007_525)
        self.assertEqual(
            MODULE.file_sha256(MODULE.CERTIFICATE),
            MODULE.EXPECTED_CERTIFICATE_SHA256,
        )
        self.assertEqual(
            MODULE.file_sha256(MODULE.SOURCE_ENGINE),
            MODULE.EXPECTED_SOURCE_ENGINE_SHA256,
        )

    def test_generated_artifact_closes_only_conditionally(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["status"], "conditional explicit-formula diagnostic complete"
        )
        self.assertFalse(data["target_hit"])
        self.assertEqual(data["curve"]["roots"], [0, 2, 136, 217, 261, 290])
        self.assertEqual(data["curve"]["parameter"], "2")
        self.assertEqual(data["curve"]["unconditional_algebraic_rank_lower_bound"], 15)
        self.assertEqual(data["curve"]["root_number"], -1)
        self.assertEqual(data["explicit_formula"]["delta"], "11/5")
        self.assertEqual(data["explicit_formula"]["support_prime_limit"], 1_007_525)
        upper = mpmath.mpf(
            data["explicit_formula"]["conservative_explicit_formula_upper"]
        )
        self.assertGreater(upper, 15)
        self.assertLess(upper, 17)
        self.assertTrue(
            data["explicit_formula"]["conservative_upper_strictly_below_17"]
        )
        self.assertIn("no unconditional rank upper bound", data["conclusion"]["unconditional"])
        self.assertEqual(
            data["conclusion"]["under_bsd_and_grh"],
            "The algebraic rank is exactly 15.",
        )

    def test_artifact_self_pins_the_single_capped_call(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["provenance"]["script_sha256"], sha256(SCRIPT))
        self.assertEqual(
            data["provenance"]["rank15_certificate_sha256"],
            sha256(MODULE.CERTIFICATE),
        )
        self.assertTrue(data["provenance"]["single_foreground_pari_call"])
        self.assertEqual(data["provenance"]["same_stage_retries"], 0)
        self.assertEqual(data["provenance"]["owned_processes_remaining"], 0)


if __name__ == "__main__":
    unittest.main()

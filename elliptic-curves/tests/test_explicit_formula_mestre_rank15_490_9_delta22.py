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
SCRIPT = CAS / "explicit_formula_mestre_rank15_490_9_delta22.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_mestre_rank15_490_9_explicit_formula_delta22.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "explicit_formula_mestre_rank15_490_9_delta22", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MestreRank15Delta22ExplicitFormulaTests(unittest.TestCase):
    def test_pinned_inputs_and_support(self) -> None:
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
        self.assertEqual(data["status"], "conditional_fixed_fiber_rank15_closure")
        self.assertFalse(data["target_hit"])
        self.assertEqual(
            data["curve"]["unconditional_algebraic_rank_lower_bound"], 15
        )
        self.assertEqual(data["curve"]["root_number"], -1)
        self.assertEqual(data["explicit_formula"]["delta"], "11/5")
        self.assertEqual(
            data["explicit_formula"]["support_prime_limit"], 1_007_525
        )
        upper = mpmath.mpf(
            data["explicit_formula"]["conservative_explicit_formula_upper"]
        )
        self.assertGreater(upper, 16)
        self.assertLess(upper, 17)
        self.assertIn("no rank upper bound", data["conclusion"]["unconditional"])
        self.assertEqual(
            data["conclusion"]["under_bsd_and_grh"],
            "The algebraic rank is exactly 15.",
        )

    def test_artifact_is_self_pinned(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["provenance"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["provenance"]["rank15_certificate_sha256"],
            hashlib.sha256(MODULE.CERTIFICATE.read_bytes()).hexdigest(),
        )
        self.assertTrue(data["provenance"]["single_foreground_pari_call"])
        self.assertTrue(data["provenance"]["no_retry"])


if __name__ == "__main__":
    unittest.main()

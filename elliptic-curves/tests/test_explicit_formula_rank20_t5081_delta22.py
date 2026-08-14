from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import unittest

import mpmath


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "explicit_formula_rank20_t5081_delta22.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank20_t5081_explicit_formula_delta22.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "explicit_formula_rank20_t5081_delta22", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Delta22ExplicitFormulaTests(unittest.TestCase):
    def test_support_and_gp_program(self) -> None:
        self.assertEqual(MODULE.prime_limit(Fraction(11, 5)), 1_007_525)
        program = MODULE.gp_program(
            MODULE.EXPECTED_MINIMAL_MODEL,
            delta=Fraction(11, 5),
            support_limit=1_007_525,
        )
        self.assertIn("D=11/5", program)
        self.assertIn("LIM=1007525", program)
        self.assertNotIn("ell2cover", program)
        self.assertNotIn("ellrank", program)

    def test_tail_bound_is_conservative_and_small(self) -> None:
        with mpmath.workdps(40):
            bound = MODULE.archimedean_tail_absolute_bound(Fraction(11, 5))
            self.assertGreater(bound, 0)
            self.assertLess(bound, mpmath.mpf("0.01"))

    def test_pinned_inputs(self) -> None:
        self.assertEqual(
            MODULE.sha256_file(MODULE.RANK20_CERTIFICATE),
            MODULE.RANK20_CERTIFICATE_SHA256,
        )
        self.assertEqual(
            MODULE.sha256_file(MODULE.DELTA2_DIAGNOSTIC),
            MODULE.DELTA2_DIAGNOSTIC_SHA256,
        )

    def test_generated_artifact(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["method"]["delta"], "11/5")
        self.assertEqual(data["method"]["prime_limit"], 1_007_525)
        self.assertTrue(data["direct_delta2_calibration"]["passes_published_value"])
        self.assertTrue(data["delta_11_over_5"]["strictly_less_than_22"])
        self.assertLess(
            mpmath.mpf(
                data["delta_11_over_5"]["conservative_explicit_formula_upper"]
            ),
            22,
        )


if __name__ == "__main__":
    unittest.main()

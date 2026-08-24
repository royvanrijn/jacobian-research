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
SCRIPT = CAS / "explicit_formula_fermigier_rank20_28917_20_delta22.py"
ARTIFACT = (
    ROOT
    / "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank20_28917_20_explicit_formula_delta22.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "explicit_formula_fermigier_rank20_28917_20_delta22", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FermigierRank20Delta22Tests(unittest.TestCase):
    def test_support_and_curve_inputs(self) -> None:
        self.assertEqual(MODULE.prime_limit(Fraction(11, 5)), 1_007_525)
        self.assertEqual(MODULE.EXPECTED_PARAMETER, "28917/20")
        self.assertEqual(MODULE.EXPECTED_MINIMAL_MODEL[:3], [1, 1, 1])
        self.assertEqual(
            MODULE.sha256_file(MODULE.IMPORTED_CERTIFICATE),
            MODULE.IMPORTED_CERTIFICATE_SHA256,
        )

    def test_imported_rank_checkpoint(self) -> None:
        certificate, _ = MODULE.load_inputs()
        self.assertEqual(certificate["point_cloud"]["selected_count"], 20)
        self.assertEqual(certificate["global_curve"]["root_number"], 1)
        self.assertTrue(certificate["global_curve"]["below_strict_log_target"])

    def test_generated_artifact(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["method"]["delta"], "11/5")
        self.assertEqual(data["method"]["prime_limit"], 1_007_525)
        self.assertEqual(
            data["delta_11_over_5"]["prime_sum"]["prime_count"], 79_057
        )
        self.assertEqual(
            data["delta_11_over_5"]["prime_sum"]["prime_power_term_count"],
            79_293,
        )
        self.assertTrue(data["delta_11_over_5"]["strictly_less_than_22"])
        self.assertLess(
            mpmath.mpf(
                data["delta_11_over_5"]["conservative_explicit_formula_upper"]
            ),
            22,
        )
        self.assertEqual(
            data["script_sha256"], MODULE.sha256_file(SCRIPT),
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_rank21_historical_finalists.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank21_historical_finalists.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_rank21_historical_finalists", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class HistoricalFinalistTests(unittest.TestCase):
    def test_source_parameters_and_factor_two_convention(self) -> None:
        self.assertEqual(
            MODULE.PAPER_PARAMETERS,
            (
                Fraction(1393, 216),
                Fraction(1649, 12),
                Fraction(6629, 348),
                Fraction(8057, 876),
            ),
        )
        self.assertEqual(
            MODULE.CONSTRUCTOR_PARAMETERS,
            (
                Fraction(1393, 108),
                Fraction(1649, 6),
                Fraction(6629, 174),
                Fraction(8057, 438),
            ),
        )

    def test_published_calibration_is_pinned(self) -> None:
        data = MODULE.validate_source_convention()
        self.assertEqual(
            MODULE.sha256_file(MODULE.CALIBRATION_ARTIFACT),
            MODULE.CALIBRATION_ARTIFACT_SHA256,
        )
        self.assertEqual(
            data["published_record_calibration"]["constructor_parameter"],
            "14721/188",
        )

    def test_generated_artifact(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(len(data["historical_finalists"]), 4)
        self.assertEqual(data["target"]["rank_at_least"], 21)
        self.assertEqual(data["bounds_and_caveats"]["stage_heights"], [50_000, 250_000, 1_000_000])
        for checkpoint in data["exact_checkpoints_stable_numerical_rank_at_least_18"]:
            self.assertGreaterEqual(checkpoint["stable_numerical_rank"], 18)
            self.assertEqual(
                checkpoint["exact_rank_certificate"]["status"], "certified"
            )


if __name__ == "__main__":
    unittest.main()

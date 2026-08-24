from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "explicit_formula_rank18_t6629_delta22.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_nagao_rank18_t6629_explicit_formula_delta22.json"
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location("explicit_formula_rank18_t6629_delta22", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class T6629ExplicitFormulaTests(unittest.TestCase):
    def test_pinned_input_and_support(self) -> None:
        checkpoint, model = MODULE.load_checkpoint()
        self.assertEqual(MODULE.PARAMETER_T, Fraction(6629, 174))
        self.assertEqual(checkpoint["exact_rank_certificate"]["certified_algebraic_rank_lower_bound"], 18)
        self.assertEqual(len(model), 5)
        self.assertEqual(MODULE.prime_limit(MODULE.DELTA), 1_007_525)

    def test_generated_artifact(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["candidate"]["constructor_parameter"], "6629/174")
        self.assertEqual(data["candidate"]["root_number"], 1)
        self.assertEqual(data["prime_power_sum"]["prime_count"], 79_057)
        self.assertEqual(data["prime_power_sum"]["prime_power_term_count"], 79_293)
        self.assertLess(
            Decimal(data["explicit_formula"]["conservative_upper_value"]),
            Decimal(20),
        )
        self.assertTrue(data["explicit_formula"]["conservative_upper_value_less_than_20"])


if __name__ == "__main__":
    unittest.main()

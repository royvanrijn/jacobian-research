from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from explicit_formula_rank20_t5081 import (  # noqa: E402
    E20_MODEL,
    PRIME_LIMIT,
    T5081_MODEL,
    calibrated_bound,
)
from explicit_formula_rank_diagnostic import DELTA, gp_program  # noqa: E402


SCRIPT = CAS / "explicit_formula_rank20_t5081.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank20_t5081_explicit_formula.json"
)


class ExplicitFormulaRank20T5081Tests(unittest.TestCase):
    def test_gp_program_has_both_exact_models(self) -> None:
        program = gp_program((("E20", E20_MODEL), ("T5081", T5081_MODEL)))
        self.assertIn(f"D={DELTA};", program)
        self.assertIn(f"LIM={PRIME_LIMIT};", program)
        self.assertIn("ROW|E20|", program)
        self.assertIn("ROW|T5081|", program)

    def test_pinned_relative_value(self) -> None:
        records = {
            "E20": {
                "log_conductor": "170.0876648422488811354314776110703548643",
                "prime_sum": "-53.5271982465901051823721048180945488965",
            },
            "T5081": {
                "log_conductor": "174.2498162285480383539049736903482307895",
                "prime_sum": "-54.6685728786628337",
            },
        }
        result = calibrated_bound(records)
        self.assertGreater(float(result["calibrated_upper_value"]), 22.20)
        self.assertLess(float(result["calibrated_upper_value"]), 22.22)

    def test_generated_artifact_matches_script_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the explicit-formula artifact has not been generated")
        artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            artifact["script_sha256"], hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertTrue(artifact["comparison"]["calibrated_upper_value_less_than_23"])
        self.assertEqual(artifact["comparison"]["candidate_root_number"], 1)
        self.assertIn("no algebraic-rank upper bound", artifact["interpretation"]["unconditional"])


if __name__ == "__main__":
    unittest.main()

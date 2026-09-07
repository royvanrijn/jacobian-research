from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_rank21_t6629_alternate_covers.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_nagao_rank21_t6629_alternate_covers.json"
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location("search_nagao_rank21_t6629_alternate_covers", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class T6629AlternateCoverTests(unittest.TestCase):
    def test_pinned_rank18_basis(self) -> None:
        path = ROOT / "artifacts/generated-results/elliptic_nagao_rank21_historical_finalists.json"
        coefficients, basis, certificate, conductor = MODULE.load_exact_basis(path)
        self.assertEqual(MODULE.PARAMETER_T, Fraction(6629, 174))
        self.assertEqual(len(coefficients), 5)
        self.assertEqual(len(basis), 18)
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 18)
        self.assertEqual(conductor["root_number"], 1)

    def test_declared_cover_budget(self) -> None:
        self.assertEqual((1 << MODULE.INPUT_CERTIFIED_RANK) - 1, 262_143)
        self.assertEqual(MODULE.SELECTED_COVER_COUNT * MODULE.CHARTS_PER_COVER, 120)

    def test_generated_artifact(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["candidate"]["constructor_parameter"], "6629/174")
        self.assertEqual(
            data["declared_budget"]["all_nonzero_certified_mod2_classes_identity_scored"],
            262_143,
        )
        self.assertEqual(data["declared_budget"]["pilot_chart_count"], 120)
        self.assertEqual(
            data["declared_budget"]["completed_search_calls"]
            + data["declared_budget"]["failed_or_timed_out_search_calls"],
            140,
        )
        self.assertEqual(
            data["results"]["combined_finite_reduction_rank"],
            data["results"]["certified_rank_lower_bound_after_search"],
        )


if __name__ == "__main__":
    unittest.main()

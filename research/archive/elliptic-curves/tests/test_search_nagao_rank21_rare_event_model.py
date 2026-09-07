from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_rank21_rare_event_model.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_rare_event_model.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_rank21_rare_event_model", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class RareEventRank21ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tables = MODULE.build_residue_tables(MODULE.MODEL_CUTOFF)

    def test_pinned_labels_preserve_theorem_status(self) -> None:
        audit = MODULE.pinned_label_audit(ROOT)
        self.assertEqual(
            [record["certified_algebraic_rank_lower_bound"] for record in audit["positives"]],
            [19, 18, 18, 18],
        )
        self.assertEqual(len(audit["controls"]), 13)
        self.assertTrue(
            all(
                "not an algebraic-rank upper bound" in record["scope_warning"]
                for record in audit["controls"]
            )
        )
        self.assertEqual(
            audit["published_calibration"]["training_role"],
            "held_out_calibration_only",
        )

    def test_local_trace_fingerprint_replays_exactly(self) -> None:
        observations = {
            item.prime: item
            for item in MODULE.local_observations(Q(6793, 64), self.tables)
        }
        expected = {
            11: (8, -5),
            19: (15, -6),
            31: (2, -9),
            41: (3, -11),
            47: (7, -12),
            59: (37, -14),
        }
        for prime, (projective_index, trace) in expected.items():
            self.assertEqual(observations[prime].projective_index, projective_index)
            self.assertEqual(observations[prime].trace, trace)
            self.assertTrue(observations[prime].good_reduction)
        self.assertEqual(
            [
                prime
                for prime, item in observations.items()
                if not item.good_reduction
            ],
            [5, 7, 13, 23, 37, 83],
        )

    def test_cross_validation_rejects_before_expensive_search(self) -> None:
        result = MODULE.build_result(ROOT)
        self.assertEqual(result["status"], "model_rejected_no_candidate_scan")
        self.assertFalse(result["selection_decision"]["scan_authorized"])
        self.assertEqual(result["selection_decision"]["accepted_models"], [])
        primary = result["models"][0]["cross_validation"]
        self.assertFalse(primary["passes_all_predeclared_thresholds"])
        self.assertEqual(primary["minimum_positive_percentile"]["exact"], "6/13")
        self.assertEqual(
            result["expensive_search"],
            {
                "broad_rational_population_scanned": 0,
                "conductors_computed": 0,
                "point_searches_launched": 0,
                "finite_reduction_certificates_triggered": 0,
            },
        )

    def test_generated_artifact_is_pinned_to_current_script(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["script_sha256"], MODULE.file_sha256(SCRIPT))
        self.assertFalse(data["selection_decision"]["scan_authorized"])
        self.assertEqual(data["local_population"]["curve_count"], 18)
        self.assertEqual(data["local_population"]["prime_count"], 44)


if __name__ == "__main__":
    unittest.main()

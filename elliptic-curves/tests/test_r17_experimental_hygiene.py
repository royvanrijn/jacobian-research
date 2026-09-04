from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
HALF_ANALYZER = (
    ROOT
    / "elkies-k3/scripts/"
    "analyze_r17_prospective_crt_half_lattice_censoring_gated.py"
)
GROUP_BUILDER = (
    ROOT / "elliptic-curves/scripts/audit_r17_training_arithmetic_groups.py"
)
GROUP_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "r17_training_arithmetic_group_audit_v1.json"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def outcome(*, eligible: bool, status: str, gain: int = 0) -> dict:
    return {
        "analysis_eligible_complete_stage_a": eligible,
        "status": status,
        "stage_a": {"certified_quotient_gain": gain},
        "largest_certified_rank_lower_bound": 17 + gain,
    }


class R17ExperimentalHygieneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.analyzer = load_module("r17_half_analyzer", HALF_ANALYZER)
        cls.group_builder = load_module("r17_group_builder", GROUP_BUILDER)
        cls.group_artifact = json.loads(GROUP_ARTIFACT.read_text())

    def test_unbalanced_censoring_suppresses_all_primary_effect_estimates(self) -> None:
        exposed = [
            outcome(eligible=True, status="COMPLETE", gain=1),
            outcome(eligible=False, status="CENSORED_TIMEOUT"),
        ]
        control = [
            outcome(eligible=True, status="COMPLETE"),
            outcome(eligible=True, status="COMPLETE"),
        ]
        comparison = self.analyzer.comparison_record("test", exposed, control)
        self.assertFalse(comparison["prospective_conclusion_authorized"])
        self.assertIsNone(comparison["risk_difference"])
        self.assertIsNone(comparison["risk_ratio"])
        self.assertIsNone(comparison["odds_ratio"])
        self.assertIsNone(comparison["fisher_exact_two_sided_p"])
        self.assertTrue(
            comparison["exposed"][
                "scheduled_denominator_is_not_an_inferential_event_rate"
            ]
        )

    def test_balanced_censoring_uses_only_complete_rows(self) -> None:
        exposed = [
            outcome(eligible=True, status="COMPLETE", gain=1),
            outcome(eligible=False, status="CENSORED_TIMEOUT"),
        ]
        control = [
            outcome(eligible=True, status="COMPLETE"),
            outcome(eligible=False, status="CENSORED_TIMEOUT"),
        ]
        comparison = self.analyzer.comparison_record("test", exposed, control)
        self.assertTrue(comparison["prospective_conclusion_authorized"])
        self.assertEqual(
            comparison["gated_complete_case_counts"],
            {
                "exposed_events": 1,
                "exposed_total": 1,
                "control_events": 0,
                "control_total": 1,
            },
        )
        self.assertEqual(comparison["risk_difference"], 1.0)

    def test_different_censor_reasons_are_not_balanced(self) -> None:
        exposed = [outcome(eligible=False, status="CENSORED_TIMEOUT")]
        control = [outcome(eligible=False, status="CENSORED_BACKEND_FAILURE")]
        comparison = self.analyzer.comparison_record("test", exposed, control)
        self.assertFalse(comparison["censoring_gate"]["balanced"])
        self.assertFalse(comparison["prospective_conclusion_authorized"])

    def test_exact_arithmetic_group_audit_replays_and_authorizes_only_v1(self) -> None:
        self.assertEqual(self.group_builder.build(), self.group_artifact)
        definition = self.group_artifact["definition"]
        self.assertEqual(
            definition["development_population"]["repeated_exact_j_group_count"], 0
        )
        self.assertEqual(
            definition["development_population"]["cross_split_exact_j_group_count"],
            0,
        )
        self.assertEqual(
            definition["prospective_holdout"][
                "twist_class_overlap_with_labelled_selection"
            ],
            0,
        )
        gate = definition["gate"]
        self.assertEqual(gate["status"], "PASS_EXACT_ISOMORPHISM_TWIST_GROUPING")
        self.assertTrue(gate["learned_score_reuse_authorized"])


if __name__ == "__main__":
    unittest.main()

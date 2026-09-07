from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-promotion-gate-v1.json"
)
SCRIPT = (
    ROOT
    / "elkies-k3/scripts/build_r17_prospective_crt_half_lattice_promotion_gate.py"
)


def canonical_hash(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


class R17ProspectiveCRTHalfLatticePromotionGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gate = json.loads(ARTIFACT.read_text())

    def test_artifact_replays_exactly(self) -> None:
        spec = importlib.util.spec_from_file_location("promotion_gate", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        self.assertEqual(module.build(), self.gate)
        body = {
            key: value
            for key, value in self.gate.items()
            if key not in {"gate_definition_sha256", "inputs", "generation"}
        }
        self.assertEqual(canonical_hash(body), self.gate["gate_definition_sha256"])

    def test_binary_endpoint_cannot_promote(self) -> None:
        self.assertEqual(
            self.gate["status"], "FAIL_CLOSED_NO_RANK32_CANDIDATE_PROMOTION"
        )
        self.assertFalse(self.gate["rank_32_candidate_promotion_authorized"])
        self.assertEqual(
            self.gate["rank_32_candidate_promotion_authorized"],
            all(self.gate["current_gate_checks"].values()),
        )
        self.assertTrue(self.gate["detector_score"]["binary_any_escape_forbidden"])
        self.assertIn(
            "at least one",
            self.gate["assessed_protocol"]["primary_event"],
        )

    def test_current_controls_do_not_validate_magnitude(self) -> None:
        calibration = self.gate["available_calibration"]
        self.assertEqual(
            calibration["sealed_holdout_distinct_displayed_jump_strata"], [12]
        )
        self.assertEqual(calibration["sealed_holdout_score_range_at_jump_12"], [3, 12])
        self.assertTrue(calibration["development_strict_order_reversals"])
        checks = self.gate["current_gate_checks"]
        self.assertFalse(
            checks[
                "score_was_validated_on_an_independent_panel_with_multiple_jump_strata"
            ]
        )
        self.assertFalse(checks["score_has_predeclared_upper_tail_acceptance_rule"])

    def test_rank32_and_arithmetic_gates_are_explicit(self) -> None:
        target = self.gate["rank_32_target"]
        self.assertEqual((target["generic_rank"], target["target_rank"]), (17, 32))
        self.assertEqual(target["required_quotient_jump"], 15)
        requirements = self.gate["promotion_requirements"]
        self.assertIn("upper-tail", requirements["tail_validation"])
        self.assertIn("residual 2-Selmer", requirements["existing_arithmetic_gate"])

    def test_prospective_effects_require_balanced_censoring(self) -> None:
        policy = self.gate["censoring_policy"]
        self.assertTrue(
            policy["scheduled_rows_may_not_be_counted_as_completed_non_events"]
        )
        self.assertIn("complete Stage-A", policy["primary_effect_denominator"])
        self.assertIn("exactly the same", policy["balance_rule"])
        self.assertIn("null", policy["failure_action"])
        self.assertTrue(
            self.gate["current_gate_checks"][
                "analyzer_requires_balanced_censoring_for_any_prospective_contrast"
            ]
        )

    def test_chart_scores_have_search_order_semantics_only(self) -> None:
        policy = self.gate["chart_order_policy"]
        self.assertEqual(
            policy["quartic_role"],
            "pointed_birational_search_chart_of_the_same_elliptic_curve",
        )
        self.assertFalse(policy["quartic_is_nontrivial_2_covering"])
        self.assertFalse(policy["quartic_represents_a_selmer_class"])
        self.assertTrue(
            all(value is False for value in policy["miss_inferences"].values())
        )
        legacy = self.gate["legacy_field_interpretation"]
        self.assertIn("not arithmetic depths", legacy["depth_fields"])
        self.assertIn("not a basis-invariant", legacy["old_deep_43"])
        self.assertIn("not a Selmer", legacy["quotient_hamming_weight"])

    def test_changed_lattice_requires_recomputed_order_and_revalidation(self) -> None:
        invalidation = self.gate["chart_order_policy"]["ordering_invalidation"]
        self.assertIn("discard the cached ordering", invalidation["required_action"])
        self.assertIn("none", invalidation["empirical_effectiveness_transfer"])
        requirement = self.gate["promotion_requirements"][
            "state_bound_chart_ordering"
        ]
        self.assertIn("After every", requirement)
        self.assertIn("lattice-state fingerprint", requirement)
        audit = self.gate["executable_chart_order_policy_audit"]
        self.assertTrue(audit["unchanged_exact_state_accepted"])
        self.assertTrue(audit["all_declared_test_changes_rejected"])
        self.assertEqual(
            audit["rejected_changes"],
            [
                "lattice_enlargement",
                "basis_change",
                "height_gram_change",
                "quotient_basis_change",
                "chart_order_change",
            ],
        )


if __name__ == "__main__":
    unittest.main()

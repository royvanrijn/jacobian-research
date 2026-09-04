from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_DIR = ROOT / "artifacts/generated-results"
MANIFEST = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
PROTOCOL = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-half-lattice-protocol-v3.json"
RUNNER = ROOT / "elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage"
ANALYZER = ROOT / "elkies-k3/scripts/analyze_r17_prospective_crt_half_lattice_experiment.py"

EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_hash(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


class R17ProspectiveCRTHalfLatticeProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST.read_text())
        cls.protocol = json.loads(PROTOCOL.read_text())

    def test_original_cohort_is_reused_unopened(self) -> None:
        self.assertEqual(
            self.manifest["commitment"]["candidate_list_sha256"],
            EXPECTED_CANDIDATE_HASH,
        )
        self.assertEqual(len(self.manifest["rows"]), 2_560)
        self.assertTrue(
            all(row["outcome_status"] == "NOT_OPENED" for row in self.manifest["rows"])
        )
        self.assertEqual(
            Counter(row["cohort"] for row in self.manifest["rows"]),
            {
                "A_356_full": 256,
                "B_385_full": 256,
                "C_matched_ordinary": 512,
                "D_two_only": 512,
                "E_odd_only": 512,
                "F_random_equal_codimension": 512,
            },
        )

    def test_protocol_hash_and_executables_are_pinned(self) -> None:
        protocol_body = {
            key: value
            for key, value in self.protocol.items()
            if key not in {"protocol_definition_sha256", "inputs", "generation"}
        }
        self.assertEqual(
            canonical_hash(protocol_body), self.protocol["protocol_definition_sha256"]
        )
        inputs = self.protocol["inputs"]
        self.assertEqual(
            inputs["elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage"],
            digest(RUNNER),
        )
        self.assertEqual(
            inputs["elkies-k3/scripts/analyze_r17_prospective_crt_half_lattice_experiment.py"],
            digest(ANALYZER),
        )

    def test_stage_a_and_stage_b_rules_are_fixed(self) -> None:
        self.assertEqual(self.protocol["native_generic_lattice"]["deepest_minimum_norm"], 12)
        self.assertEqual(self.protocol["native_generic_lattice"]["deepest_class_count"], 43)
        self.assertEqual(
            self.protocol["native_generic_lattice"]["exact_minimum_norm_histogram"],
            {"0": 1, "4": 1311, "6": 26672, "8": 63925, "10": 39120, "12": 43},
        )
        self.assertEqual(self.protocol["stage_a"]["covers_per_fibre"], 43)
        self.assertEqual(
            self.protocol["stage_b"]["gate"],
            "at least one Stage-A point has a full exact finite-reduction independence certificate beyond specialized MW17",
        )
        self.assertTrue(self.protocol["stage_b"]["no_other_escalation_rule"])

    def test_search_and_acceptance_limits_are_fail_closed(self) -> None:
        pipeline = self.protocol["cover_pipeline"]
        self.assertEqual(pipeline["height_bound_each_cover"], 100_000)
        self.assertEqual(
            pipeline["wall_timeout_seconds_each_cover_including_minimize_reduce_search"],
            15,
        )
        self.assertEqual(pipeline["retries"], 0)
        acceptance = self.protocol["point_acceptance"]
        self.assertTrue(acceptance["exact_original_short_curve_equation_required"])
        self.assertTrue(
            acceptance["combined_mod2_rank_must_equal_17_plus_all_counted_directions"]
        )
        self.assertTrue(acceptance["uncertified_returned_points_do_not_count_as_escapes"])

    def test_external_plus12_controls_pass_before_freeze(self) -> None:
        controls = self.protocol["positive_control_acceptance"]
        self.assertEqual(
            (
                controls["curve356-rank29"]["stage_a_exact_quotient_rank"],
                controls["curve356-rank29"]["stage_b_union_exact_quotient_rank"],
            ),
            (12, 12),
        )
        self.assertEqual(
            (
                controls["curve385-rank29"]["stage_a_exact_quotient_rank"],
                controls["curve385-rank29"]["stage_b_union_exact_quotient_rank"],
            ),
            (3, 4),
        )

    def test_primary_estimand_and_conditional_stage_b_are_predeclared(self) -> None:
        primary = self.protocol["primary_outcome"]
        self.assertIn("scheduled fibre", primary["primary_estimand"])
        self.assertTrue(
            primary["stage_b_is_conditional_recovery_depth_not_an_unconditional_event_comparison"]
        )


if __name__ == "__main__":
    unittest.main()

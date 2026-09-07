from __future__ import annotations

from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_DIR = ROOT / "artifacts/generated-results"
PHASE1 = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-local-stability-v1.json"
MANIFEST = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
FEATURES = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz"
PROTOCOL = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-search-protocol-v2.json"
POINTS = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-point-search-ledger-v2.json"
SENSITIVITY = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-search-sensitivity-v1.json"
ANALYSIS = ARTIFACT_DIR / "elkies-k3-r17-prospective-crt-statistical-analysis-v1.json"

EXPECTED_CYLINDER_HASH = "500dc6931c5aeaf3d6d9982bb994286d7aee36e7c87b9e414e8b7e0ef8aef15c"
EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
EXPECTED_PROTOCOL_HASH = "63d6b9e83f52bc7208b9057298e05941dfcedc85d53f5681186c953498947d4b"


def canonical_hash(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


class R17ProspectiveCRTExperimentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phase1 = json.loads(PHASE1.read_text())
        cls.manifest = json.loads(MANIFEST.read_text())
        with gzip.open(FEATURES, "rt") as source:
            cls.features = json.load(source)
        cls.protocol = json.loads(PROTOCOL.read_text())
        cls.points = json.loads(POINTS.read_text())
        cls.sensitivity = json.loads(SENSITIVITY.read_text())
        cls.analysis = json.loads(ANALYSIS.read_text())

    def test_refined_cylinders_are_frozen_but_empirical(self) -> None:
        self.assertEqual(
            self.phase1["frozen_cylinder_definition_sha256"], EXPECTED_CYLINDER_HASH
        )
        self.assertEqual(
            canonical_hash(self.phase1["frozen_cylinder_definition"]), EXPECTED_CYLINDER_HASH
        )
        cylinders = {
            row["anchor_curve_id"]: row
            for row in self.phase1["frozen_cylinder_definition"]["cylinders"]
        }
        self.assertEqual(
            [row["exponent"] for row in cylinders[356]["prime_power_conditions"]],
            [15, 5, 3, 3, 3],
        )
        self.assertEqual(
            [row["exponent"] for row in cylinders[385]["prime_power_conditions"]],
            [18, 4, 3, 3, 3],
        )
        self.assertTrue(all(row["empirical_only"] for row in cylinders.values()))
        self.assertFalse(
            self.phase1["frozen_cylinder_definition"]["selection_uses_point_search_outcomes"]
        )

    def test_candidate_commitment_replays_and_remains_unopened(self) -> None:
        fields = (
            "sample_id",
            "match_set_id",
            "anchor_curve_id",
            "cohort",
            "parameter",
            "projective_pair",
            "cylinder_residue",
            "cylinder_modulus",
        )
        payload = [{key: row[key] for key in fields} for row in self.manifest["rows"]]
        self.assertEqual(canonical_hash(payload), EXPECTED_CANDIDATE_HASH)
        self.assertEqual(
            self.manifest["commitment"]["candidate_list_sha256"], EXPECTED_CANDIDATE_HASH
        )
        self.assertEqual(len(self.manifest["rows"]), 2560)
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
        self.assertTrue(
            all(row["outcome_status"] == "NOT_OPENED" for row in self.manifest["rows"])
        )

    def test_presearch_panel_covers_every_frozen_row_and_keeps_selmer_open(self) -> None:
        self.assertEqual(self.features["candidate_list_sha256"], EXPECTED_CANDIDATE_HASH)
        self.assertEqual(len(self.features["rows"]), 2560)
        self.assertEqual(self.features["summary"]["bounded_search_authorizations"], 2560)
        self.assertEqual(self.features["summary"]["complete_two_selmer_groups_computed"], 0)
        self.assertEqual(
            self.features["summary"]["finite_proved_residual_upper_bounds_computed"], 0
        )
        survival = self.features["summary"]["anchor_fingerprint_matched_place_count_by_cohort"]
        self.assertEqual(survival["A_356_full"], {"5": 256})
        self.assertEqual(survival["B_385_full"], {"5": 256})
        self.assertEqual(survival["E_odd_only"], {"4": 512})
        for row in self.features["rows"]:
            self.assertIsNone(row["selmer_measurement_status"]["proved_residual_upper_bound"])

    def test_backend_amendment_preceded_any_point_outcome(self) -> None:
        self.assertEqual(self.protocol["protocol_definition_sha256"], EXPECTED_PROTOCOL_HASH)
        self.assertTrue(
            self.protocol["frozen_before_any_point_search_call_completed_or_returned_points"]
        )
        self.assertFalse(self.protocol["selection_or_rebalancing_changed"])
        canary = self.protocol["infeasible_v1_canary"]
        self.assertFalse(canary["point_search_call_reached"])
        self.assertFalse(canary["points_returned_or_inspected"])

    def test_complete_point_ledger_preserves_zero_events_and_all_failures(self) -> None:
        self.assertEqual(self.points["candidate_list_sha256"], EXPECTED_CANDIDATE_HASH)
        self.assertEqual(self.points["search_protocol_sha256"], EXPECTED_PROTOCOL_HASH)
        self.assertEqual(len(self.points["records"]), 2560)
        self.assertEqual(
            Counter(row["status"] for row in self.points["records"]),
            {"BOUNDED_PROTOCOL_NO_ESCAPE_FOUND": 2560},
        )
        self.assertEqual(self.points["summary"]["certified_escape_rows"], 0)
        self.assertEqual(self.points["summary"]["certified_extra_directions"], 0)

    def test_zero_event_analysis_does_not_claim_rank_or_selmer_bounds(self) -> None:
        primary = self.analysis["analysis"]["primary_comparison"]
        self.assertEqual(primary["exposed"]["certified_escape_rows"], 0)
        self.assertEqual(primary["control"]["certified_escape_rows"], 0)
        self.assertEqual(primary["risk_difference"], 0.0)
        self.assertIsNone(primary["risk_ratio"])
        self.assertIsNone(primary["odds_ratio"])
        self.assertEqual(primary["fisher_exact_two_sided_p"], 1.0)
        self.assertEqual(
            self.analysis["analysis"]["predictor_fit"]["status"],
            "NOT_FIT_ZERO_OUTCOME_VARIANCE",
        )
        boundary = " ".join(self.analysis["claim_boundary"])
        self.assertIn("not rank 17", boundary)
        self.assertIn("no Selmer upper bound", boundary)

    def test_post_experiment_positive_controls_show_detector_limit(self) -> None:
        self.assertTrue(self.sensitivity["post_experiment_diagnostic_only"])
        self.assertFalse(
            self.sensitivity["changes_frozen_candidates_protocol_or_contrasts"]
        )
        self.assertEqual(
            self.sensitivity["status"],
            "FAILED_TO_REDETECT_BOTH_KNOWN_PLUS12_POSITIVE_CONTROLS",
        )
        self.assertEqual(
            [(row["curve_id"], row["certified_escape_count"]) for row in self.sensitivity["controls"]],
            [(356, 0), (385, 0)],
        )
        self.assertEqual(
            self.analysis["analysis"]["chain_assessment"]["prospective_enrichment_claim"],
            "NO_EVIDENCE_DETECTOR_LIMITED_AT_THE_FROZEN_SEARCH_BOUND",
        )


if __name__ == "__main__":
    unittest.main()

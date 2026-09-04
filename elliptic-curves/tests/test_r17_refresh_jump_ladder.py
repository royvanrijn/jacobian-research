from __future__ import annotations

from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
PROTOCOL1 = ART / "r17_refresh_jump_ladder_protocol_v1.json"
BLIND1 = ART / "r17_refresh_jump_ladder_blind_v1.json"
PROTOCOL2 = ART / "r17_refresh_jump_ladder_protocol_v2.json"
BLIND2 = ART / "r17_refresh_jump_ladder_blind_v2.json"
VERIFIED2 = ART / "r17_refresh_jump_ladder_verification_v2.json"
ANALYSIS2 = ART / "r17_refresh_jump_ladder_analysis_v2.json"

EXPECTED_HASHES = {
    DATA: "ee4333d0c97737a2beed4487fa215b18ac9e18b5711efba42dcc7761446ae037",
    PROTOCOL1: "c348d9134290c354fe962dd822b5b0f2355e8a1c747168f95e880447129a287a",
    BLIND1: "15eae6f0c65132841befd161ec3bb8657c0a086dd224d47eafd1e5568f305e53",
    PROTOCOL2: "d080c990d32b24e76695db393db276d29b620462478f1840efc0bb122537d737",
    BLIND2: "0699b53c2bc7d77673231bc0d377dc725880efd26a18d1eb2af613d28578c165",
    VERIFIED2: "06d3e0a097a40bad50d6a0698fd7dbb4a06021149b570f43e583152eddf3eaa4",
    ANALYSIS2: "67889bc834aa6b9e02d31e982b521e56bc7ad5960911c7dc2a4188e3a8a8b034",
}

EXPECTED_RESPONSE = {
    478: 6,
    498: 6,
    531: 11,
    532: 3,
    534: 11,
    535: 10,
    536: 11,
    537: 10,
    538: 5,
    539: 6,
    540: 8,
    541: 8,
    543: 12,
    544: 0,
    545: 11,
    546: 8,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class R17RefreshJumpLadderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(DATA.read_text())
        cls.protocol1 = json.loads(PROTOCOL1.read_text())
        cls.blind1 = json.loads(BLIND1.read_text())
        cls.protocol2 = json.loads(PROTOCOL2.read_text())
        cls.blind2 = json.loads(BLIND2.read_text())
        cls.verified2 = json.loads(VERIFIED2.read_text())
        cls.analysis2 = json.loads(ANALYSIS2.read_text())

    def test_frozen_bytes(self):
        for path, expected in EXPECTED_HASHES.items():
            self.assertEqual(digest(path), expected, str(path))

    def test_redacted_input_boundary(self):
        self.assertEqual(self.data["status"], "FROZEN_MW17_ONLY_NO_PUBLIC_COMPLEMENT")
        self.assertEqual(self.data["case_count"], 16)
        self.assertFalse(self.data["redaction"]["contains_displayed_complement_coordinates"])
        self.assertFalse(self.data["redaction"]["contains_displayed_rank_or_jump"])
        self.assertEqual(self.data["pre_search_exclusion"]["curve_id"], 499)
        permitted = {
            "curve_id",
            "representative_class",
            "native_chart",
            "short_model",
            "generic_points",
            "generic_height_gram",
            "generic_height_gram_determinant",
        }
        self.assertTrue(all(set(row) == permitted for row in self.data["cases"]))

    def test_cross_class_amendment_is_disclosed(self):
        self.assertEqual(self.protocol1["status"], "FROZEN_BEFORE_BLIND_RECOVERY")
        self.assertEqual(len(self.blind1["results"]), 1)
        self.assertEqual(self.blind1["results"][0]["curve_id"], 478)
        amendment = self.protocol2["cross_class_amendment"]
        self.assertTrue(amendment["curve478_selected_set_and_order_unchanged"])
        self.assertFalse(amendment["remaining_fifteen_recovery_outcomes_known_at_amendment"])
        self.assertTrue(amendment["confirmatory_endpoints_and_acceptance_rules_unchanged"])

    def test_exact_blind_response_and_equal_caps(self):
        self.assertEqual(
            self.blind2["status"],
            "PASS_COMPLETE_BLIND_RECOVERY_BEFORE_PUBLIC_COMPLEMENT",
        )
        observed = {
            int(row["curve_id"]): int(
                row["exact_quotient_rank_recovered_before_public_complement"]
            )
            for row in self.blind2["results"]
        }
        self.assertEqual(observed, EXPECTED_RESPONSE)
        for row in self.blind2["results"]:
            initial = int(row["initial"]["exact_quotient_rank_recovered"])
            adaptive = int(
                (row.get("adaptive") or {}).get(
                    "exact_incremental_quotient_rank_recovered", 0
                )
            )
            self.assertEqual(
                initial + adaptive,
                row["exact_quotient_rank_recovered_before_public_complement"],
            )
            if row["curve_id"] == 544:
                self.assertEqual(row["attempted_chart_count"], 43)
                self.assertEqual(initial, 0)
            else:
                self.assertEqual(row["attempted_chart_count"], 344)
                self.assertEqual(row["timeout_chart_count"], 0)
                self.assertEqual(row["pari_failure_chart_count"], 0)

    def test_public_complement_opened_after_blind_hash(self):
        boundary = self.verified2["phase_boundary"]
        self.assertTrue(
            boundary["public_complement_opened_only_after_blind_status_and_hash_were_fixed"]
        )
        self.assertEqual(
            boundary["blind_artifact_sha256_before_public_complement_import"],
            EXPECTED_HASHES[BLIND2],
        )
        failures = {
            int(row["curve_id"])
            for row in self.verified2["results"]
            if not row["all_final_blind_basis_points_in_opened_public_subgroup"]
        }
        self.assertEqual(failures, {478, 539})

    def test_predeclared_confirmatory_endpoints_pass(self):
        self.assertEqual(self.analysis2["status"], "PASS_USABLE_EXTREME_JUMP_DETECTOR")
        ordinal = self.analysis2["confirmatory"]["ordinal_association"]
        self.assertTrue(ordinal["passed"])
        self.assertEqual(ordinal["concordance_minus_discordance"], 78)
        self.assertEqual((ordinal["exact_p_numerator"], ordinal["exact_p_denominator"]), (60852, 2421619200))
        self.assertAlmostEqual(ordinal["kendall_tau_b"], 0.7503122325921043)
        tail = self.analysis2["confirmatory"]["upper_tail_enrichment"]
        self.assertTrue(tail["passed"])
        self.assertEqual(
            tail["table"],
            {
                "true_tail_detector_positive": 7,
                "true_tail_detector_negative": 1,
                "non_tail_detector_positive": 0,
                "non_tail_detector_negative": 8,
            },
        )
        self.assertEqual((tail["exact_p_numerator"], tail["exact_p_denominator"]), (1, 1430))
        self.assertTrue(self.analysis2["confirmatory"]["joint_decision"]["passed"])

    def test_analysis_replays(self):
        source = ROOT / "elliptic-curves/cas/analyze_r17_refresh_jump_ladder_v2.py"
        module = SourceFileLoader("r17_jump_ladder_analysis_v2_test", str(source)).load_module()
        self.assertEqual(module.build(), self.analysis2)

    def test_post_freeze_fibration_and_j_class_sensitivity(self):
        split = self.analysis2["descriptive"][
            "post_freeze_fibration_and_j_class_sensitivity"
        ]
        self.assertEqual(split["role"], "descriptive_post_freeze_not_confirmatory")
        q11 = split["pooled_q_at_least_11"]
        self.assertEqual(
            q11["table"],
            {
                "true_tail_detector_positive": 6,
                "true_tail_detector_negative": 1,
                "non_tail_detector_positive": 1,
                "non_tail_detector_negative": 8,
            },
        )
        self.assertEqual((q11["exact_p_numerator"], q11["exact_p_denominator"]), (4, 715))

        frames = {row["stratum"]: row for row in split["by_fibration"]}
        self.assertEqual(frames["published-R17"]["case_count"], 13)
        self.assertEqual(frames["alternate-Q80"]["case_count"], 3)
        self.assertFalse(
            frames["alternate-Q80"]["high_S_enrichment"]["q_at_least_10"][
                "estimable"
            ]
        )

        j_classes = {row["stratum"]: row for row in split["by_j_map_class"]}
        class_08234 = j_classes["norm12-orbit-08234"]
        self.assertEqual(class_08234["case_count"], 7)
        self.assertEqual(
            (
                class_08234["ordinal_association"]["exact_p_numerator"],
                class_08234["ordinal_association"]["exact_p_denominator"],
            ),
            (2, 210),
        )
        self.assertEqual(
            class_08234["high_S_enrichment"]["q_at_least_10"]["table"],
            {
                "true_tail_detector_positive": 5,
                "true_tail_detector_negative": 0,
                "non_tail_detector_positive": 0,
                "non_tail_detector_negative": 2,
            },
        )
        conditional = split["conditional_on_j_map_class"]
        self.assertEqual(
            (
                conditional["ordinal_association"]["exact_p_numerator"],
                conditional["ordinal_association"]["exact_p_denominator"],
            ),
            (2, 2520),
        )
        self.assertEqual(
            (
                conditional["high_S_enrichment"]["q_at_least_10"][
                    "exact_p_numerator"
                ],
                conditional["high_S_enrichment"]["q_at_least_10"][
                    "exact_p_denominator"
                ],
            ),
            (3, 63),
        )
        self.assertEqual(
            (
                conditional["high_S_enrichment"]["q_at_least_11"][
                    "exact_p_numerator"
                ],
                conditional["high_S_enrichment"]["q_at_least_11"][
                    "exact_p_denominator"
                ],
            ),
            (9, 63),
        )


if __name__ == "__main__":
    unittest.main()

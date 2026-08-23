from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys
import unittest

import pypdf


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from audit_rank21_conductor_public_data import (  # noqa: E402
    EXPECTED_EXACT_CONDUCTORS,
    SOURCE_SPECS,
    STRICT_LOG_CONDUCTOR_TARGET,
    verify_payload,
)


REPOSITORY = Path(__file__).resolve().parents[2]
ARTIFACT = (
    REPOSITORY
    / "artifacts/generated-results/"
    "elliptic_rank21_conductor_public_data_audit.json"
)
SCRIPT = CAS_DIRECTORY / "audit_rank21_conductor_public_data.py"


class Rank21ConductorPublicDataAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_pinned_source_hashes_and_software_metadata(self) -> None:
        report = self.report
        inventories = report["source_inventories"]
        self.assertEqual(set(inventories), set(SOURCE_SPECS))
        for name, spec in SOURCE_SPECS.items():
            self.assertEqual(inventories[name]["sha256"], spec.sha256)
            self.assertEqual(inventories[name]["url"], spec.url)
            self.assertTrue(inventories[name]["matches_pinned_revision"])
        self.assertEqual(
            report["reproduction"]["script_sha256"],
            sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertTrue(
            report["reproduction"]["all_source_hashes_verified_before_write"]
        )
        # This is a pinned 2026-08-14 source snapshot, not a rolling artifact.
        # Requiring the recorded interpreter to equal the current environment
        # makes an unchanged historical audit fail after every Python update.
        self.assertEqual(report["software"]["python"], "3.14.6")
        self.assertEqual(report["software"]["pypdf"], pypdf.__version__)
        self.assertRegex(report["software"]["pari_gp"], r"^\[\d+, \d+, \d+\]$")

    def test_all_printed_rank21_through_rank29_points_are_exactly_on_curve(self) -> None:
        pages = self.report["record_page_replays"]
        self.assertEqual(set(pages), {str(rank) for rank in range(21, 30)})
        for rank in range(21, 30):
            page = pages[str(rank)]
            self.assertEqual(page["rank_lower_bound_printed_by_source"], rank)
            self.assertEqual(page["printed_point_count"], rank)
            self.assertEqual(page["exact_membership_checks_passed"], rank)
            self.assertEqual(page["printed_point_indices"], list(range(1, rank + 1)))
            self.assertFalse(page["independence_reproved_from_this_page_alone"])
            self.assertEqual(len(page["weierstrass_coefficients"]), 5)

    def test_exact_local_conductor_replays_and_factorizations(self) -> None:
        by_rank = {
            item["rank_lower_bound_printed_by_source"]: item
            for item in self.report["candidate_inventory"]
        }
        self.assertEqual(set(by_rank), set(range(21, 30)))
        for rank, expected_conductor in EXPECTED_EXACT_CONDUCTORS.items():
            conductor = by_rank[rank]["conductor_evidence"]
            self.assertEqual(conductor["status"], "exact_local_pari_gp_replay")
            self.assertEqual(conductor["conductor"], str(expected_conductor))
            self.assertTrue(conductor["factor_product_verified_before_gp"])
            self.assertTrue(conductor["all_factor_bases_proved_prime_by_gp"])
            product = 1
            for factor in conductor["complete_discriminant_factorization"]:
                product *= int(factor["prime"]) ** factor["exponent"]
            self.assertEqual(product, abs(int(conductor["minimal_discriminant"])))
            self.assertGreaterEqual(
                float(conductor["log_conductor"]),
                float(STRICT_LOG_CONDUCTOR_TARGET),
            )
            self.assertFalse(conductor["below_strict_log_conductor_target"])

        for rank in (26, 27, 28):
            conductor = by_rank[rank]["conductor_evidence"]
            self.assertEqual(conductor["status"], "dujella_machine_table_log_only")
            self.assertGreater(float(conductor["margin_above_strict_target"]), 65.0)

    def test_current_authoritative_database_scopes_have_no_hidden_hit(self) -> None:
        results = self.report["source_results"]
        torsion = results["dujella_torsion_frontier"]
        self.assertEqual(torsion["trivial_torsion_frontier"], 29)
        self.assertEqual(torsion["largest_nontrivial_torsion_frontier"], 20)
        self.assertFalse(torsion["any_nontrivial_torsion_frontier_at_least_21"])

        history = results["dujella_rank_history"]
        self.assertEqual(history["highest_linked_record_rank"], 29)
        self.assertFalse(history["rank30_page_linked"])

        machine = results["dujella_machine_tables"]
        self.assertEqual(machine["coefficient_vector_count"], 1179)
        self.assertEqual(machine["conductor_log_count"], 1180)
        self.assertTrue(machine["known_z6_row_count_mismatch_in_pinned_download"])
        self.assertEqual(
            [item["rank"] for item in machine["z1_current_and_predecessor_records"]],
            [29, 28, 27, 26, 25, 24],
        )

        lmfdb = results["lmfdb"]
        self.assertEqual(lmfdb["curve_count"], 3_824_372)
        self.assertEqual(lmfdb["isogeny_class_count"], 2_917_287)
        self.assertEqual(lmfdb["largest_conductor_in_database"], 299_996_953)
        self.assertEqual(lmfdb["largest_rank_displayed"], 5)
        self.assertFalse(lmfdb["rank_at_least_21_entry_present"])

        elkies = results["elkies_author_small_conductor_table"]
        self.assertEqual(elkies["displayed_rank_range"], [0, 11])
        self.assertFalse(elkies["rank21_in_scope"])
        watkins = results["primary_pdfs"]["elkies_watkins"]
        self.assertEqual(watkins["abstract_search_rank_range"], [5, 11])
        self.assertFalse(watkins["rank21_in_scope"])

    def test_fermigier_primary_paper_has_no_adjacent_exact_specializations(self) -> None:
        fermigier = self.report["source_results"]["primary_pdfs"]["fermigier"]
        self.assertEqual(fermigier["page_count"], 5)
        self.assertEqual(
            fermigier["rational_t_specializations_printed"], ["19754/39"]
        )
        self.assertTrue(fermigier["only_printed_specialization_is_record_fiber"])
        self.assertFalse(
            fermigier["adjacent_E19_E20_E21_exact_models_or_parameters_printed"]
        )
        self.assertTrue(fermigier["adjacent_rows_are_score_only"])
        self.assertEqual(
            set(fermigier["comparison_score_rows"]), {"E19", "E20", "E21", "E22"}
        )
        self.assertTrue(
            all(len(row) == 8 for row in fermigier["comparison_score_rows"].values())
        )
        adjacent = self.report["fermigier_adjacent_specialization_result"]
        self.assertEqual(adjacent["exact_new_adjacent_specializations_recovered"], [])
        self.assertIn("does not publish", adjacent["archive_blocker"])

    def test_local_duplicate_filter_and_final_assessment(self) -> None:
        inventory = self.report["candidate_inventory"]
        local_ranks = {
            item["rank_lower_bound_printed_by_source"]
            for item in inventory
            if item["excluded_as_already_local"]
        }
        self.assertEqual(local_ranks, {21, 22, 29})
        new_ranks = {
            item["rank_lower_bound_printed_by_source"]
            for item in inventory
            if item["new_to_local_artifacts_at_audit_time"]
        }
        self.assertEqual(new_ranks, {23, 24, 25, 26, 27, 28})
        self.assertTrue(
            all(not item["plausible_target_candidate_after_conductor_gate"] for item in inventory)
        )
        for item in inventory:
            for match in item["local_artifact_isomorphism_matches"]:
                path = REPOSITORY / match["path"]
                self.assertEqual(match["sha256"], sha256(path.read_bytes()).hexdigest())

        assessment = self.report["assessment"]
        self.assertFalse(assessment["new_public_model_meeting_both_rank_and_conductor_gates"])
        self.assertFalse(assessment["finite_reduction_certificate_triggered_for_new_curve"])
        self.assertFalse(assessment["breakthrough_result_found"])
        self.assertEqual(
            assessment["closest_public_certified_near_miss"]["rank_lower_bound"], 22
        )
        self.assertEqual(
            assessment["closest_public_certified_near_miss"]["excess_over_strict_target"],
            "0.004910950637428796108330351524",
        )
        self.assertIn("bounded negative", assessment["epistemic_status"])

    def test_hash_guard_rejects_unpinned_payload(self) -> None:
        with self.assertRaisesRegex(ValueError, "changed"):
            verify_payload("dujella_rank21", b"not the pinned rank-21 page")


if __name__ == "__main__":
    unittest.main()

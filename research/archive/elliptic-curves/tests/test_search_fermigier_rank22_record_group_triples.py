#!/usr/bin/env python3
"""Focused checks for the exact weight-three record-group slice tranche."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_fermigier_rank22_record_group_triples.py"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_rank22_record_group_triples.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "f28040e5643ccc9f3336f60de23e2d67802052a5f234bcc2e2b8567c3c2a2a8c"
)
EXPECTED_ARTIFACT_SHA256 = (
    "2803b1fa276c80eccceac5ce83215f8678d9fb771abdccb8e5043a9962b1ed36"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierRank22RecordGroupTriplesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.slice_rows = [
            row
            for direction in cls.data["direction_searches"]
            for row in direction["slice_searches"]
        ]

    def test_pinned_files_and_exact_lower_level_frontier(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        source = self.data["source"]
        self.assertEqual(
            source["published_accidental_preimage_sha256"],
            "6224da9ce4db3150a197a2cf1d9bc6c1a7d0cc6f01245b3f834945f76775ab15",
        )
        self.assertEqual(
            source["transport_source_sha256"],
            "8650d548b7c101514d3403a46aa8ea49cef071b231d86df2e3950be84a868011",
        )
        self.assertEqual(
            source["certified_published_rank22_basis_point_sha256"],
            "4b2b89f3f432be8599b6ab5109c1221af52d1aa96c05eba1229983d929c9e727",
        )
        self.assertEqual(
            source["previous_l1_le_2_artifact_sha256"],
            "4928c44e27cada74b7a558dd97edfba20a554a508ac0e98ca051df7dea66a3c1",
        )
        self.assertEqual(
            source["previous_exact_slice_result_sha256"],
            "867ebd2e8c2b1a3608ea5577dbaa356de6c1b8e2574c45290fcbbbbf2ea4baf5",
        )
        prior_x = self.data["prior_record_quartic_abscissas"]
        self.assertEqual(prior_x["weight2_direction_count"], 462)
        self.assertEqual(prior_x["full_prior_x_count"], 494)
        self.assertEqual(
            prior_x["full_prior_x_sha256"],
            "62fcbbda7493f940b8001582929fedf27e822fc6d40f6ddf15adf8b4d820dc82",
        )

    def test_all_6160_exact_triples_are_unique_and_digest_pinned(self) -> None:
        population = self.data["full_triple_population"]
        self.assertTrue(population["global_sign_quotient"])
        self.assertEqual(population["coefficient_alphabet"], [-1, 0, 1])
        self.assertEqual(population["exact_l1_norm"], 3)
        self.assertEqual(population["full_vector_count"], 6_160)
        self.assertEqual(population["exceptional_inverse_count"], 0)
        self.assertEqual(population["prior_quartic_abscissa_exclusions"], 0)
        self.assertEqual(population["duplicate_new_abscissas"], 0)
        self.assertEqual(population["genuinely_new_unique_abscissa_count"], 6_160)
        self.assertEqual(
            population["full_vector_direction_sha256"],
            "76bf3cfc5de64af12865209492bf45616470a3e0f5d19e5527551989e1fb24d3",
        )
        self.assertEqual(
            population["genuinely_new_direction_sha256"],
            "4673b556e0a60943d86b23c7f293bfc3a9952f6acd25ef5e152297e84a106455",
        )
        self.assertEqual(population["minimum_new_abscissa_projective_height"], 69_550_814)

    def test_rank_blind_selection_is_exact_and_honestly_scoped(self) -> None:
        selection = self.data["rank_blind_selection"]
        self.assertFalse(selection["selection_uses_specialized_rank"])
        self.assertFalse(selection["selection_uses_conductor"])
        self.assertEqual(selection["lowest_height_keep"], 192)
        self.assertEqual(selection["highest_modular_keep"], 192)
        self.assertEqual(selection["per_first_support_index_keep"], 4)
        self.assertEqual(selection["per_relative_sign_pattern_keep"], 16)
        self.assertEqual(selection["selected_direction_count"], 399)
        self.assertEqual(selection["unselected_direction_count"], 5_761)
        self.assertEqual(
            selection["full_modular_score_sha256"],
            "13b003c942375155b95f8f3cf19d6a39518aee002db8c73552ed77af30602915",
        )
        self.assertEqual(
            selection["selected_direction_sha256"],
            "f22867994c80eecd7f33b48b3bf1d788cf128a9d86d77352187cbc6a9df9815e",
        )
        self.assertEqual(
            selection["selected_population_and_strata_sha256"],
            "5d5a1c4f0d2557ddd2ea43eab1c4bb09fca85e85f25b1f2de9c304191fd54134",
        )
        audit = selection["postrun_independent_selector_audit"]
        self.assertEqual(audit["selected_slope_usable_prime_count_minimum"], 6)
        self.assertEqual(audit["selected_slope_usable_prime_count_maximum"], 10)
        self.assertEqual(audit["selected_basis_participation_minimum"], 37)
        self.assertEqual(audit["selected_basis_participation_maximum"], 69)
        self.assertEqual(
            audit["selected_relative_sign_pattern_counts"],
            {"+1,+1": 58, "+1,-1": 155, "-1,+1": 118, "-1,-1": 68},
        )
        self.assertEqual(len(audit["audit_caveats"]), 4)
        provenance = self.data["execution_provenance"]
        self.assertEqual(
            provenance["bounded_execution_script_sha256"],
            "1147c06cecec8cccd7439357ccffc98acda0574c89031bb4d80236bb40107d4f",
        )
        self.assertEqual(
            provenance["stable_replay_script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        self.assertTrue(provenance["metadata_normalization_only"])
        self.assertEqual(provenance["bounded_searches_rerun_during_normalization"], 0)
        selected = self.data["selected_directions"]
        self.assertEqual(len(selected), 399)
        self.assertEqual(len({row["direction_id"] for row in selected}), 399)
        self.assertTrue(all(row["selection_strata"] for row in selected))
        for row in selected:
            vector = row["coefficient_vector"]
            self.assertEqual(sum(abs(value) for value in vector), 3)
            self.assertEqual(next(value for value in vector if value), 1)
        outcome = self.data["outcome"]
        self.assertEqual(outcome["full_genuinely_new_triple_direction_count"], 6_160)
        self.assertEqual(outcome["searched_selected_direction_count"], 399)
        self.assertEqual(outcome["unsearched_triple_direction_count"], 5_761)
        self.assertFalse(outcome["unsearched_triples_claimed_negative"])

    def test_all_798_selected_slices_completed_once(self) -> None:
        self.assertEqual(len(self.data["direction_searches"]), 399)
        self.assertEqual(len(self.slice_rows), 798)
        self.assertEqual(
            Counter(row["search"]["status"] for row in self.slice_rows),
            Counter({"completed": 798}),
        )
        self.assertTrue(all(not row["search"]["retried"] for row in self.slice_rows))
        self.assertTrue(
            all(row["record_T0_calibration_count"] == 1 for row in self.slice_rows)
        )
        self.assertEqual(
            Counter(row["search"]["signed_point_count"] for row in self.slice_rows),
            Counter({2: 797, 4: 1}),
        )
        outcome = self.data["outcome"]
        self.assertEqual(outcome["slice_calls_attempted"], 798)
        self.assertEqual(outcome["slice_calls_completed"], 798)
        self.assertEqual(outcome["slice_calls_timed_out_or_errored"], 0)
        self.assertEqual(outcome["record_T0_calibrated_slices"], 798)
        self.assertEqual(
            outcome["incidence_classification_counts"],
            {
                "genuinely-new-forced-fibre": 1,
                "record-source-calibration-excluded": 798,
            },
        )
        self.assertEqual(
            outcome["exact_selected_slice_result_sha256"],
            "dab194537af8ce6d252ea8d4634de343e35cc89285ebed763ad3922786e22f75",
        )

    def test_new_subtarget_fibre_is_exact_but_rank12(self) -> None:
        self.assertEqual(len(self.data["candidates"]), 1)
        candidate = self.data["candidates"][0]
        self.assertEqual(candidate["parameter_t"], "29771/78")
        self.assertEqual(candidate["signed_parameters"], ["-29771/78"])
        self.assertEqual(candidate["source_direction_ids"], ["p02_p06_m12"])
        self.assertEqual(candidate["distinct_forced_quartic_abscissa_count"], 1)
        conductor = candidate["conductor_probe"]
        self.assertEqual(conductor["status"], "completed")
        self.assertTrue(conductor["below_strict_log_conductor_target"])
        self.assertEqual(
            conductor["log_conductor"],
            "158.467530623289488846827352371228153283265542653031750991322",
        )
        self.assertEqual(conductor["root_number"], 1)
        forced = candidate["forced_points"][0]
        parameter = Fraction(candidate["parameter_t"])
        x_value = Fraction(forced["quartic_x"])
        z_value = Fraction(forced["quartic_z"])
        # Import locally to keep the hash-only tests lightweight.
        import sys

        sys.path.insert(0, str(ROOT / "elliptic-curves" / "cas"))
        from fermigier_mestre import FermigierMestreFamily

        self.assertEqual(z_value**2, FermigierMestreFamily.quartic_value(parameter, x_value))
        triage = candidate["rank_triage"]
        self.assertEqual(triage["generic_seed_stable_numerical_rank"], 12)
        self.assertEqual(triage["forced_augmented_stable_numerical_rank"], 12)
        self.assertEqual(triage["full_pool_stable_numerical_rank"], 12)
        self.assertEqual(triage["numerical_rank_gain_from_forced_directions"], 0)
        self.assertEqual(triage["numerical_rank_gain_after_H50000"], 0)
        self.assertEqual(len(triage["new_search_group_pullbacks"]), 2)
        self.assertNotIn("finite_reduction_attempt", triage)

    def test_clean_negative_target_handoff(self) -> None:
        outcome = self.data["outcome"]
        self.assertEqual(outcome["genuinely_new_forced_fibres"], 1)
        self.assertEqual(outcome["completed_conductors"], 1)
        self.assertEqual(outcome["subtarget_conductors"], 1)
        self.assertEqual(outcome["rank_triage_count"], 1)
        self.assertEqual(outcome["maximum_stable_numerical_rank"], 12)
        self.assertFalse(self.data["target"]["hit"])
        self.assertEqual(self.data["execution"]["owned_processes_remaining"], 0)
        self.assertTrue(self.data["parameters"]["no_retries"])
        prior = self.data["prior_parameter_decontamination"]
        self.assertEqual(prior["terminal_prior_parameter_count"], 1_239)
        self.assertEqual(
            prior["terminal_prior_parameter_sha256"],
            "9482e61650aa8bb1fd45c3765e5db92c1474090faee8d831e0d73cee4fc864c4",
        )


if __name__ == "__main__":
    unittest.main()

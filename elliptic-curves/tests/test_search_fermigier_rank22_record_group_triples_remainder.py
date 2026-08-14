#!/usr/bin/env python3
"""Focused checks for the exhaustive weight-three remainder tranche."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
SCRIPT = CAS / "search_fermigier_rank22_record_group_triples_remainder.py"
ARTIFACT = (
    GENERATED / "elliptic_fermigier_rank22_record_group_triples_remainder.json"
)
STREAM = (
    GENERATED
    / "elliptic_fermigier_rank22_record_group_triples_remainder_stream.jsonl"
)
SELECTED_ARTIFACT = (
    GENERATED / "elliptic_fermigier_rank22_record_group_triples.json"
)
AUXILIARY_ARTIFACT = (
    GENERATED / "elliptic_fermigier_rank22_auxiliary_orbits.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "c3cfe703bf2659e8ff5deea574aabdcf025a72d4e91276b8eb06837a88968f82"
)
EXPECTED_ARTIFACT_SHA256 = (
    "829bf44e50fb8d3583592190732ce4c9057c1bc54c9ababbea4e6cc002e6e028"
)
EXPECTED_STREAM_SHA256 = (
    "66e69019bb53b28310bc1a4fa0d40989fe2fd95e677b53bca55feaafa2f3b5de"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direction_digest(records: list[dict]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda item: item["direction_id"]):
        digest.update(
            (
                f"{record['direction_id']}|{record['quartic_x']}|"
                f"{record['quartic_z']}|"
                f"{','.join(map(str, record['coefficient_vector']))}\n"
            ).encode()
        )
    return digest.hexdigest()


class FermigierRank22RecordGroupTripleRemainderTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.rows = [json.loads(line) for line in STREAM.read_text().splitlines()]
        cls.slices = [
            row
            for direction in cls.rows
            for row in direction["slice_searches"]
        ]

    def test_pinned_files_and_honest_execution_provenance(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(sha256(STREAM), EXPECTED_STREAM_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        provenance = self.data["execution_provenance"]
        self.assertEqual(
            provenance["bounded_execution_script_sha256"],
            "bb11ad0b5460d1b10518cf2d8620da255447a3b3f616e7c184a27f73703f3c71",
        )
        self.assertEqual(
            provenance[
                "bounded_execution_artifact_sha256_before_metadata_normalization"
            ],
            "f8aad4a0ed562c71ed6f6f8a15b8d229dac2d2cb60601d501aeeffd49c3df891",
        )
        self.assertEqual(
            provenance["bounded_execution_auxiliary_orbit_artifact_sha256"],
            "1008336232ac65bb2bace6ff7008ffb18c1b491c6f0295e16fda14949d5b94d6",
        )
        self.assertEqual(
            provenance["stable_replay_script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        self.assertTrue(provenance["metadata_normalization_only"])
        self.assertEqual(
            provenance["bounded_slice_searches_rerun_during_normalization"], 0
        )

    def test_exact_5761_remainder_closes_the_full_6160_population(self) -> None:
        population = self.data["population"]
        self.assertEqual(population["full_exact_weight3_direction_count"], 6_160)
        self.assertEqual(population["previously_searched_direction_count"], 399)
        self.assertEqual(population["declared_remainder_direction_count"], 5_761)
        self.assertEqual(population["declared_remainder_slice_call_count"], 11_522)
        self.assertFalse(population["score_used_for_pruning"])
        self.assertEqual(
            population["remainder_direction_sha256"],
            "68a5922267e96f8da52f61235174c2423ce3b98e3c83d8785b00cf15d5336047",
        )
        self.assertEqual(len(self.rows), 5_761)
        self.assertEqual(len({row["direction_id"] for row in self.rows}), 5_761)
        self.assertEqual(direction_digest(self.rows), population["remainder_direction_sha256"])
        self.assertTrue(
            all(
                (left["projective_height"], left["quartic_x"], left["direction_id"])
                <= (right["projective_height"], right["quartic_x"], right["direction_id"])
                for left, right in zip(self.rows, self.rows[1:])
            )
        )
        selected = json.loads(SELECTED_ARTIFACT.read_text())["selected_directions"]
        selected_ids = {row["direction_id"] for row in selected}
        remainder_ids = {row["direction_id"] for row in self.rows}
        self.assertEqual(len(selected_ids), 399)
        self.assertTrue(selected_ids.isdisjoint(remainder_ids))
        self.assertEqual(len(selected_ids | remainder_ids), 6_160)
        self.assertEqual(
            direction_digest([*selected, *self.rows]),
            "4673b556e0a60943d86b23c7f293bfc3a9952f6acd25ef5e152297e84a106455",
        )
        for row in self.rows:
            vector = row["coefficient_vector"]
            self.assertEqual(sum(abs(value) for value in vector), 3)
            self.assertEqual(next(value for value in vector if value), 1)
            self.assertTrue(row["exact_auxiliary_inverse_checked"])
            self.assertTrue(row["exact_short_group_combination_checked"])

    def test_all_11522_slices_completed_once_with_only_calibration(self) -> None:
        self.assertEqual(len(self.slices), 11_522)
        self.assertEqual(
            Counter(row["slope"] for row in self.slices),
            Counter({-1: 5_761, 1: 5_761}),
        )
        self.assertEqual(
            Counter(row["search"]["status"] for row in self.slices),
            Counter({"completed": 11_522}),
        )
        self.assertTrue(all(not row["search"]["retried"] for row in self.slices))
        self.assertTrue(
            all(row["record_T0_calibration_count"] == 1 for row in self.slices)
        )
        self.assertTrue(
            all(row["qualifying_new_parameter_count"] == 0 for row in self.slices)
        )
        self.assertEqual(
            Counter(row["search"]["signed_point_count"] for row in self.slices),
            Counter({2: 11_522}),
        )
        self.assertEqual(
            Counter(len(row["incidences"]) for row in self.slices),
            Counter({1: 11_522}),
        )
        self.assertTrue(
            all(
                row["incidences"][0]["classification"]
                == "record-source-calibration-excluded"
                for row in self.slices
            )
        )
        self.assertEqual(
            self.data["stream"]["exact_slice_result_sha256"],
            "6818781aed81ca1bfd6822e3f276190d79b42da48b6762f3be9d5638500cb37e",
        )

    def test_current_auxiliary_artifact_replays_the_stable_parameter_set(self) -> None:
        prior = self.data["prior_parameter_decontamination"]
        predecessor = prior["predecessor_parameter_manifest"]
        self.assertEqual(
            predecessor["auxiliary_orbit_artifact_sha256_observed"],
            sha256(AUXILIARY_ARTIFACT),
        )
        self.assertFalse(
            predecessor["auxiliary_orbit_artifact_sha256_used_as_replay_gate"]
        )
        self.assertEqual(predecessor["auxiliary_orbit_extracted_parameter_count"], 666)
        self.assertEqual(
            predecessor["auxiliary_orbit_extracted_parameter_sha256"],
            "22e296780827d722bed88acb678323349b8f68abf7b298603901dbfdb49a8be1",
        )
        sys.path.insert(0, str(CAS))
        from search_fermigier_published_pair_fiber_products import (
            extract_parameter_values,
            rational_digest,
        )

        values = extract_parameter_values(json.loads(AUXILIARY_ARTIFACT.read_text()))
        self.assertEqual(len(values), 666)
        self.assertEqual(
            rational_digest(sorted(values)),
            predecessor["auxiliary_orbit_extracted_parameter_sha256"],
        )
        self.assertEqual(predecessor["terminal_prior_parameter_count"], 1_239)
        self.assertEqual(
            predecessor["terminal_prior_parameter_sha256"],
            "9482e61650aa8bb1fd45c3765e5db92c1474090faee8d831e0d73cee4fc864c4",
        )
        self.assertEqual(prior["terminal_prior_parameter_count"], 1_240)
        self.assertEqual(
            prior["terminal_prior_parameter_sha256"],
            "ca5be373dbe7934de0f7ee680203d6ae799f57fb55f712e73ad3777947b20cdd",
        )

    def test_clean_exhaustive_negative_handoff(self) -> None:
        outcome = self.data["outcome"]
        self.assertTrue(outcome["full_remainder_exhausted"])
        self.assertEqual(outcome["completed_direction_count"], 5_761)
        self.assertEqual(outcome["open_remainder_direction_count"], 0)
        self.assertFalse(outcome["open_remainder_claimed_negative"])
        self.assertEqual(outcome["slice_calls_attempted"], 11_522)
        self.assertEqual(outcome["slice_calls_completed"], 11_522)
        self.assertEqual(outcome["slice_calls_timed_out_or_errored"], 0)
        self.assertEqual(outcome["record_T0_calibrated_slices"], 11_522)
        self.assertEqual(
            outcome["incidence_classification_counts"],
            {"record-source-calibration-excluded": 11_522},
        )
        self.assertEqual(outcome["genuinely_new_forced_fibres"], 0)
        self.assertEqual(outcome["completed_conductors"], 0)
        self.assertEqual(outcome["subtarget_conductors"], 0)
        self.assertEqual(outcome["rank_triage_count"], 0)
        self.assertIsNone(outcome["maximum_stable_numerical_rank"])
        self.assertEqual(self.data["candidates"], [])
        self.assertFalse(self.data["target"]["hit"])
        self.assertFalse(self.data["execution"]["wall_cap_triggered"])
        self.assertEqual(self.data["execution"]["owned_processes_remaining"], 0)
        self.assertTrue(self.data["parameters"]["no_retries"])


if __name__ == "__main__":
    unittest.main()

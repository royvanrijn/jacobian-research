#!/usr/bin/env python3
"""Focused checks for the selected exact weight-four record-group tranche."""

from __future__ import annotations

from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
SCRIPT = CAS / "search_fermigier_rank22_record_group_quads.py"
ARTIFACT = GENERATED / "elliptic_fermigier_rank22_record_group_quads.json"
STREAM = GENERATED / "elliptic_fermigier_rank22_record_group_quads_stream.jsonl"
EXPECTED_SCRIPT_SHA256 = (
    "96213612a758efb60961bcaf5ba9882ffc51cf8c78f19c5bc704a3ddd32c972a"
)
EXPECTED_ARTIFACT_SHA256 = (
    "ef9a78aca5db699a7fd724599365690186b6105d5567b2a5956690a9de3eea2b"
)
EXPECTED_STREAM_SHA256 = (
    "b23514bacd6660dc6c9829bdc117e5467c245e7c8ede8d6108ff2d835b5a47c6"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direction_digest(records: list[dict]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda row: row["direction_id"]):
        digest.update(
            (
                f"{record['direction_id']}|{record['quartic_x']}|"
                f"{record['quartic_z']}|"
                f"{','.join(map(str, record['coefficient_vector']))}\n"
            ).encode()
        )
    return digest.hexdigest()


class FermigierRank22RecordGroupQuadsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.rows = [json.loads(line) for line in STREAM.read_text().splitlines()]
        cls.slices = [
            row for direction in cls.rows for row in direction["slice_searches"]
        ]

    def test_pinned_files_and_execution_provenance(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(sha256(STREAM), EXPECTED_STREAM_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(self.data["stream"]["sha256"], EXPECTED_STREAM_SHA256)
        self.assertEqual(
            self.data["stream"]["exact_slice_result_sha256"],
            "ab043bff4621955ebbb42532b822450691f68526d4f9596b164080f24daaa709",
        )
        self.assertEqual(self.data["execution"]["owned_processes_remaining"], 0)
        self.assertFalse(self.data["execution"]["wall_cap_triggered"])
        self.assertTrue(self.data["parameters"]["no_retries"])
        self.assertEqual(self.data["parameters"]["parallel_workers"], 4)

    def test_full_exact_weight4_population_and_prior_exclusion(self) -> None:
        population = self.data["full_weight4_population"]
        self.assertTrue(population["global_sign_quotient"])
        self.assertEqual(population["coefficient_alphabet"], [-1, 0, 1])
        self.assertEqual(population["exact_l1_norm"], 4)
        self.assertEqual(population["full_vector_count"], 58_520)
        self.assertEqual(population["exceptional_inverse_count"], 0)
        self.assertEqual(population["duplicate_new_abscissas"], 0)
        self.assertEqual(population["prior_quartic_abscissa_exclusions"], 1)
        self.assertEqual(population["genuinely_new_unique_abscissa_count"], 58_519)
        self.assertEqual(
            population["full_vector_classified_direction_sha256"],
            "4cbf7b6fd718f16c3919fe82c5cf5d49ba7453cdb03e8bcb6f0b9f0033a71781",
        )
        self.assertEqual(
            population["genuinely_new_direction_sha256"],
            "88ae41d291aedb71ca1f2b3516b62d7c417a1cd854dcc1a8e272f058fa54f7ba",
        )
        self.assertEqual(population["minimum_new_abscissa_projective_height"], 296_064_404)
        self.assertEqual(
            population["prior_quartic_abscissa_exclusion_records"],
            [
                {
                    "direction_id": "p02_p06_m08_m12",
                    "quartic_x": "618754/195",
                    "quartic_z": "5092764629261876/494325",
                }
            ],
        )
        lower = self.data["lower_weight_direction_exclusion"]
        self.assertEqual(lower["full_lower_weight_x_count"], 6_654)
        self.assertEqual(
            lower["full_lower_weight_x_sha256"],
            "741591c009f0e3c8a09232e2e73a66c7155b38a0d43bc44d44e1a020fcac4796",
        )
        vectors = []
        for support in itertools.combinations(range(22), 4):
            for relative_signs in itertools.product((-1, 1), repeat=3):
                vector = [0] * 22
                vector[support[0]] = 1
                for index, sign in zip(
                    support[1:], relative_signs, strict=True
                ):
                    vector[index] = sign
                vectors.append(tuple(vector))
        self.assertEqual(len(vectors), 58_520)
        self.assertEqual(len(set(vectors)), 58_520)
        for vector in vectors:
            self.assertEqual(sum(abs(value) for value in vector), 4)
            self.assertEqual(next(value for value in vector if value), 1)

    def test_blind_selector_and_honest_lower_weight_calibration(self) -> None:
        selection = self.data["rank_and_conductor_blind_selection"]
        self.assertFalse(selection["selection_uses_weight4_specialized_rank"])
        self.assertFalse(selection["selection_uses_weight4_conductor"])
        self.assertFalse(selection["selection_uses_weight4_slice_outcomes"])
        self.assertEqual(selection["selected_direction_count"], 5_179)
        self.assertEqual(selection["unselected_direction_count"], 53_340)
        self.assertEqual(
            selection["modular_lookup_sha256"],
            "74c56623dc6c06ea59ef015647cc5a3939dcc92819ed28806383504c37d09a99",
        )
        self.assertEqual(
            selection["full_modular_profile_sha256"],
            "6a3dddb0c94c6b8e48136451178bbb96bc3c404506a96ead56ba8ca11ab8080a",
        )
        self.assertEqual(
            selection["selected_direction_sha256"],
            "f701b3780e34b1fcee7f613ea78672221c1b92187788a21b0b37a7bc55114b3d",
        )
        self.assertEqual(
            selection["selection_and_strata_sha256"],
            "ced62382651dd937ad5eef4439d72cc564343d17edf03a0636a89cb9c74f4733",
        )
        selected = self.data["selected_directions"]
        self.assertEqual(len(selected), 5_179)
        self.assertEqual(len({row["direction_id"] for row in selected}), 5_179)
        self.assertEqual(direction_digest(selected), selection["selected_direction_sha256"])
        self.assertTrue(all(row["selection_strata"] for row in selected))

        calibration = self.data["lower_weight_selector_calibration"]
        self.assertEqual(calibration["fully_searched_weight2_population"], 462)
        self.assertEqual(calibration["selector_weight2_retained_count"], 78)
        self.assertEqual(calibration["known_weight2_second_fibre_count"], 0)
        self.assertEqual(calibration["fully_searched_weight3_population"], 6_160)
        self.assertEqual(calibration["selector_weight3_retained_count"], 560)
        self.assertEqual(calibration["known_weight3_second_fibre_count"], 1)
        self.assertTrue(calibration["known_hit_retained"])
        self.assertEqual(calibration["known_hit_height_rank_one_based"], 1)
        self.assertEqual(
            calibration["known_hit_descending_modular_rank_one_based"], 5_502
        )
        self.assertIn("in-sample", calibration["selection_leakage"])

    def test_all_selected_slices_completed_once_with_only_calibration(self) -> None:
        self.assertEqual(len(self.rows), 5_179)
        self.assertEqual(direction_digest(self.rows), "f701b3780e34b1fcee7f613ea78672221c1b92187788a21b0b37a7bc55114b3d")
        self.assertEqual(len(self.slices), 10_358)
        self.assertEqual(
            Counter(row["slope"] for row in self.slices),
            Counter({-1: 5_179, 1: 5_179}),
        )
        self.assertEqual(
            Counter(row["search"]["status"] for row in self.slices),
            Counter({"completed": 10_358}),
        )
        self.assertEqual(
            Counter(row["search"]["signed_point_count"] for row in self.slices),
            Counter({2: 10_358}),
        )
        self.assertTrue(all(not row["search"]["retried"] for row in self.slices))
        self.assertTrue(
            all(row["record_T0_calibration_count"] == 1 for row in self.slices)
        )
        self.assertTrue(
            all(row["qualifying_new_parameter_count"] == 0 for row in self.slices)
        )
        self.assertTrue(
            all(
                len(row["incidences"]) == 1
                and row["incidences"][0]["classification"]
                == "record-source-calibration-excluded"
                for row in self.slices
            )
        )

    def test_clean_negative_handoff_and_complete_prior_decontamination(self) -> None:
        outcome = self.data["outcome"]
        self.assertTrue(outcome["selected_tranche_exhausted"])
        self.assertEqual(outcome["completed_selected_direction_count"], 5_179)
        self.assertEqual(outcome["open_selected_direction_count"], 0)
        self.assertEqual(outcome["slice_calls_attempted"], 10_358)
        self.assertEqual(outcome["slice_calls_completed"], 10_358)
        self.assertEqual(outcome["slice_calls_timed_out_or_errored"], 0)
        self.assertEqual(outcome["record_T0_calibrated_slices"], 10_358)
        self.assertEqual(
            outcome["incidence_classification_counts"],
            {"record-source-calibration-excluded": 10_358},
        )
        self.assertEqual(outcome["genuinely_new_forced_fibres"], 0)
        self.assertEqual(outcome["completed_conductors"], 0)
        self.assertEqual(outcome["rank_triage_count"], 0)
        self.assertIsNone(outcome["maximum_stable_numerical_rank"])
        self.assertEqual(outcome["unselected_direction_count"], 53_340)
        self.assertFalse(outcome["unselected_directions_claimed_negative"])
        self.assertEqual(self.data["candidates"], [])
        self.assertFalse(self.data["target"]["hit"])
        prior = self.data["prior_parameter_decontamination"]
        self.assertEqual(prior["unique_prior_parameter_count"], 1_686)
        self.assertEqual(
            prior["prior_parameter_sha256"],
            "f2e060657e16e4b0b57f3ae210afba8c5147b69ed23c496416f2899339da3548",
        )
        self.assertEqual(len(prior["artifact_sources"]), 28)


if __name__ == "__main__":
    unittest.main()

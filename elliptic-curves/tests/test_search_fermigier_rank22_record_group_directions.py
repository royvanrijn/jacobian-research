#!/usr/bin/env python3
"""Focused checks for the rank-22 record-quartic group-direction search."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))

from search_fermigier_rank22_record_group_directions import (  # noqa: E402
    RecordQuarticAuxiliary,
)


SCRIPT = CAS / "search_fermigier_rank22_record_group_directions.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_rank22_record_group_directions.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "921aeb1eb1193d9f4dba98d7b381d17a28004e15ff114bc06919c42f6a622e28"
)
EXPECTED_ARTIFACT_SHA256 = (
    "4928c44e27cada74b7a558dd97edfba20a554a508ac0e98ca051df7dea66a3c1"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierRank22RecordGroupDirectionsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.slice_rows = [
            row
            for direction in cls.data["direction_searches"]
            for row in direction["slice_searches"]
        ]

    def test_pinned_files_and_exact_rank22_source(self) -> None:
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
        known = source["known_record_abscissas"]
        self.assertEqual(known["H1000000_abscissa_count"], 27)
        self.assertEqual(known["published_preimage_abscissa_count"], 22)
        self.assertEqual(known["union_count"], 32)
        self.assertEqual(
            known["union_sha256"],
            "a692902016852128dac0997705d5463a4d94f48a0c4089d324f0c110f8ac8167",
        )

    def test_explicit_birational_map_round_trips_the_basis(self) -> None:
        map_record = self.data["pointed_quartic_map"]
        self.assertEqual(map_record["record_parameter_t"], "39508/39")
        self.assertEqual(
            map_record["origin_quartic_point"],
            {"x": "-39508/39", "z": "1293950112938/507"},
        )
        auxiliary = RecordQuarticAuxiliary.construct()
        self.assertEqual(
            [str(value) for value in auxiliary.weierstrass_coefficients],
            map_record["generalized_weierstrass_coefficients_a1_a2_a3_a4_a6"],
        )
        basis = self.data["transported_rank22_basis"]
        self.assertEqual(len(basis), 22)
        self.assertEqual([row["label"] for row in basis], [f"P{i}" for i in range(1, 23)])
        for row in basis:
            quartic_point = (
                Fraction(row["quartic_preimage"]["x"]),
                Fraction(row["quartic_preimage"]["z"]),
            )
            auxiliary_point = (
                Fraction(row["auxiliary_group_point"]["x"]),
                Fraction(row["auxiliary_group_point"]["y"]),
            )
            self.assertEqual(auxiliary.forward(quartic_point), auxiliary_point)
            self.assertEqual(auxiliary.inverse(auxiliary_point), quartic_point)
            self.assertTrue(row["exact_short_group_coordinate_checked"])

    def test_sign_quotiented_orbit_manufactures_462_new_directions(self) -> None:
        orbit = self.data["orbit_generation"]
        self.assertTrue(orbit["global_sign_quotient"])
        self.assertEqual(orbit["coefficient_alphabet"], [-1, 0, 1])
        self.assertEqual(orbit["maximum_l1_norm"], 2)
        self.assertEqual(orbit["sign_quotiented_vector_count"], 484)
        self.assertEqual(orbit["singleton_transport_calibrations"], 22)
        self.assertEqual(orbit["known_abscissa_exclusions"], 22)
        self.assertEqual(orbit["exceptional_inverse_count"], 0)
        self.assertEqual(orbit["unique_orbit_abscissa_count"], 484)
        self.assertEqual(orbit["genuinely_new_pair_combination_abscissa_count"], 462)
        self.assertEqual(orbit["minimum_new_abscissa_projective_height"], 60_479_122)
        self.assertEqual(
            orbit["maximum_new_abscissa_projective_height"],
            1_535_173_523_922_167_228_392_881_239_832_529_034,
        )
        self.assertEqual(
            orbit["new_direction_sha256"],
            "2819df851b876fe430465e2d7fbb838c34b60d597aea2ee48c8d912b07c19939",
        )
        directions = self.data["direction_searches"]
        self.assertEqual(len(directions), 462)
        self.assertEqual(len({row["direction_id"] for row in directions}), 462)
        self.assertEqual(len({row["quartic_x"] for row in directions}), 462)
        for row in directions:
            vector = row["coefficient_vector"]
            self.assertEqual(len(vector), 22)
            self.assertEqual(sum(abs(value) for value in vector), 2)
            self.assertEqual(next(value for value in vector if value), 1)
            self.assertTrue(row["exact_auxiliary_inverse_checked"])
            self.assertTrue(row["exact_short_group_combination_checked"])

    def test_every_prior_parameter_is_pinned_before_search(self) -> None:
        prior = self.data["prior_decontamination"]
        self.assertEqual(prior["base_manifest_parameter_count"], 590)
        self.assertEqual(
            prior["base_manifest_parameter_sha256"],
            "64c09a13b427938a44251a91f74a116f7f9e685aed07c6159550e7ec3ea51291",
        )
        self.assertEqual(prior["auxiliary_orbit_extracted_parameter_count"], 666)
        self.assertEqual(
            prior["H50000_seen_parameters"],
            ["39508/39", "42058/25", "48363/26", "23317/6"],
        )
        self.assertEqual(prior["terminal_prior_parameter_count"], 1_239)
        self.assertEqual(
            prior["terminal_prior_parameter_sha256"],
            "9482e61650aa8bb1fd45c3765e5db92c1474090faee8d831e0d73cee4fc864c4",
        )

    def test_all_924_slices_complete_with_only_the_source_calibration(self) -> None:
        self.assertEqual(len(self.slice_rows), 924)
        self.assertEqual(
            Counter(row["search"]["status"] for row in self.slice_rows),
            Counter({"completed": 924}),
        )
        self.assertEqual(
            Counter(row["search"]["signed_point_count"] for row in self.slice_rows),
            Counter({2: 924}),
        )
        self.assertTrue(all(not row["search"]["retried"] for row in self.slice_rows))
        self.assertTrue(
            all(row["record_T0_calibration_count"] == 1 for row in self.slice_rows)
        )
        self.assertTrue(all(len(row["incidences"]) == 1 for row in self.slice_rows))
        self.assertTrue(
            all(
                row["incidences"][0]["classification"]
                == "record-source-calibration-excluded"
                for row in self.slice_rows
            )
        )
        outcome = self.data["outcome"]
        self.assertEqual(outcome["slice_calls_attempted"], 924)
        self.assertEqual(outcome["slice_calls_completed"], 924)
        self.assertEqual(outcome["slice_calls_timed_out_or_errored"], 0)
        self.assertEqual(outcome["record_T0_calibrated_slices"], 924)
        self.assertEqual(
            outcome["incidence_classification_counts"],
            {"record-source-calibration-excluded": 924},
        )
        self.assertEqual(
            outcome["exact_slice_result_sha256"],
            "867ebd2e8c2b1a3608ea5577dbaa356de6c1b8e2574c45290fcbbbbf2ea4baf5",
        )

    def test_negative_target_scope_and_clean_process_exit(self) -> None:
        self.assertEqual(self.data["candidates"], [])
        outcome = self.data["outcome"]
        self.assertEqual(outcome["genuinely_new_forced_fibres"], 0)
        self.assertEqual(outcome["completed_conductors"], 0)
        self.assertEqual(outcome["subtarget_conductors"], 0)
        self.assertEqual(outcome["rank_triage_count"], 0)
        self.assertIsNone(outcome["maximum_stable_numerical_rank"])
        self.assertFalse(self.data["target"]["hit"])
        self.assertTrue(self.data["parameters"]["no_retries"])
        self.assertEqual(self.data["execution"]["owned_processes_remaining"], 0)


if __name__ == "__main__":
    unittest.main()

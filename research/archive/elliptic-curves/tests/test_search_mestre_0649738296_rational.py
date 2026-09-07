#!/usr/bin/env python3
"""Exact replay tests for the rank-14 Mestre-family rational scan."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from search_mestre_root_tuple_scale import (  # noqa: E402
    finite_reduction_attempt,
    point_digest,
    point_on_short_curve,
)
from search_mestre_0649738296_rational import (  # noqa: E402
    CALIBRATION_PARAMETER,
    DISCOVERY_PRIMES,
    FROZEN_COMPLETE_ARTIFACT_SHA256,
    HELD_PRIMES,
    ROOTS,
    SIGNAL_CERTIFICATE_PRIME_BOUND,
    candidate_record,
    discriminant_features,
    exact_discriminant,
    family_coefficients,
    result_digest,
    run_scanner,
    select_population,
    stable_digest,
)


SCRIPT = CAS / "search_mestre_0649738296_rational.py"
SCANNER = CAS / "scan_mestre_0649738296_rational.cpp"
ARTIFACT = GENERATED / "elliptic_mestre_0649738296_rational.json"
COMPLETE_ARTIFACT = (
    GENERATED / "elliptic_mestre_root_tuple_scale_max100_complete.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "16eb2bfd0e43afc2fb801dc6bf319c0f604ae060230eafe12f4c8b93c0b2ac8c"
)
EXPECTED_SCANNER_SHA256 = (
    "af3d4c09c78115a6962edc77968d04c0bd2fdc40dbb3ff3e0767c6386e64eb8f"
)
EXPECTED_ARTIFACT_SHA256 = (
    "1c0e001c6c03e557722e16897f66ad2c90c93aa7d88f7cbbdd286700c66eaa78"
)
EXPECTED_SCANNER_STDOUT_SHA256 = (
    "ccaa84c4b52719e960fe9ad366ed6d1c068d22d46138f4db86a6e50194a4f448"
)
EXPECTED_FEATURE_SHA256 = (
    "dd44e7d33e7b8a4fbc8000ab7649352b510eb0257a09105243073bc7f7df17ae"
)
EXPECTED_SELECTED_SHA256 = (
    "24cd3d511403d82de1cc24dea2e45e4bfc1f3bb3f314c13b67b93f3b2b98730f"
)
CALIBRATION_BASIS_SHA256 = (
    "e177eb60564f0c4d80f862c501409a60672f66603e673efa9a5b7e6fed83cb10"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Mestre0649738296RationalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.records = cls.data["conductor_first_screen"]["records"]
        cls.by_identifier = {record["identifier"]: record for record in cls.records}

    def test_files_and_frozen_rank14_calibration_are_pinned(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(SCANNER), EXPECTED_SCANNER_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(sha256(COMPLETE_ARTIFACT), FROZEN_COMPLETE_ARTIFACT_SHA256)
        provenance = self.data["provenance"]
        self.assertEqual(provenance["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(provenance["scanner_sha256"], EXPECTED_SCANNER_SHA256)
        self.assertTrue(provenance["all_external_processes_foreground_and_capped"])
        self.assertTrue(provenance["whole_process_groups_killed_and_reaped_on_timeout"])
        self.assertTrue(provenance["no_retries"])
        self.assertEqual(provenance["owned_processes_remaining"], 0)

    def test_exact_even_family_and_factored_discriminant_replay(self) -> None:
        construction = SixRootMestreConstruction(
            tuple(Fraction(root) for root in ROOTS)
        )
        self.assertEqual(construction.quartic_condition, 0)
        self.assertFalse(construction.is_reflection_symmetric)
        self.assertEqual(construction.quartic_square_scale, 8)
        for parameter in (
            Fraction(1),
            Fraction(2),
            Fraction(7, 3),
            Fraction(29, 11),
            Fraction(-13, 7),
        ):
            coefficients = construction.primitive_jacobian_coefficients(parameter)
            self.assertEqual(coefficients, family_coefficients(parameter))
            self.assertEqual(coefficients, family_coefficients(-parameter))
            self.assertEqual(
                construction.primitive_quartic_coefficients(parameter),
                construction.primitive_quartic_coefficients(-parameter),
            )
            self.assertEqual(
                construction.primitive_quartic_discriminant(parameter),
                exact_discriminant(parameter),
            )
        self.assertEqual(exact_discriminant(Fraction(45, 2)), 0)
        family = self.data["family"]
        self.assertEqual(tuple(family["roots"]), ROOTS)
        self.assertEqual(family["primitive_quartic_square_scale"], 8)
        self.assertIn("T=a/b", family["search_quotient"])

    def test_full_reduced_scanner_and_exact_feature_selection_replay(self) -> None:
        scan = run_scanner(
            SCANNER,
            numerator_bound=30_000,
            denominator_bound=1_000,
            keep=8_192,
            compile_timeout=30,
            scan_timeout=30,
        )
        self.assertEqual(scan.primitive_population, 18_244_819)
        self.assertEqual(scan.prior_excluded, 8)
        self.assertEqual(scan.evaluated_population, 18_244_811)
        self.assertEqual(len(scan.candidates), 8_192)
        self.assertEqual(scan.stdout_sha256, EXPECTED_SCANNER_STDOUT_SHA256)
        self.assertEqual(tuple(self.data["modular_scan"]["discovery_primes"]), DISCOVERY_PRIMES)
        self.assertEqual(tuple(self.data["modular_scan"]["held_primes"]), HELD_PRIMES)
        self.assertFalse(set(DISCOVERY_PRIMES) & set(HELD_PRIMES))
        self.assertEqual(
            scan.discovery_table_digest,
            self.data["modular_scan"]["discovery_table_digest"],
        )
        self.assertEqual(
            scan.held_table_digest,
            self.data["modular_scan"]["held_table_digest"],
        )

        features = {
            candidate.identifier: discriminant_features(candidate)
            for candidate in scan.candidates
        }
        feature_records = [
            candidate_record(candidate, features[candidate.identifier])
            for candidate in scan.candidates
        ]
        self.assertEqual(stable_digest(feature_records), EXPECTED_FEATURE_SHA256)
        self.assertEqual(
            self.data["modular_scan"]["retained_candidate_feature_sha256"],
            EXPECTED_FEATURE_SHA256,
        )
        selection = select_population(scan.candidates, features)
        selected_ids = sorted(candidate.identifier for candidate, _ in selection)
        self.assertEqual(len(selection), 208)
        self.assertEqual(
            hashlib.sha256("\n".join(selected_ids).encode()).hexdigest(),
            EXPECTED_SELECTED_SHA256,
        )
        self.assertEqual(set(selected_ids), set(self.by_identifier))

    def test_T2_exact_rank14_calibration_replays_and_is_excluded(self) -> None:
        calibration = self.data["calibration_T2"]
        self.assertEqual(Fraction(calibration["parameter"]), CALIBRATION_PARAMETER)
        self.assertTrue(calibration["excluded_from_selection"])
        self.assertNotIn("t2_1", self.by_identifier)
        self.assertEqual(calibration["basis_point_sha256"], CALIBRATION_BASIS_SHA256)
        complete = json.loads(COMPLETE_ARTIFACT.read_text())
        source = next(
            record
            for record in complete["fiber_records"]
            if record["identifier"] == "r0_6_49_73_82_96_t2"
        )
        basis = tuple(
            (Fraction(point["jacobian_x"]), Fraction(point["jacobian_y"]))
            for point in source["H5000_phase"]["small_prime_saturation"][
                "saturated_basis"
            ]
        )
        coefficients = family_coefficients(CALIBRATION_PARAMETER)
        self.assertEqual(point_digest(basis), CALIBRATION_BASIS_SHA256)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in basis))
        replay = finite_reduction_attempt(
            coefficients,
            basis,
            prime_bound=SIGNAL_CERTIFICATE_PRIME_BOUND,
        )
        self.assertEqual(replay, calibration["finite_reduction_certificate"])
        self.assertEqual(replay["certified_algebraic_rank_lower_bound"], 14)
        self.assertEqual(replay["combined_exact_rank_over_F2"], 14)
        self.assertEqual(replay["two_torsion_certificate_prime"], 19)

    def test_conductor_closure_and_all_H5000_calls_are_accounted(self) -> None:
        population = self.data["conductor_first_screen"]["population"]
        self.assertTrue(
            self.data["conductor_first_screen"][
                "population_closed_before_any_point_or_rank_call"
            ]
        )
        self.assertEqual(population["selected"], 208)
        self.assertEqual(population["completed"], 197)
        self.assertEqual(population["subtarget"], 62)
        self.assertEqual(
            population["status_histogram"],
            {
                "completed exact PARI minimal-model/conductor computation": 197,
                "timeout-no-retry": 11,
            },
        )
        self.assertEqual(len(self.records), 208)
        completed = [
            record
            for record in self.records
            if record["conductor_phase"]["status"].startswith("completed")
        ]
        self.assertEqual(len(completed), 197)
        self.assertTrue(
            all("H5000" in record.get("point_stages", {}) for record in completed)
        )
        self.assertTrue(
            all(
                record["point_stages"]["H5000"]["status"].startswith("completed")
                for record in completed
            )
        )
        self.assertEqual(
            sum(
                Decimal(record["conductor_phase"]["log_conductor"])
                < Decimal("182.72")
                for record in completed
            ),
            62,
        )

    def test_staged_rank_frontier_is_closed_without_broadening(self) -> None:
        protocol = self.data["point_search_protocol"]
        self.assertEqual(
            protocol["completed_stage_calls"],
            {"H5000": 197, "H50000": 32, "H250000": 8},
        )
        self.assertEqual(
            protocol["stage_status_histograms"],
            {
                "H5000": {
                    "completed exact point membership and numerical height triage": 197
                },
                "H50000": {
                    "completed exact point membership and numerical height triage": 32
                },
                "H250000": {
                    "completed exact point membership and numerical height triage": 8
                },
                "H1000000": {"timeout-no-retry": 2},
            },
        )
        self.assertEqual(
            protocol["stable_numerical_rank_histograms"],
            {
                "H5000": {"9": 1, "10": 194, "11": 2},
                "H50000": {"10": 30, "11": 2},
                "H250000": {"10": 5, "11": 3},
                "H1000000": {},
            },
        )
        self.assertEqual(protocol["maximum_stable_numerical_rank"], 11)
        self.assertFalse(protocol["adaptive_broadening_performed"])
        self.assertEqual(
            protocol["finite_reduction_trigger_stable_numerical_rank"], 18
        )
        self.assertEqual(self.data["target"]["hits"], [])
        leaders = self.data["numerical_leaders"]
        self.assertEqual(leaders[0]["identifier"], "t2745_38")
        self.assertEqual(leaders[0]["stable_numerical_rank"], 11)
        self.assertEqual(leaders[0]["stage"], "H250000")
        self.assertEqual(result_digest(self.data), self.data["result_sha256"])


if __name__ == "__main__":
    unittest.main()

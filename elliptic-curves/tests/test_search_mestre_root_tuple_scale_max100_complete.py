#!/usr/bin/env python3
"""Focused exact checks for the complete max-root-100 Mestre fiber panel."""

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
    tuple_digest,
)
from search_mestre_root_tuple_scale_max100_complete import (  # noqa: E402
    PARAMETERS,
    SIGNAL_CERTIFICATE_PRIME_BOUND,
    fiber_identifier,
    identifier_digest,
    result_digest,
)


SCRIPT = CAS / "search_mestre_root_tuple_scale_max100_complete.py"
ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale_max100_complete.json"
MAX50_ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale.json"
MAX100_ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale_max100.json"
EXPECTED_SCRIPT_SHA256 = (
    "5fefa372dca3563bedef08a55b1adad8c7db26c5161db38faa89a2e12d59ef6c"
)
EXPECTED_ARTIFACT_SHA256 = (
    "c2cb2e68a54cc1625224cf1aed1cce0196582c16d21373bdbc86c67a1e91d24c"
)
EXPECTED_FAMILY_SHA256 = (
    "e92e9cd0be8fc8006275797df2752b714df0237ae27ce2b3ba4829c988681973"
)
EXPECTED_ADMISSIBLE_SHA256 = (
    "8761c92ea0dd703fb3dee53ffed330429c54755edcaa6af452b00b3ed99bcc7d"
)
EXPECTED_INADMISSIBLE_SHA256 = (
    "fcdbaa833b715198b328ac14910b2886ad8023907540608fcb410f23fcdacaf1"
)
LEADER = "r0_6_49_73_82_96_t2"
LEADER_ROOTS = (0, 6, 49, 73, 82, 96)
LEADER_PARAMETER = 2
LEADER_CONDUCTOR = "220490760420432505357569342958539032520287541627870"
LEADER_LOG_CONDUCTOR = (
    "115.919940254867962972869938621968140350952714539550596915297"
)
LEADER_BASIS_SHA256 = (
    "e177eb60564f0c4d80f862c501409a60672f66603e673efa9a5b7e6fed83cb10"
)
LEADER_CERTIFICATE_PRIMES = (13, 29, 37, 47, 53, 59, 61, 67, 71, 73, 107)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRootTupleScaleMax100CompleteTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.max50 = json.loads(MAX50_ARTIFACT.read_text())
        cls.max100 = json.loads(MAX100_ARTIFACT.read_text())
        cls.families = tuple(tuple(roots) for roots in cls.data["scope"]["families"])
        cls.records = cls.data["fiber_records"]
        cls.by_identifier = {record["identifier"]: record for record in cls.records}

    def test_script_artifact_and_frozen_inputs_are_pinned(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            self.data["provenance"]["script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        self.assertEqual(
            self.data["frozen_inputs"]["max50_artifact_sha256"],
            sha256(MAX50_ARTIFACT),
        )
        self.assertEqual(
            self.data["frozen_inputs"]["max100_artifact_sha256"],
            sha256(MAX100_ARTIFACT),
        )
        self.assertTrue(self.data["frozen_inputs"]["all_frozen_files_read_only"])
        self.assertTrue(
            self.data["provenance"]["all_external_processes_foreground_and_capped"]
        )
        self.assertTrue(self.data["provenance"]["no_retries"])
        self.assertEqual(self.data["provenance"]["owned_processes_remaining"], 0)

    def test_complete_family_and_exact_admissible_populations_replay(self) -> None:
        scope = self.data["scope"]
        population = self.data["population"]
        self.assertEqual(len(self.families), 235)
        self.assertEqual(scope["complete_family_count"], 235)
        self.assertEqual(tuple_digest(self.families), EXPECTED_FAMILY_SHA256)
        self.assertEqual(scope["family_sha256"], EXPECTED_FAMILY_SHA256)
        self.assertEqual(scope["integer_parameters"], list(PARAMETERS))
        self.assertEqual(scope["proposed_fiber_count"], 1_880)
        self.assertEqual(
            scope["proposed_fiber_identifier_sha256"],
            identifier_digest(
                fiber_identifier(roots, parameter)
                for roots in self.families
                for parameter in PARAMETERS
            ),
        )

        admissible = []
        inadmissible = []
        for roots in self.families:
            construction = SixRootMestreConstruction(
                tuple(Fraction(root) for root in roots)
            )
            for parameter in PARAMETERS:
                parameter_q = Fraction(parameter)
                degeneracy = construction.visible_point_degeneracy(parameter_q)
                identifier = fiber_identifier(roots, parameter)
                is_admissible = (
                    construction.quartic_discriminant(parameter_q) != 0
                    and degeneracy.collision_loss == 0
                    and degeneracy.zero_ordinates == 0
                )
                (admissible if is_admissible else inadmissible).append(identifier)
        self.assertEqual(len(admissible), 1_426)
        self.assertEqual(len(inadmissible), 454)
        self.assertEqual(identifier_digest(admissible), EXPECTED_ADMISSIBLE_SHA256)
        self.assertEqual(identifier_digest(inadmissible), EXPECTED_INADMISSIBLE_SHA256)
        self.assertEqual(population["admissible_fiber_count"], len(admissible))
        self.assertEqual(population["inadmissible_fiber_count"], len(inadmissible))
        self.assertEqual(
            population["admissible_fiber_identifier_sha256"],
            EXPECTED_ADMISSIBLE_SHA256,
        )
        self.assertEqual(
            population["inadmissible_fiber_identifier_sha256"],
            EXPECTED_INADMISSIBLE_SHA256,
        )
        self.assertEqual(len(self.records), 1_426)
        self.assertEqual(len(self.by_identifier), 1_426)
        self.assertEqual(set(self.by_identifier), set(admissible))

    def test_phase_specific_frozen_exclusions_replay_exactly(self) -> None:
        exclusions = self.data["exact_phase_exclusions"]
        max50_records = self.max50["specialization_screen"][
            "conductor_first_fiber_records"
        ]
        max100_records = self.max100["specialization_screen"]["conductor_records"]
        max50_h5000 = self.max50["specialization_screen"]["point_search_finalists"]
        max100_h5000 = self.max100["specialization_screen"]["h5000_records"]
        populations = (
            (
                max50_records,
                202,
                "frozen_max50_visible_and_conductor_identifier_sha256",
            ),
            (max100_records, 262, "frozen_max100_conductor_identifier_sha256"),
            (max50_h5000, 25, "frozen_max50_H5000_identifier_sha256"),
            (max100_h5000, 64, "frozen_max100_H5000_identifier_sha256"),
        )
        for records, expected_count, digest_key in populations:
            identifiers = [record["identifier"] for record in records]
            self.assertEqual(len(identifiers), expected_count)
            self.assertEqual(len(set(identifiers)), expected_count)
            self.assertEqual(
                identifier_digest(identifiers), exclusions[digest_key]
            )
        self.assertEqual(exclusions["visible_fresh_count"], 1_224)
        self.assertEqual(exclusions["conductor_fresh_count"], 962)
        self.assertEqual(exclusions["H5000_fresh_count"], 1_337)
        self.assertTrue(exclusions["exclusions_are_phase_specific"])
        self.assertTrue(exclusions["no_phase_retry"])

    def test_every_conductor_and_H5000_phase_closed_below_target(self) -> None:
        population = self.data["population"]
        self.assertEqual(population["subtarget_conductor_count"], 1_426)
        self.assertTrue(
            all(
                record["conductor_phase"].get(
                    "below_strict_log_conductor_target_numerically"
                )
                and Decimal(record["conductor_phase"]["log_conductor"])
                < Decimal("182.72")
                for record in self.records
            )
        )
        self.assertEqual(
            population["conductor_status_histogram"],
            {
                "completed fresh exact PARI minimal-model/conductor computation": 962,
                "reused complete frozen max-root-100 conductor phase": 262,
                "reused complete frozen max-root-50 conductor phase": 202,
            },
        )
        self.assertEqual(
            population["H5000_status_histogram"],
            {
                "completed exact bounded point checks and numerical height triage": 1_337,
                "reused complete frozen max-root-100 H5000 point/height phase": 64,
                "reused complete frozen max-root-50 H5000 point/height phase": 25,
            },
        )
        self.assertEqual(
            population["H5000_stable_numerical_rank_histogram"],
            {
                "4": 5,
                "5": 23,
                "6": 94,
                "7": 164,
                "8": 175,
                "9": 304,
                "10": 381,
                "11": 211,
                "12": 61,
                "13": 7,
                "14": 1,
            },
        )
        self.assertEqual(population["maximum_H5000_stable_numerical_rank"], 14)
        self.assertEqual(population["maximum_any_stage_stable_numerical_rank"], 14)
        self.assertEqual(
            population["maximum_certified_algebraic_rank_lower_bound"], 14
        )
        self.assertEqual(population["distinct_certified_signal_count"], 1)
        self.assertEqual(self.data["target"]["hits"], [])

    def test_rank14_frontier_certificate_replays_exactly(self) -> None:
        record = self.by_identifier[LEADER]
        self.assertEqual(tuple(record["roots"]), LEADER_ROOTS)
        self.assertEqual(record["parameter"], LEADER_PARAMETER)
        conductor = record["conductor_phase"]
        self.assertEqual(conductor["conductor"], LEADER_CONDUCTOR)
        self.assertEqual(conductor["log_conductor"], LEADER_LOG_CONDUCTOR)
        self.assertEqual(conductor["root_number"], 1)
        phase = record["H5000_phase"]
        self.assertEqual(phase["stable_numerical_rank"], 14)
        self.assertEqual(phase["signed_points_returned"], 44)
        self.assertEqual(phase["distinct_nonzero_ordinate_abscissas"], 22)
        self.assertEqual(phase["pool_point_count_modulo_inverse"], 22)
        self.assertEqual(
            phase["small_prime_saturation"]["saturated_basis_sha256"],
            LEADER_BASIS_SHA256,
        )

        construction = SixRootMestreConstruction(
            tuple(Fraction(root) for root in LEADER_ROOTS)
        )
        coefficients = construction.primitive_jacobian_coefficients(
            Fraction(LEADER_PARAMETER)
        )
        basis = tuple(
            (Fraction(point["jacobian_x"]), Fraction(point["jacobian_y"]))
            for point in phase["small_prime_saturation"]["saturated_basis"]
        )
        self.assertEqual(len(basis), 14)
        self.assertEqual(point_digest(basis), LEADER_BASIS_SHA256)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in basis))
        replay = finite_reduction_attempt(
            coefficients,
            basis,
            prime_bound=SIGNAL_CERTIFICATE_PRIME_BOUND,
        )
        self.assertEqual(replay, phase["finite_reduction_attempt"])
        self.assertEqual(replay["status"], "certified")
        self.assertEqual(replay["certified_algebraic_rank_lower_bound"], 14)
        self.assertEqual(replay["combined_exact_rank_over_F2"], 14)
        self.assertEqual(
            tuple(replay["certificate_primes"]), LEADER_CERTIFICATE_PRIMES
        )
        self.assertEqual(replay["two_torsion_certificate_prime"], 19)
        self.assertEqual(
            self.data["certified_frontier"],
            [
                {
                    "basis_point_sha256": LEADER_BASIS_SHA256,
                    "certificate_primes": list(LEADER_CERTIFICATE_PRIMES),
                    "certified_algebraic_rank_lower_bound": 14,
                    "combined_exact_rank_over_F2": 14,
                    "conductor": conductor,
                    "identifier": LEADER,
                    "parameter": LEADER_PARAMETER,
                    "primitive_jacobian_coefficients": [str(value) for value in coefficients],
                    "roots": list(LEADER_ROOTS),
                    "stage": "H5000",
                    "two_torsion_certificate_prime": 19,
                }
            ],
        )

    def test_staged_escalation_and_internal_digest_are_closed(self) -> None:
        stages = self.data["escalation_records"]
        self.assertEqual({stage: len(records) for stage, records in stages.items()}, {
            "H50000": 16,
            "H250000": 4,
            "H1000000": 1,
        })
        for stage in ("H50000", "H250000"):
            leader = next(record for record in stages[stage] if record["identifier"] == LEADER)
            self.assertEqual(leader["phase"]["stable_numerical_rank"], 14)
            self.assertEqual(
                leader["phase"]["finite_reduction_attempt"][
                    "certified_algebraic_rank_lower_bound"
                ],
                14,
            )
        deepest = stages["H1000000"][0]
        self.assertEqual(deepest["identifier"], LEADER)
        self.assertEqual(deepest["input_stable_numerical_rank"], 14)
        self.assertEqual(deepest["phase"]["status"], "timeout-no-retry")
        self.assertEqual(deepest["phase"]["point_timeout_seconds"], 25.0)
        self.assertEqual(result_digest(self.data), self.data["result_sha256"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import sys
import unittest


ELLIPTIC_DIRECTORY = Path(__file__).resolve().parents[1]
CAS_DIRECTORY = ELLIPTIC_DIRECTORY / "cas"
REPO_ROOT = ELLIPTIC_DIRECTORY.parent
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from search_mestre_0430313946_frontier import (  # noqa: E402
    A_COEFFICIENTS,
    B_COEFFICIENTS,
    CALIBRATION_PARAMETER,
    DISCOVERY_PRIMES,
    FROZEN_SCALE_ARTIFACT_SHA256,
    FROZEN_SCALE_RESULT_SHA256,
    HELD_PRIMES,
    ROOTS,
    exact_local_table_digest,
    family_coefficients,
    point_on_short_curve,
    run_scanner_strata,
)
from search_mestre_root_tuple_scale import capped_minimal_curve_data  # noqa: E402


ARTIFACT = (
    REPO_ROOT
    / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_0430313946_frontier.json"
)
SCALE_ARTIFACT = (
    REPO_ROOT
    / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale.json"
)
SCRIPT = CAS_DIRECTORY / "search_mestre_0430313946_frontier.py"
SCANNER = CAS_DIRECTORY / "scan_mestre_0430313946.cpp"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Mestre0430313946GeometryTests(unittest.TestCase):
    def test_exact_even_polynomial_model_is_pinned(self) -> None:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
        self.assertEqual(construction.quartic_condition, 0)
        self.assertFalse(construction.is_reflection_symmetric)
        self.assertEqual(construction.quartic_square_scale, 40)
        self.assertEqual(len(A_COEFFICIENTS), 9)
        self.assertEqual(len(B_COEFFICIENTS), 13)
        self.assertTrue(all(A_COEFFICIENTS[index] == 0 for index in range(1, 9, 2)))
        self.assertTrue(all(B_COEFFICIENTS[index] == 0 for index in range(1, 13, 2)))
        for parameter in (Q(1), Q(5), Q(7, 3), Q(-11, 4)):
            self.assertEqual(
                construction.primitive_jacobian_coefficients(parameter),
                family_coefficients(parameter),
            )
            self.assertEqual(
                construction.primitive_quartic_coefficients(parameter),
                construction.primitive_quartic_coefficients(-parameter),
            )

    def test_exact_local_table_digests_are_pinned(self) -> None:
        self.assertFalse(set(DISCOVERY_PRIMES) & set(HELD_PRIMES))
        self.assertEqual(
            exact_local_table_digest(DISCOVERY_PRIMES), "5434355788148175162"
        )
        self.assertEqual(
            exact_local_table_digest(HELD_PRIMES), "15815491295577555910"
        )

    @unittest.skipUnless(shutil.which("c++"), "a C++17 compiler is required")
    def test_small_compiled_scanner_replays_exactly(self) -> None:
        result = run_scanner_strata(
            SCANNER,
            compiler="c++",
            compile_timeout=30,
            scan_timeout=30,
            denominator_bound=20,
            strata=(("tiny", 100, 10),),
        )[0]
        self.assertEqual(result.primitive_population, 1249)
        self.assertEqual(result.prior_excluded, 8)
        self.assertEqual(result.evaluated_population, 1241)
        self.assertEqual(len(result.candidates), 10)
        self.assertEqual(result.calibration.parameter, CALIBRATION_PARAMETER)
        self.assertEqual(result.calibration.discovery_good_primes, 19)
        self.assertEqual(result.calibration.held_good_primes, 21)
        self.assertNotIn(CALIBRATION_PARAMETER, {item.parameter for item in result.candidates})


class Mestre0430313946ArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = json.loads(ARTIFACT.read_text())

    def test_frozen_scale_boundary_and_broad_population_are_exact(self) -> None:
        artifact = self.artifact
        boundary = artifact["frozen_scale_boundary"]
        scan = artifact["modular_scan"]
        self.assertEqual(sha256(SCALE_ARTIFACT), FROZEN_SCALE_ARTIFACT_SHA256)
        self.assertEqual(boundary["artifact_sha256"], FROZEN_SCALE_ARTIFACT_SHA256)
        self.assertEqual(boundary["result_sha256"], FROZEN_SCALE_RESULT_SHA256)
        self.assertEqual(boundary["exact_prior_parameters"], [str(i) for i in range(1, 9)])
        self.assertTrue(scan["bands_disjoint"])
        self.assertEqual(
            scan["global_box_exhausted"],
            {
                "numerator": [1, 30000],
                "denominator": [1, 1000],
                "primitive_positive_rationals": 18_244_819,
                "exact_prior_excluded": 8,
                "evaluated": 18_244_811,
            },
        )
        strata = {record["name"]: record for record in scan["strata"]}
        self.assertEqual(strata["global"]["stdout_sha256"], "ce1066d2ad457b2721179e3e1df9eca8090a4cc5baef124dbb6db066fc4f6f5a")
        self.assertEqual(strata["medium"]["stdout_sha256"], "171afb5cbab50b296bd18f6edd5ab1629de38cab86eaa79681e636e11daf0df2")
        self.assertEqual(strata["low"]["stdout_sha256"], "ed1fafd20674c14a13ce91ab72328a0bd839fcf43115a0a538d3c9fc83e57f0b")

    def test_T5_has_an_unconditional_rank_twelve_certificate(self) -> None:
        calibration = self.artifact["calibration_T5"]
        certificate = calibration["finite_reduction_certificate"]
        coefficients = tuple(Q(value) for value in calibration["short_weierstrass_coefficients"])
        points = tuple(
            (Q(record["x"]), Q(record["y"]))
            for record in certificate["saturation"]["saturated_basis"]
        )
        self.assertEqual(len(points), 12)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in points))
        signatures = find_mod2_reduction_certificate(
            coefficients, points, prime_bound=certificate["certificate_prime_bound"]
        )
        self.assertEqual(
            [signature.prime for signature in signatures],
            certificate["certificate_primes"],
        )
        self.assertEqual(combined_mod2_rank(signatures, len(points)), 12)
        self.assertEqual(
            find_two_torsion_certificate_prime(coefficients, prime_bound=200), 37
        )
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 12)
        self.assertTrue(
            calibration["exact_log_conductor_bound"]["strict_target_proved_exactly"]
        )
        self.assertEqual(
            calibration["conductor_replay"],
            {
                "minimal_model": [0, 0, 0, -236126386510083, 1672228464226272736082],
                "conductor": "42267144357590204581790184193541040",
                "log_conductor": "79.7293181239103885582859998012530932967366034437835865568055",
                "minimal_discriminant": "-365437713046325163441114479198250077699558400",
                "root_number": 1,
            },
        )

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is required")
    def test_T5_exact_conductor_replay_is_stable(self) -> None:
        expected = self.artifact["calibration_T5"]["conductor_replay"]
        replay = capped_minimal_curve_data(
            family_coefficients(CALIBRATION_PARAMETER),
            timeout=8,
            stack_bytes=256_000_000,
        )
        self.assertEqual(replay, expected)

    def test_conductor_first_and_staged_frontier_are_pinned(self) -> None:
        artifact = self.artifact
        conductor = artifact["conductor_first_screen"]
        protocol = artifact["point_search_protocol"]
        self.assertTrue(conductor["population_closed_before_point_or_rank_triage"])
        self.assertEqual(conductor["selected_population"], 128)
        self.assertEqual(conductor["completed"], 118)
        self.assertEqual(conductor["timeouts"], 10)
        self.assertEqual(conductor["errors"], 0)
        self.assertEqual(conductor["singular_rejections"], 0)
        self.assertEqual(protocol["same_height_retries"], 0)
        self.assertEqual(protocol["completed_stage_calls"], 32)
        self.assertEqual(protocol["maximum_stable_numerical_rank"], 13)
        self.assertEqual(artifact["target"]["hits"], [])
        lead = next(
            record
            for record in conductor["records"]
            if record["parameter"] == "151/40"
        )
        self.assertEqual(
            lead["conductor_phase"]["log_conductor"],
            "131.095611084338516107896471786977870032999001111887998478298",
        )
        self.assertEqual(lead["conductor_phase"]["root_number"], -1)
        self.assertEqual(lead["point_stages"]["H250000"]["stable_numerical_rank"], 13)
        self.assertEqual(lead["point_stages"]["H1000000"]["stable_numerical_rank"], 13)
        self.assertEqual(
            lead["point_stages"]["H1000000"]["finite_reduction_attempt"],
            {"status": "not triggered", "trigger_stable_numerical_rank": 21},
        )

    def test_provenance_and_result_are_pinned(self) -> None:
        artifact = self.artifact
        provenance = artifact["provenance"]
        self.assertEqual(provenance["script_sha256"], sha256(SCRIPT))
        self.assertEqual(provenance["scanner_sha256"], sha256(SCANNER))
        self.assertEqual(sha256(SCRIPT), "00722f6cce4aeb99d87d20970933abccf29310afea44ecd1febcce5ba3f9f099")
        self.assertEqual(sha256(SCANNER), "010a3c50b058d9f2256a50e68347c86fe8e59a52e088022d7205334459fe883a")
        self.assertEqual(
            artifact["result_sha256"],
            "c73afd3abc9af90bedb98b4faca38375bcb81f638fcb3bfa42191b8782e2c87e",
        )
        self.assertEqual(sha256(ARTIFACT), "546cfc676b28f6956808b2698260d3bab4f9490dab5f2efc195f487ab6a2e514")


if __name__ == "__main__":
    unittest.main()

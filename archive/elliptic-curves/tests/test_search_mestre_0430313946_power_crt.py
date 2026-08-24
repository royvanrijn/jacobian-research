from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


ELLIPTIC_DIRECTORY = Path(__file__).resolve().parents[1]
CAS_DIRECTORY = ELLIPTIC_DIRECTORY / "cas"
REPO_ROOT = ELLIPTIC_DIRECTORY.parent
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from search_mestre_0430313946_power_crt import (  # noqa: E402
    ROOTS,
    SELECTED_POWER_PRIMES,
    combine_selected_balls,
    derive_exact_geometry,
    discover_local_profiles,
    enumerate_lattice_population,
)


SCRIPT = CAS_DIRECTORY / "search_mestre_0430313946_power_crt.py"
ARTIFACT = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_mestre_0430313946_power_crt.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Mestre0430313946PowerCRTGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
        cls.core, cls.geometry = derive_exact_geometry(construction)
        cls.profiles, cls.groups = discover_local_profiles(cls.core)

    def test_exact_discriminant_factor_and_fixed_divisor(self) -> None:
        geometry = self.geometry
        self.assertEqual(geometry["primitive_quartic_discriminant_degree"], 20)
        self.assertEqual(
            geometry["primitive_quartic_discriminant_coefficient_content"], 529_200
        )
        self.assertEqual(
            geometry["factorization"]["text"],
            "(2*T-21)^2*(2*T+21)^2*Q16(T)",
        )
        self.assertEqual(
            geometry["primitive_pair_homogeneous_fixed_divisor"]["value"], 28
        )
        self.assertEqual(
            geometry["short_jacobian"]["weierstrass_discriminant_core_multiplier"],
            4_499_817_235_200,
        )

    def test_complete_p4_discovery_and_selected_clean_unions_are_pinned(self) -> None:
        self.assertEqual(len(self.profiles), 46)
        expected = {
            5: [(2, 5, 4, "split multiplicative"), (3, 5, 4, "split multiplicative")],
            11: [
                (3, 11, 4, "split multiplicative"),
                (8, 11, 4, "split multiplicative"),
                (39, 121, 4, "split multiplicative"),
                (50, 121, 5, "split multiplicative"),
                (71, 121, 5, "split multiplicative"),
                (82, 121, 4, "split multiplicative"),
            ],
            13: [
                (0, 13, 4, "split multiplicative"),
                (4, 13, 5, "split multiplicative"),
                (9, 13, 5, "split multiplicative"),
            ],
            37: [
                (674, 1369, 4, "nonsplit multiplicative"),
                (695, 1369, 4, "nonsplit multiplicative"),
            ],
            43: [
                (914, 1849, 4, "split multiplicative"),
                (935, 1849, 4, "split multiplicative"),
            ],
        }
        self.assertEqual(tuple(sorted(self.groups)), SELECTED_POWER_PRIMES)
        for prime, records in expected.items():
            actual = [
                (
                    item["residue"],
                    item["modulus"],
                    item["forced_core_valuation"],
                    item["reduction"]["kind"],
                )
                for item in self.groups[prime]
            ]
            self.assertEqual(actual, records)

    def test_exact_gauss_box_is_outside_the_active_scan_box(self) -> None:
        classes = combine_selected_balls(self.groups)
        population, metadata = enumerate_lattice_population(
            classes, self.groups, coefficient_radius=24
        )
        self.assertEqual(len(classes), 144)
        self.assertEqual(len(population), 34_885)
        self.assertEqual(metadata["unique_sign_quotiented_representatives"], 34_889)
        self.assertEqual(metadata["active_global_box"]["exact_overlap_count"], 4)
        self.assertEqual(
            metadata["active_global_box"]["exact_overlap_parameters"],
            ["21/2", "46323/2456", "73801/2098", "91309/2973"],
        )
        self.assertTrue(
            all(
                numerator > 100_000 or denominator > 5_000
                for numerator, denominator in population
            )
        )
        self.assertEqual(
            metadata["outside_population_sha256"],
            "1fbe2e65788aa02546995c21d65d7fbd8f16a659dfc3d696acd7dcee6a41bb82",
        )


class Mestre0430313946PowerCRTArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = json.loads(ARTIFACT.read_text())

    def test_conductor_population_is_closed_and_complete(self) -> None:
        artifact = self.artifact
        self.assertTrue(
            artifact["ranking"]["selection"][
                "population_closed_before_conductor_calls"
            ]
        )
        self.assertEqual(artifact["ranking"]["scored_nonsingular_population"], 34_885)
        self.assertEqual(artifact["outcome"]["conductor_completed"], 20)
        self.assertEqual(artifact["outcome"]["conductor_timeouts"], 0)
        self.assertEqual(artifact["outcome"]["conductor_errors"], 0)
        self.assertEqual(artifact["outcome"]["subtarget_conductors"], 1)
        self.assertEqual(artifact["outcome"]["subtarget_triaged"], 1)
        self.assertEqual(artifact["outcome"]["target_hits"], [])

    def test_new_subtarget_fiber_and_bounded_point_frontier_are_pinned(self) -> None:
        record = next(
            item
            for item in self.artifact["conductor_first_records"]
            if item["t"] == "209001/3868"
        )
        conductor = record["conductor_phase"]
        self.assertEqual(
            conductor["conductor"],
            "731876369999526151436443857289480698073586382504666247966616525641014396190770",
        )
        self.assertTrue(
            conductor["decimal_log_recomputed_from_exact_conductor"].startswith(
                "179.2894935808073638872297092"
            )
        )
        self.assertEqual(conductor["root_number"], 1)
        self.assertTrue(conductor["below_strict_log_conductor_target"])
        triage = record["point_rank_triage"]
        self.assertEqual(triage["visible_stable_numerical_rank"], 10)
        self.assertEqual(triage["maximum_stable_numerical_rank"], 10)
        self.assertEqual(
            [item["height_bound"] for item in triage["searches"]],
            [50_000, 250_000, 1_000_000],
        )
        self.assertEqual(
            [item["accumulated_distinct_abscissas"] for item in triage["searches"]],
            [12, 12, 12],
        )
        self.assertIsNone(triage["finite_reduction_certificate_attempt"])

    def test_T5_calibration_is_recovered_and_excluded(self) -> None:
        calibration = self.artifact["calibration"]
        self.assertEqual(calibration["parameter"], "5")
        self.assertTrue(calibration["prior_values_recovered"])
        self.assertTrue(
            calibration["excluded_from_lattice_and_conductor_candidate_population"]
        )
        self.assertEqual(calibration["visible_stable_numerical_rank"], 10)
        self.assertTrue(
            calibration["conductor"]["log_conductor"].startswith(
                "79.729318123910"
            )
        )

    def test_default_artifact_provenance_is_pinned(self) -> None:
        artifact = self.artifact
        self.assertEqual(artifact["provenance"]["script_sha256"], sha256(SCRIPT))
        self.assertEqual(
            sha256(SCRIPT),
            "65016b98307317acde67b6c53eb8f26aa5050e85185d187609ea80fbcc7dfd5f",
        )
        self.assertEqual(
            artifact["result_sha256"],
            "811883853a5e2246947bf173619ad64e8c45291b7624aa5e14f9a36d083b5222",
        )
        self.assertEqual(
            sha256(ARTIFACT),
            "3f0cc34e75cb67d2a2482d56115a2c91f5cce6097c0fd834875f9d673ffaf9b8",
        )


if __name__ == "__main__":
    unittest.main()

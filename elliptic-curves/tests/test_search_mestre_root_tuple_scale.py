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
from search_mestre_root_tuple_scale import (  # noqa: E402
    NAGAO_NORMALIZED_ROOTS,
    classify_nonreflection,
    compiled_enumeration,
    max14_calibration,
    point_on_short_curve,
    primitive_visible_points,
    quartic_point_to_jacobian,
    tuple_digest,
    verify_enumerator_records,
)


ARTIFACT = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_mestre_root_tuple_scale.json"
)
SCRIPT = CAS_DIRECTORY / "search_mestre_root_tuple_scale.py"
CPP_SOURCE = CAS_DIRECTORY / "enumerate_mestre_root_tuples_scale.cpp"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipUnless(shutil.which("c++"), "a C++17 compiler is required")
class MestreRootTupleScaleGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.enumeration = compiled_enumeration(50)
        cls.calibration_enumeration = compiled_enumeration(14)
        verify_enumerator_records(cls.enumeration)
        verify_enumerator_records(cls.calibration_enumeration)
        cls.nonsingular, cls.singular, cls.witnesses = classify_nonreflection(
            cls.enumeration
        )

    def test_max_root_fifty_enumeration_is_pinned(self) -> None:
        enumeration = self.enumeration
        self.assertEqual(enumeration.normalized_count, 1_032_506)
        self.assertEqual(enumeration.obstruction_count, 4_099)
        self.assertEqual(enumeration.reflection_count, 3_966)
        self.assertEqual(enumeration.nonreflection_count, 133)
        self.assertEqual(
            tuple_digest(enumeration.obstruction_roots),
            "4093a3d9a203f067e02bd2126f3d287c33aad4a79a190fb8d2e96014b75b37b4",
        )
        self.assertEqual(
            tuple_digest(enumeration.nonreflection_roots),
            "d142a01ada2aff0114e8ac44b4c41bef2d5660ea8a27319c9085d3a7f70ce995",
        )

    def test_exact_nonsingularity_filter_and_nagao_calibration(self) -> None:
        self.assertEqual(len(self.nonsingular), 44)
        self.assertEqual(len(self.singular), 89)
        self.assertEqual(set(self.witnesses.values()), {1})
        self.assertIn(NAGAO_NORMALIZED_ROOTS, self.nonsingular)
        self.assertEqual(
            tuple_digest(self.nonsingular),
            "a892b1824af10c5e2dd428478778e8cc8db69fa057cfde03dbc10f1ab273433e",
        )

    def test_old_max_root_fourteen_boundary_is_recovered(self) -> None:
        calibration = max14_calibration(
            self.enumeration, self.calibration_enumeration
        )
        self.assertEqual(calibration["affine_normalized_root_tuples"], 1023)
        self.assertEqual(calibration["obstruction_zero"], 68)
        self.assertEqual(calibration["generically_nonsingular"], 59)
        self.assertEqual(calibration["generically_nonsingular_nonreflection"], 2)
        self.assertEqual(
            calibration["nonsingular_nonreflection_tuples"],
            [[0, 1, 7, 8, 9, 11], [0, 2, 8, 9, 11, 14]],
        )

    def test_general_covariant_map_checks_twelve_points_exactly(self) -> None:
        construction = SixRootMestreConstruction(
            tuple(Q(root) for root in (0, 4, 30, 31, 39, 46))
        )
        parameter = Q(5)
        quartic_points = primitive_visible_points(construction, parameter)
        jacobian_points = tuple(
            quartic_point_to_jacobian(construction, parameter, point)
            for point in quartic_points
        )
        coefficients = construction.primitive_jacobian_coefficients(parameter)
        self.assertEqual(len(quartic_points), 12)
        self.assertEqual(len({point[0] for point in quartic_points}), 12)
        self.assertEqual(len(jacobian_points), 12)
        self.assertTrue(
            all(point_on_short_curve(coefficients, point) for point in jacobian_points)
        )


class MestreRootTupleScaleArtifactTests(unittest.TestCase):
    def test_default_artifact_records_the_closed_leakage_free_screen(self) -> None:
        artifact = json.loads(ARTIFACT.read_text())
        enumeration = artifact["enumeration"]
        population = artifact["specialization_screen"]["population"]
        screen = artifact["specialization_screen"]
        self.assertEqual(enumeration["nonreflection_generically_nonsingular_count"], 44)
        self.assertTrue(enumeration["nagao_calibration"]["recovered"])
        self.assertEqual(enumeration["families_beyond_prior_bound"], 42)
        self.assertEqual(
            enumeration["new_families_excluding_known_nagao_calibration"], 41
        )
        self.assertEqual(population["proposed_integer_fibers"], 336)
        self.assertEqual(population["admissible_fibers"], 202)
        self.assertEqual(population["conductor_completed"], 202)
        self.assertEqual(population["conductor_timeouts"], 0)
        self.assertEqual(population["conductor_errors"], 0)
        self.assertEqual(population["visible_triage_completed"], 202)
        self.assertEqual(population["maximum_visible_stable_numerical_rank"], 10)
        self.assertEqual(screen["maximum_augmented_stable_numerical_rank"], 12)
        self.assertEqual(screen["maximum_escalated_stable_numerical_rank"], 12)
        self.assertEqual(screen["target_hits"], [])
        self.assertEqual(artifact["target"]["hits"], [])
        self.assertTrue(
            screen["protocol"][
                "conductor_population_closed_before_point_or_rank_triage"
            ]
        )

    def test_default_artifact_provenance_is_pinned(self) -> None:
        artifact = json.loads(ARTIFACT.read_text())
        provenance = artifact["provenance"]
        self.assertEqual(provenance["script_sha256"], sha256(SCRIPT))
        self.assertEqual(provenance["compiled_source_sha256"], sha256(CPP_SOURCE))
        self.assertEqual(
            sha256(SCRIPT),
            "5e7228b95ae995019fbc50b9f7667de41e06a86b4490f0feacff5702bb5cc174",
        )
        self.assertEqual(
            sha256(CPP_SOURCE),
            "31650333800698201819eddc91bf228089824bca026c629c9360683324a69eb5",
        )
        self.assertEqual(
            artifact["result_sha256"],
            "7c3f451a92f208d241955d2500cdcf416d772e919bb54f7181f5c40fd8f53def",
        )
        self.assertEqual(
            sha256(ARTIFACT),
            "fd2dccb1fd08aad70857df7ca19df77bd521e2be017b98f5579a748fd26cfc14",
        )


if __name__ == "__main__":
    unittest.main()

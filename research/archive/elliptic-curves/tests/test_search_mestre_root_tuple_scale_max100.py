#!/usr/bin/env python3
"""Focused replay checks for the standalone max-root-100 Mestre census."""

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
TESTS = ROOT / "elliptic-curves" / "tests"
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from search_mestre_root_tuple_scale import (  # noqa: E402
    classify_nonreflection,
    point_on_short_curve,
    tuple_digest,
    verify_enumerator_records,
)
from search_mestre_root_tuple_scale_max100 import (  # noqa: E402
    EXPECTED_MAX100_COUNTS,
    EXPECTED_MAX100_NONREFLECTION_SHA256,
    EXPECTED_MAX100_NONSINGULAR_SHA256,
    EXPECTED_MAX100_OBSTRUCTION_SHA256,
    EXPECTED_MAX50_NONSINGULAR_SHA256,
    EXPECTED_NEW_FAMILY_SHA256,
    compiled_enumeration_max100,
    family_feature_record,
    predeclare_h5000_fibers,
    result_digest,
    select_family_tranche,
    stable_json_digest,
)


SCRIPT = CAS / "search_mestre_root_tuple_scale_max100.py"
ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale_max100.json"
FROZEN_CPP = CAS / "enumerate_mestre_root_tuples_scale.cpp"
FROZEN_DRIVER = CAS / "search_mestre_root_tuple_scale.py"
FROZEN_TEST = TESTS / "test_search_mestre_root_tuple_scale.py"
FROZEN_ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale.json"
EXPECTED_SCRIPT_SHA256 = (
    "34677c38be30aa15e99b3239a6d487a51c158fa33326826d37ceead310555600"
)
EXPECTED_ARTIFACT_SHA256 = (
    "63dcd39555ad8b39c7b584a16663164bf73e6c6c59906b6a230bfa9b9f65a3bb"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRootTupleScaleMax100Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.enumeration, _ = compiled_enumeration_max100(
            FROZEN_CPP, compile_timeout=30, enumeration_timeout=30
        )
        verify_enumerator_records(cls.enumeration)
        cls.nonsingular, cls.singular, cls.witnesses = classify_nonreflection(
            cls.enumeration
        )
        cls.old_families = tuple(
            roots for roots in cls.nonsingular if roots[-1] <= 50
        )
        cls.new_families = tuple(
            roots for roots in cls.nonsingular if roots[-1] > 50
        )
        cls.features = [
            family_feature_record(roots) for roots in cls.new_families
        ]
        cls.selected, cls.selection = select_family_tranche(cls.features)
        cls.h5000_ids, cls.h5000_selection = predeclare_h5000_fibers(
            cls.selected
        )

    def test_frozen_inputs_and_new_files_are_pinned(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            self.data["provenance"]["script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        frozen = self.data["frozen_max50_inputs"]
        self.assertTrue(frozen["all_frozen_files_read_only"])
        self.assertEqual(frozen["compiled_source_sha256"], sha256(FROZEN_CPP))
        self.assertEqual(frozen["driver_sha256"], sha256(FROZEN_DRIVER))
        self.assertEqual(frozen["test_sha256"], sha256(FROZEN_TEST))
        self.assertEqual(frozen["artifact_sha256"], sha256(FROZEN_ARTIFACT))
        self.assertEqual(
            frozen["compiled_source_sha256"],
            "31650333800698201819eddc91bf228089824bca026c629c9360683324a69eb5",
        )
        self.assertEqual(
            frozen["driver_sha256"],
            "5e7228b95ae995019fbc50b9f7667de41e06a86b4490f0feacff5702bb5cc174",
        )
        self.assertEqual(
            frozen["test_sha256"],
            "a3930892e7e574161c0713c6c9b7c9f5aee0aa74e8e7acb89c250e2f9975d7c3",
        )
        self.assertEqual(
            frozen["artifact_sha256"],
            "fd2dccb1fd08aad70857df7ca19df77bd521e2be017b98f5579a748fd26cfc14",
        )

    def test_complete_exact_max100_census_and_max50_exclusion(self) -> None:
        enumeration = self.enumeration
        self.assertEqual(
            (
                enumeration.normalized_count,
                enumeration.obstruction_count,
                enumeration.reflection_count,
                enumeration.nonreflection_count,
            ),
            EXPECTED_MAX100_COUNTS,
        )
        self.assertEqual(
            tuple_digest(enumeration.obstruction_roots),
            EXPECTED_MAX100_OBSTRUCTION_SHA256,
        )
        self.assertEqual(
            tuple_digest(enumeration.nonreflection_roots),
            EXPECTED_MAX100_NONREFLECTION_SHA256,
        )
        self.assertEqual(len(self.nonsingular), 235)
        self.assertEqual(len(self.singular), 542)
        self.assertEqual(
            tuple_digest(self.nonsingular), EXPECTED_MAX100_NONSINGULAR_SHA256
        )
        self.assertEqual(set(self.witnesses.values()), {1})
        self.assertEqual(len(self.old_families), 44)
        self.assertEqual(
            tuple_digest(self.old_families), EXPECTED_MAX50_NONSINGULAR_SHA256
        )
        self.assertEqual(len(self.new_families), 191)
        self.assertEqual(tuple_digest(self.new_families), EXPECTED_NEW_FAMILY_SHA256)
        census = self.data["census"]
        self.assertEqual(census["complete_diameter_prefix"], [5, 100])
        self.assertEqual(census["open_diameter_remainder"], [])
        self.assertEqual(
            census["affine_normalized_primitive_reflection_quotient_count"],
            36_475_792,
        )
        self.assertEqual(census["degree_five_obstruction_zero_count"], 33_945)
        self.assertEqual(census["reflection_obstruction_zero_count"], 33_168)
        self.assertEqual(census["nonreflection_obstruction_zero_count"], 777)
        self.assertEqual(census["nonreflection_generically_nonsingular_count"], 235)
        self.assertEqual(census["nonreflection_generically_singular_count"], 542)
        self.assertEqual(census["genuinely_new_diameter_51_to_100_family_count"], 191)
        self.assertTrue(census["exact_obstruction_replayed_in_python"])

    def test_rank_blind_family_and_fiber_selections_replay_exactly(self) -> None:
        selection = self.data["rank_blind_selection"]
        self.assertFalse(selection["selection_uses_conductor"])
        self.assertFalse(selection["selection_uses_point_search"])
        self.assertFalse(selection["selection_uses_numerical_or_algebraic_rank"])
        self.assertEqual(selection["selection_eligible_family_count"], 165)
        self.assertEqual(selection["local_coverage_ineligible_family_count"], 26)
        self.assertEqual(
            stable_json_digest(self.features),
            selection["full_new_family_feature_sha256"],
        )
        self.assertEqual(self.selected, selection["selected_families"])
        self.assertEqual(
            self.selection["selected_family_sha256"],
            "d6b9e04b19f119fc67ab9308384fd633e751641ab3470dca9298e4f7e683b2d0",
        )
        self.assertEqual(
            Counter(record["selection_stratum"] for record in self.selected),
            Counter(
                {
                    "top-20 fixed-panel local score": 20,
                    "geometry-51-60": 4,
                    "geometry-61-70": 4,
                    "geometry-71-80": 4,
                    "geometry-81-90": 4,
                    "geometry-91-100": 4,
                }
            ),
        )
        self.assertEqual(self.h5000_ids, selection["predeclared_h5000_identifiers"])
        self.assertEqual(
            self.h5000_selection["predeclared_fiber_sha256"],
            "368967bbe417851c76a6d69af722b059c2612fa98849a20a35b41987299cb8c0",
        )
        self.assertEqual(len(self.h5000_ids), 64)
        self.assertEqual(len(set(self.h5000_ids)), 64)

    def test_conductor_phase_is_closed_before_all_point_triage(self) -> None:
        screen = self.data["specialization_screen"]
        protocol = screen["protocol"]
        population = screen["population"]
        self.assertTrue(
            protocol["conductor_population_closed_before_any_point_or_rank_call"]
        )
        self.assertTrue(protocol["all_admissible_selected_family_fibers_receive_conductor"])
        self.assertTrue(protocol["h5000_population_predeclared_before_conductor_phase"])
        self.assertTrue(protocol["no_retries"])
        self.assertEqual(population["proposed_integer_fibers"], 320)
        self.assertEqual(population["inadmissible_fibers"], 58)
        self.assertEqual(population["admissible_conductor_records"], 262)
        self.assertEqual(population["conductor_completed"], 262)
        self.assertEqual(population["conductor_timeouts"], 0)
        self.assertEqual(population["conductor_errors"], 0)
        self.assertEqual(population["subtarget_conductors"], 262)
        conductor_records = screen["conductor_records"]
        self.assertEqual(len(conductor_records), 262)
        self.assertEqual(
            len({record["identifier"] for record in conductor_records}), 262
        )
        self.assertTrue(
            all(
                record["conductor_phase"]["status"]
                == "completed exact PARI minimal-model/conductor computation"
                for record in conductor_records
            )
        )
        self.assertTrue(
            all(
                record["conductor_phase"][
                    "below_strict_log_conductor_target_numerically"
                ]
                for record in conductor_records
            )
        )

    def test_all_h5000_records_are_exact_and_frontier_is_rank13(self) -> None:
        screen = self.data["specialization_screen"]
        population = screen["population"]
        records = screen["h5000_records"]
        self.assertEqual(len(records), 64)
        self.assertEqual(population["h5000_completed"], 64)
        self.assertEqual(population["h5000_timeouts"], 0)
        self.assertEqual(population["h5000_errors"], 0)
        self.assertEqual(population["maximum_stable_numerical_rank"], 13)
        self.assertEqual(
            population["stable_numerical_rank_histogram"],
            {"6": 3, "7": 7, "8": 6, "9": 6, "10": 23, "11": 17, "12": 1, "13": 1},
        )
        for record in records:
            triage = record["point_triage"]
            self.assertEqual(
                triage["status"],
                "completed exact H5000 point checks and numerical height triage",
            )
            roots = tuple(record["roots"])
            construction = SixRootMestreConstruction(
                tuple(Fraction(root) for root in roots)
            )
            coefficients = construction.primitive_jacobian_coefficients(
                Fraction(record["parameter"])
            )
            subset = tuple(
                (Fraction(point["x"]), Fraction(point["y"]))
                for point in triage["numerical_subset"]
            )
            self.assertEqual(len(subset), triage["stable_numerical_rank"])
            self.assertTrue(
                all(point_on_short_curve(coefficients, point) for point in subset)
            )
            self.assertEqual(
                record["finite_reduction_attempt"]["status"], "not triggered"
            )
        leader = max(
            records, key=lambda record: record["point_triage"]["stable_numerical_rank"]
        )
        self.assertEqual(leader["identifier"], "r0_6_47_55_70_80_t8")
        self.assertEqual(leader["point_triage"]["stable_numerical_rank"], 13)
        self.assertEqual(
            leader["conductor_phase"]["log_conductor"],
            "82.3515440580102287639528785405710145305597823305329922776592",
        )
        self.assertEqual(leader["conductor_phase"]["root_number"], -1)
        self.assertEqual(leader["point_triage"]["pool_point_count_modulo_inverse"], 61)
        self.assertEqual(self.data["target"]["hits"], [])
        self.assertEqual(
            self.data["result_sha256"], "fac27d39205588488cb9cdefb7a8eb6929ba9f1bed6a5df8d0561039ae273107"
        )
        self.assertEqual(result_digest(self.data), self.data["result_sha256"])
        self.assertEqual(self.data["provenance"]["owned_processes_remaining"], 0)


if __name__ == "__main__":
    unittest.main()

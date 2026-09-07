#!/usr/bin/env python3
"""Focused checks for the standalone max-root-200 Mestre continuation."""

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
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from search_mestre_root_tuple_scale import (  # noqa: E402
    point_on_short_curve,
    tuple_digest,
)
from search_mestre_root_tuple_scale_max100 import (  # noqa: E402
    EXPECTED_MAX100_COUNTS,
    EXPECTED_MAX100_NONREFLECTION_SHA256,
    EXPECTED_MAX100_OBSTRUCTION_SHA256,
    stable_json_digest,
)
from search_mestre_root_tuple_scale_max200 import (  # noqa: E402
    compiled_enumeration_max200,
    gf_l_rank_and_pivots,
    independent_normalized_count,
    mod3_independence_certificate,
    screen_result_digest,
    verify_enumerator_records_fast,
)


SOURCE = CAS / "enumerate_mestre_root_tuples_scale_max200.cpp"
SCRIPT = CAS / "search_mestre_root_tuple_scale_max200.py"
CENSUS = GENERATED / "elliptic_mestre_root_tuple_scale_max200_census.json"
ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale_max200.json"
EXPECTED_SOURCE_SHA256 = (
    "56f5111765315fefea45628066a2971894bb963469d48823f08b717ba91c0c3a"
)
EXPECTED_SCRIPT_SHA256 = (
    "405a2b9f7653c89af0e3e6caf2e77765cb4bfc88fccf88edffa67d3435aebf24"
)
EXPECTED_CENSUS_SHA256 = (
    "7270769007f9c130fce8b1813164373de9c6a5eb1c6d86cfe71b8c96fada161b"
)
EXPECTED_ARTIFACT_SHA256 = (
    "5e1b53e187520735efba46fc8fd9cbdd4dfd4284545a815f6416baf3be84f342"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRootTupleScaleMax200Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.census = json.loads(CENSUS.read_text())
        cls.data = json.loads(ARTIFACT.read_text())

    def test_files_and_frozen_checkpoint_are_pinned(self) -> None:
        self.assertEqual(sha256(SOURCE), EXPECTED_SOURCE_SHA256)
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(CENSUS), EXPECTED_CENSUS_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            self.census["provenance"]["compiled_source_sha256"],
            EXPECTED_SOURCE_SHA256,
        )
        self.assertEqual(
            self.census["provenance"]["script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        self.assertEqual(
            self.data["provenance"]["script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        self.assertEqual(
            self.data["input"]["census_checkpoint_sha256"],
            EXPECTED_CENSUS_SHA256,
        )
        self.assertEqual(
            self.census["frozen_max100_inputs"]["artifact_sha256"],
            "63dcd39555ad8b39c7b584a16663164bf73e6c6c59906b6a230bfa9b9f65a3bb",
        )

    def test_v2_enumerator_exactly_recovers_frozen_max100_prefix(self) -> None:
        enumeration, _, _ = compiled_enumeration_max200(
            SOURCE,
            100,
            compile_timeout=30,
            enumeration_timeout=30,
        )
        verify_enumerator_records_fast(enumeration)
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
        self.assertEqual(
            independent_normalized_count(100)["primitive_reflection_orbit_count"],
            enumeration.normalized_count,
        )
        prefix = self.census["exact_max100_prefix_recovery"]
        self.assertTrue(prefix["record_for_record_equal_to_frozen_enumerator"])

    def test_complete_exact_max200_census_and_generic_classification(self) -> None:
        census = self.census["census"]
        self.assertEqual(
            census["affine_normalized_primitive_reflection_quotient_count"],
            1_225_592_478,
        )
        self.assertEqual(census["degree_five_obstruction_zero_count"], 275_387)
        self.assertEqual(census["reflection_obstruction_zero_count"], 271_600)
        self.assertEqual(census["nonreflection_obstruction_zero_count"], 3_787)
        self.assertEqual(census["nonreflection_generically_nonsingular_count"], 1_038)
        self.assertEqual(census["nonreflection_generically_singular_count"], 2_749)
        self.assertEqual(
            census["genuinely_new_diameter_101_to_200_family_count"], 803
        )
        self.assertEqual(
            census["obstruction_tuple_sha256"],
            "745332fd1280f0775820892566e886ae6f6a566781375dedad45bf359ab10c9e",
        )
        self.assertEqual(
            census["nonsingular_nonreflection_tuple_sha256"],
            "fce8f6dd09f1a4fe91c0adbd392d77dfb670627d5a41ba02da593a2e018d8d25",
        )
        new_roots = tuple(
            tuple(roots)
            for roots in self.census["tuple_populations"][
                "genuinely_new_nonsingular_roots"
            ]
        )
        self.assertEqual(len(new_roots), 803)
        self.assertTrue(all(101 <= roots[-1] <= 200 for roots in new_roots))
        self.assertEqual(
            tuple_digest(new_roots),
            "b75e13d200ec8d20043b1387aafb0710e53570901fac9066e4b9b11465cc7f45",
        )
        self.assertEqual(census["all_nonsingularity_witness_parameters"], [1])
        self.assertTrue(
            census[
                "all_obstruction_records_replayed_by_independent_python_integer_formula"
            ]
        )
        self.assertEqual(
            census["independent_burnside_mobius_count"],
            {
                "primitive_unquotiented_count": 2_450_913_356,
                "primitive_reflection_fixed_count": 271_600,
                "primitive_reflection_orbit_count": 1_225_592_478,
            },
        )

    def test_every_admissible_panel_fiber_has_an_internal_exact_rank_certificate(self) -> None:
        panel = self.data["complete_panel_screen"]
        population = panel["population"]
        self.assertEqual(population["new_family_count"], 803)
        self.assertEqual(population["proposed_panel_fiber_count"], 6_424)
        self.assertEqual(population["admissible_panel_fiber_count"], 5_708)
        self.assertEqual(population["inadmissible_panel_fiber_count"], 716)
        self.assertEqual(population["exact_mod3_certificate_count"], 5_708)
        self.assertEqual(population["maximum_visible_certified_rank_lower_bound"], 11)
        self.assertEqual(
            population["visible_rank_lower_bound_histogram"],
            {"5": 10, "6": 168, "7": 592, "8": 124, "9": 829, "10": 3064, "11": 921},
        )
        families = panel["family_records"]
        self.assertEqual(len(families), 803)
        self.assertEqual(stable_json_digest(families), panel["family_records_sha256"])
        admissible_count = 0
        for family in families:
            self.assertEqual(len(family["fibers"]), 8)
            for fiber in family["fibers"]:
                if not fiber["admissible"]:
                    continue
                admissible_count += 1
                certificate = fiber["mod3_finite_reduction_certificate"]
                rows = [
                    row
                    for signature in certificate["signatures"]
                    for row in signature["rows"]
                ]
                rank, pivots = gf_l_rank_and_pivots(
                    rows, certificate["point_count"], 3
                )
                self.assertEqual(rank, certificate["combined_exact_rank_over_F3"])
                self.assertEqual(
                    [pivot + 1 for pivot in pivots],
                    certificate["independent_subset_indices_one_based"],
                )
                self.assertNotEqual(
                    certificate["rational_3_torsion_exclusion"]["group_order"] % 3,
                    0,
                )
        self.assertEqual(admissible_count, 5_708)

    def test_rank_aware_diversity_selection_and_followup_are_closed(self) -> None:
        selection = self.data["rank_aware_diversity_selection"]
        self.assertTrue(selection["population_closed_before_conductor_calls"])
        self.assertEqual(selection["selected_family_count"], 64)
        self.assertEqual(selection["global_rank_aware_family_keep"], 34)
        self.assertEqual(
            selection["diversity_decile_counts"],
            {f"{lower}-{lower + 9}": 3 for lower in range(101, 201, 10)},
        )
        self.assertEqual(
            selection["selected_identifier_sha256"],
            "707a1bd8ca3a8888ced34384a504964c937ad4d86c26eeeffc7009e5c740ac46",
        )
        followup = self.data["leader_followup"]
        protocol = followup["protocol"]
        population = followup["population"]
        self.assertTrue(protocol["conductor_population_closed_before_any_point_or_height_call"])
        self.assertEqual(population["selected_leaders"], 64)
        self.assertEqual(population["conductor_completed"], 64)
        self.assertEqual(population["conductor_timeouts"], 0)
        self.assertEqual(population["conductor_errors"], 0)
        self.assertEqual(population["subtarget_conductors"], 64)
        self.assertEqual(population["point_search_completed"], 64)
        self.assertEqual(population["point_search_timeouts"], 0)
        self.assertEqual(population["point_search_errors"], 0)
        self.assertEqual(population["maximum_stable_numerical_rank"], 14)
        self.assertEqual(
            population["stable_numerical_rank_histogram"],
            {"11": 23, "12": 25, "13": 13, "14": 3},
        )
        self.assertEqual(population["immediate_exact_gain_attempts"], 3)

    def test_three_rank14_gains_replay_as_exact_certificates(self) -> None:
        records = [
            record
            for record in self.data["leader_followup"]["records"]
            if record["point_triage"]["stable_numerical_rank"] == 14
        ]
        self.assertEqual(
            [record["identifier"] for record in records],
            [
                "r0_17_142_145_162_200_t7",
                "r0_25_57_104_116_148_t1",
                "r0_7_121_128_183_194_t1",
            ],
        )
        expected_logs = {
            "r0_17_142_145_162_200_t7": (
                "105.858106121602461999446914582229034611901464156493645384544"
            ),
            "r0_25_57_104_116_148_t1": (
                "98.0874192183937578879903215252112681448598558379041886577325"
            ),
            "r0_7_121_128_183_194_t1": (
                "104.040597421688051105809637910994658460739729372335191396216"
            ),
        }
        for record in records:
            construction = SixRootMestreConstruction(
                tuple(Fraction(root) for root in record["roots"])
            )
            coefficients = construction.primitive_jacobian_coefficients(
                Fraction(record["parameter"])
            )
            subset = tuple(
                (Fraction(point["x"]), Fraction(point["y"]))
                for point in record["point_triage"]["numerical_subset"]
            )
            self.assertEqual(len(subset), 14)
            self.assertTrue(
                all(point_on_short_curve(coefficients, point) for point in subset)
            )
            replay = mod3_independence_certificate(
                coefficients, subset, prime_bound=499
            )
            self.assertEqual(
                json.loads(json.dumps(replay)),
                record["immediate_exact_gain_attempt"]["mod3"],
            )
            self.assertEqual(replay["certified_algebraic_rank_lower_bound"], 14)
            self.assertEqual(
                record["immediate_exact_gain_attempt"][
                    "best_certified_algebraic_rank_lower_bound"
                ],
                14,
            )
            self.assertEqual(
                record["conductor_phase"]["log_conductor"],
                expected_logs[record["identifier"]],
            )
            self.assertEqual(record["conductor_phase"]["root_number"], 1)
            self.assertTrue(
                record["conductor_phase"][
                    "below_strict_log_conductor_target_numerically"
                ]
            )
        self.assertEqual(self.data["target"]["hits"], [])
        self.assertEqual(
            self.data["result_sha256"],
            "546907f506326b4878c070ff33f0bea46716d2443cd758797a821aef71b09139",
        )
        self.assertEqual(screen_result_digest(self.data), self.data["result_sha256"])
        self.assertTrue(self.data["provenance"]["worker_pool_joined_before_write"])
        self.assertEqual(self.data["provenance"]["owned_processes_remaining"], 0)


if __name__ == "__main__":
    unittest.main()

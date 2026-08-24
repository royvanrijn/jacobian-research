from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_section7_auxiliary_group_orbits.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_section7_auxiliary_group_orbits.json"
)
STREAM_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_section7_auxiliary_group_orbit_stream.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_section7_auxiliary_group_orbits", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class Section7AuxiliaryGroupOrbitTests(unittest.TestCase):
    def test_declared_slices_and_saturated_bases_are_exact(self) -> None:
        artifact, priority = MODULE.load_inputs(ROOT)
        self.assertEqual(len(MODULE.ORBIT_SPECIFICATIONS), 2)
        self.assertEqual(
            {(item.slope, item.intercept) for item in MODULE.ORBIT_SPECIFICATIONS},
            {(Q(1), Q(-4471, 339)), (Q(-1), Q(154687, 447))},
        )
        for specification in MODULE.ORBIT_SPECIFICATIONS:
            item = priority[specification.priority_index_zero_based]
            self.assertEqual(item.identifier, specification.label)
            self.assertEqual(len(specification.expected_saturated_basis), 7)
            self.assertEqual(len(specification.expected_height_subset_indices), 7)
            known = MODULE.h200000_parameters(
                artifact, specification.priority_index_zero_based
            )
            self.assertGreaterEqual(len(known), 9)
            program = MODULE.build_gp_orbit_program(specification, item, known)
            self.assertIn("ellfromeqn", program)
            self.assertIn("ellsaturation(E,K,20)", program)
            self.assertIn("P0,ellmul(E,Q0,2)", program)

    def test_exact_parser_and_generic_decontamination(self) -> None:
        parsed = MODULE.parse_gp_exact("[[1/2,[-1,0,1]],[-3,[0,0,0]]]")
        self.assertEqual(parsed[0][0], Q(1, 2))
        self.assertEqual(parsed[1][0], Q(-3))
        parameter = Q(5081, 47)
        visible_x = Q(MODULE.SECTION7_ROOTS[0]) + parameter
        self.assertIn(
            "visible-00-plus",
            MODULE.generic_abscissa_labels(parameter, visible_x),
        )
        companion = MODULE.SECTION7_LINEAR_COMPANION_SECTIONS[0]
        companion_x = companion.slope * parameter + companion.intercept
        self.assertIn(
            f"linear-{companion.label}",
            MODULE.generic_abscissa_labels(parameter, companion_x),
        )

    def test_artifact_replays_bounded_population_and_negative_gate(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"], MODULE.file_sha256(SCRIPT)
        )
        self.assertEqual(data["status"], "bounded_exact_auxiliary_group_orbit_complete")
        self.assertEqual(data["generation"]["raw_pullbacks"], 2 * 3**7)
        self.assertEqual(
            data["generation"]["unique_parameters_after_all_exact_exclusions"],
            4251,
        )
        self.assertEqual(data["generation"]["proxy_below_190_count"], 0)
        self.assertEqual(data["exact_conductors"]["attempted"], 32)
        self.assertEqual(
            data["exact_conductors"]["completed"]
            + len(data["exact_conductors"]["failures"]),
            32,
        )
        self.assertEqual(data["exact_conductors"]["subtarget_count"], 0)
        for record in data["orbit_records"]:
            self.assertEqual(record["stable_numerical_auxiliary_rank"], 7)
            self.assertEqual(record["saturation_index_improvement"], 64)
            self.assertEqual(record["coefficient_vector_count"], 3**7)
            self.assertEqual(record["raw_pullback_count"], 3**7)

    def test_complete_exclusion_stream_matches_digest(self) -> None:
        data = json.loads(STREAM_ARTIFACT.read_text(encoding="utf-8"))
        parameters = tuple(Q(value) for value in data["parameters"])
        self.assertEqual(len(parameters), 4251)
        self.assertEqual(
            MODULE.stream_sha256(parameters), data["parameter_stream_sha256"]
        )
        self.assertEqual(
            data["parameter_stream_sha256"],
            "b3209a0bf363b66367bbd1765b7e869f9e95260a5002f8911870313a81242b35",
        )
        self.assertTrue(
            all(MODULE.projective_height(value) > 200000 for value in parameters)
        )


if __name__ == "__main__":
    unittest.main()

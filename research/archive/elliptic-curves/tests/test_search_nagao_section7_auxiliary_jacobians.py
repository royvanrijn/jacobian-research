from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_section7_auxiliary_jacobians.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_section7_auxiliary_jacobians.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_section7_auxiliary_jacobians", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class Section7AuxiliaryJacobianTests(unittest.TestCase):
    def test_four_slices_are_exact_quartics_through_T0(self) -> None:
        self.assertEqual(len(MODULE.SLICE_SPECIFICATIONS), 4)
        for specification in MODULE.SLICE_SPECIFICATIONS:
            auxiliary = MODULE.make_auxiliary_slice(specification)
            self.assertEqual(len(auxiliary.quartic_coefficients), 5)
            self.assertEqual(
                auxiliary.quartic_value(MODULE.T0), specification.base_y**2
            )
            self.assertEqual(
                specification.slope * MODULE.T0 + specification.intercept,
                specification.base_x,
            )

    def test_birational_map_and_general_group_law_replay(self) -> None:
        auxiliary = MODULE.make_auxiliary_slice(MODULE.SLICE_SPECIFICATIONS[0])
        parameter, _ = MODULE.generic_intersection_parameters(auxiliary)[0]
        ordinate = MODULE.rational_square_root(auxiliary.quartic_value(parameter))
        assert ordinate is not None
        image = auxiliary.forward((parameter, ordinate))
        self.assertEqual(auxiliary.inverse(image), (parameter, ordinate))
        doubled = MODULE.weierstrass_multiply(
            auxiliary.weierstrass_coefficients, image, 2
        )
        self.assertEqual(
            doubled,
            MODULE.weierstrass_add(
                auxiliary.weierstrass_coefficients, image, image
            ),
        )
        inverse = auxiliary.inverse(doubled)
        self.assertIsNotNone(inverse)
        assert inverse is not None
        self.assertEqual(inverse[1] ** 2, auxiliary.quartic_value(inverse[0]))

    def test_generic_intersections_and_vector_bound_are_exact(self) -> None:
        for specification in MODULE.SLICE_SPECIFICATIONS:
            auxiliary = MODULE.make_auxiliary_slice(specification)
            intersections = MODULE.generic_intersection_parameters(auxiliary)
            self.assertEqual(len(intersections), 12)
        vectors = MODULE.coefficient_vectors(7)
        self.assertEqual(len(vectors), 6986)
        self.assertTrue(
            all(
                max(map(abs, vector)) <= MODULE.MAX_ABSOLUTE_COEFFICIENT
                and sum(map(abs, vector)) <= MODULE.MAX_L1_NORM
                for vector in vectors
            )
        )

    def test_generated_artifact_respects_nonduplication_gate(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["script_sha256"], MODULE.file_sha256(SCRIPT))
        self.assertTrue(data["method"]["no_naive_hyperellratpoints_calls"])
        self.assertEqual(
            data["method"]["naive_projective_height_exclusion"], 200000
        )
        self.assertEqual(
            data["generation"]["group_elements_attempted"],
            data["generation"]["vector_count_per_slice"]
            * data["generation"]["slice_count"],
        )
        for record in data["proxy_selection"]["top_proxy_records"]:
            self.assertGreater(record["projective_height"], 200000)


if __name__ == "__main__":
    unittest.main()

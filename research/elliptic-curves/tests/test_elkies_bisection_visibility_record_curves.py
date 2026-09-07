from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest

from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_irreducible_p


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "elliptic-curves/scripts"
    / "analyze_elkies_bisection_visibility_and_record_curves.py"
)
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_visibility_record_curves_v1.json"
)

SPEC = importlib.util.spec_from_file_location("bisection_visibility", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class BisectionVisibilityRecordCurvesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = json.loads(ARTIFACT.read_text())

    def test_pinned_artifact_is_current(self) -> None:
        self.assertEqual(
            MODULE.canonical_json(MODULE.build_artifact()),
            ARTIFACT.read_text(),
        )

    def test_rank28_complement_and_curve394_visibility(self) -> None:
        visibility = {
            item["parameter"]: item for item in self.artifact["visibility"]
        }
        rank28 = visibility["-9529/5471"]
        self.assertEqual(rank28["visible_span_dimension"], 1)
        self.assertEqual(rank28["invisible_quotient_dimension"], 10)
        self.assertEqual(
            rank28["canonical_complement_columns_zero_based"],
            [0, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        )
        self.assertEqual(visibility["3/8"]["invisible_quotient_dimension"], 0)

    def test_record_curve_irreducibility_witnesses(self) -> None:
        expected = {
            "rank29": 461,
            "curve273": 367,
            "curve302": 397,
            "curve398": 1009,
            "curve399": 83,
            "curve400": 157,
        }
        for label, prime in expected.items():
            record = self.artifact["record_curve_provenance"][label]
            witness = record["irreducible_mod_prime_witness"]
            self.assertEqual(record["primitive_polynomial_degree"], 24)
            self.assertEqual(witness["prime"], prime)
            self.assertEqual(len(witness["coefficients_high_to_low"]), 25)
            self.assertTrue(
                gf_irreducible_p(
                    witness["coefficients_high_to_low"], prime, ZZ
                )
            )
            self.assertEqual(record["rational_affine_parameters"], [])
            self.assertFalse(record["parameter_at_infinity"])

    def test_rank28_positive_control(self) -> None:
        control = self.artifact["rank28_positive_control"]
        self.assertEqual(control["factor_degrees"], [1, 23])
        self.assertEqual(control["linear_factor"], "5471*t+9529")
        self.assertEqual(control["recovered_parameter"], "-9529/5471")
        witness = control["irreducible_degree23_cofactor_mod_prime_witness"]
        self.assertEqual(witness["prime"], 197)
        self.assertTrue(
            gf_irreducible_p(
                witness["coefficients_high_to_low"], witness["prime"], ZZ
            )
        )


if __name__ == "__main__":
    unittest.main()

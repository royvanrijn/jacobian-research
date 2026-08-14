from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (  # noqa: E402
    SECTION7_CONSTRUCTION,
    SECTION7_CONSTRUCTOR_PARAMETER,
    SECTION7_JACOBIAN_RELATIONS,
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_PAPER_PARAMETER,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_JACOBIAN_RELATIONS,
    SECTION7_ROOTS,
    section7_primitive_quartic_coefficients,
)
from verify_nagao_section7_linear_sections import (  # noqa: E402
    classify_abscissae,
    classify_cubic_abscissae,
    classify_quintic_abscissae,
    classify_quadratic_abscissae,
    verify_k3_fiber_geometry,
)


Q = Fraction
SCRIPT = CAS / "verify_nagao_section7_linear_sections.py"
DATA = CAS / "nagao_1994_section7.py"
MOD_L = CAS / "mod_l_reduction_independence.py"
ARTIFACT = GENERATED / "elliptic_nagao_section7_linear_sections.json"


def add_short(
    coefficients: tuple[Fraction, ...],
    left: tuple[Fraction, Fraction] | None,
    right: tuple[Fraction, Fraction] | None,
) -> tuple[Fraction, Fraction] | None:
    if left is None:
        return right
    if right is None:
        return left
    _, _, _, coefficient_a, _ = coefficients
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1**2 + coefficient_a) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope**2 - x1 - x2
    return x3, -y1 + slope * (x1 - x3)


class NagaoSection7LinearSectionTests(unittest.TestCase):
    def test_parameters_and_explicit_quartic_normalization(self) -> None:
        self.assertEqual(SECTION7_ROOTS, (346, 260, 255, 146, 55, 0))
        self.assertEqual(SECTION7_PAPER_PARAMETER, Q(5081, 94))
        self.assertEqual(SECTION7_CONSTRUCTOR_PARAMETER, Q(5081, 47))
        for parameter in (Q(1), Q(5081, 47), Q(-23, 7)):
            self.assertEqual(
                section7_primitive_quartic_coefficients(parameter),
                primitive_quartic_coefficients(
                    SECTION7_CONSTRUCTION, parameter
                ),
            )
        geometry = verify_k3_fiber_geometry()
        self.assertEqual(geometry["surface_type"], "elliptic K3")
        self.assertEqual(geometry["infinity_fiber"]["type"], "I4")
        self.assertEqual(geometry["geometric_generic_rank_upper_bound"], 15)

    def test_nine_companion_formulas_and_relations_specialize_exactly(self) -> None:
        parameter = Q(37, 5)
        companions = {
            section.label: section.jacobian_point(parameter)
            for section in SECTION7_LINEAR_COMPANION_SECTIONS
        }
        quadratics = {
            section.label: section.jacobian_point(parameter)
            for section in SECTION7_QUADRATIC_COMPANION_SECTIONS
        }
        visible = tuple(
            quartic_point_to_short_jacobian(
                SECTION7_CONSTRUCTION, parameter, point
            )
            for point in primitive_visible_points(
                SECTION7_CONSTRUCTION, parameter
            )
        )
        basis = visible[:11] + (companions["plus-7/27"],)
        targets = {
            "visible-11": visible[11],
            **{
                label: point
                for label, point in companions.items()
                if label != "plus-7/27"
            },
            **quadratics,
        }
        coefficients = short_jacobian_coefficients(
            SECTION7_CONSTRUCTION, parameter
        )
        relations = {
            **SECTION7_JACOBIAN_RELATIONS,
            **SECTION7_QUADRATIC_JACOBIAN_RELATIONS,
        }
        for label, relation in relations.items():
            total = None
            for point, coefficient in zip(basis, relation):
                if coefficient:
                    summand = point if coefficient == 1 else (point[0], -point[1])
                    total = add_short(coefficients, total, summand)
            self.assertEqual(total, targets[label])

    def test_exact_linear_abscissa_classification(self) -> None:
        result = classify_abscissae()
        self.assertEqual(len(result["recovered_companion_sections"]), 6)
        self.assertIn("27*m - 43", result["nonsingular_slope_resultant_gcd"])
        expected = "n*(n - 346)*(n - 260)*(n - 255)*(n - 146)*(n - 55)"
        self.assertEqual(
            result["slope_plus_or_minus_one_intercept_eliminants"],
            {"-1": expected, "1": expected},
        )
        quadratic = classify_quadratic_abscissae()
        self.assertEqual(
            quadratic["genuine_quadratic_n_nonzero_groebner_basis"], ["1"]
        )
        self.assertEqual(len(quadratic["recovered_quadratic_sections"]), 3)
        self.assertEqual(
            quadratic["genuine_quadratic_n_zero_groebner_basis"],
            [
                "5373*k - 1389190",
                "(5373*m - 56)*(5373*m + 22)*(5373*m + 34)",
            ],
        )
        cubic = classify_cubic_abscissae()
        self.assertTrue(cubic["leading_cubic_coefficient_forced_zero"])
        self.assertEqual(cubic["new_cubic_sections"], 0)
        self.assertEqual(
            cubic["singular_reduced_groebner_basis"],
            [
                "G[1]=a",
                "G[2]=c",
                "G[3]=5373d-1389190",
                "G[4]=155113830117b3-12830724b-41888",
            ],
        )
        quintic = classify_quintic_abscissae()
        self.assertTrue(quintic["leading_quintic_coefficient_forced_zero"])
        self.assertEqual(quintic["new_quintic_sections"], 0)
        self.assertEqual(
            quintic["singular_reduced_groebner_basis"],
            [
                "G[1]=a5",
                "G[2]=a4",
                "G[3]=a3",
                "G[4]=a1",
                "G[5]=5373*a0-1389190",
                "G[6]=155113830117*a2^3-12830724*a2-41888",
            ],
        )

    def test_generated_theorem_artifact_and_hashes(self) -> None:
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["proved_consequences"]["generic_rank_at_least"], 12)
        self.assertEqual(
            data["proved_consequences"]["geometric_generic_rank_at_most"], 15
        )
        self.assertEqual(
            data["proved_consequences"][
                "polynomial_abscissa_degree_at_most_five_sections_classified"
            ],
            21,
        )
        self.assertTrue(
            data["quintic_abscissa_classification"][
                "leading_quintic_coefficient_forced_zero"
            ]
        )
        self.assertEqual(
            data["exact_specialization_independence"]["combined_column_rank"],
            12,
        )
        self.assertFalse(data["target_hit"])
        self.assertEqual(
            data["script_sha256"], hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(
            data["data_sha256"], hashlib.sha256(DATA.read_bytes()).hexdigest()
        )
        self.assertEqual(
            data["mod_l_engine_sha256"],
            hashlib.sha256(MOD_L.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import shutil
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from fermigier_mestre import (  # noqa: E402
    FermigierMestreFamily,
    ROOTS,
)
from mestre_root_tuples import (  # noqa: E402
    SixRootMestreConstruction,
    affine_normalized_integer_root_tuples,
    mestre_quartic_condition,
    normalize_integer_root_tuple,
)
from survey_mestre_root_tuples import p_adic_profile  # noqa: E402
from survey_mestre_root_tuples import (  # noqa: E402
    discriminant_degree_profile,
    nagao_base_change_profile,
)
from pari_bridge import minimal_curve_data  # noqa: E402


class MestreRootTupleTests(unittest.TestCase):
    def test_affine_integer_normalization(self) -> None:
        self.assertEqual(
            normalize_integer_root_tuple((30, 12, 10, 28, 26, 14)),
            (0, 1, 2, 8, 9, 10),
        )
        self.assertEqual(
            normalize_integer_root_tuple((11, 10, 4, 3, 2, 0)),
            (0, 1, 7, 8, 9, 11),
        )

    def test_bounded_normalized_enumerator_is_pinned(self) -> None:
        tuples = affine_normalized_integer_root_tuples(12)
        self.assertEqual(len(tuples), 412)
        quartic_count = sum(
            mestre_quartic_condition(
                SixRootMestreConstruction(tuple(Q(root) for root in roots)).polynomial
            )
            == 0
            for roots in tuples
        )
        self.assertEqual(quartic_count, 40)

    def test_max_root_fourteen_has_only_two_nonreflection_nonsingular_families(self) -> None:
        candidates = []
        for roots in affine_normalized_integer_root_tuples(14):
            construction = SixRootMestreConstruction(
                tuple(Q(root) for root in roots)
            )
            if (
                construction.is_quartic_family
                and not construction.is_reflection_symmetric
                and any(
                    construction.primitive_discriminant_polynomial
                )
            ):
                candidates.append(roots)
        self.assertEqual(
            candidates,
            [(0, 1, 7, 8, 9, 11), (0, 2, 8, 9, 11, 14)],
        )

    def test_degree_five_obstruction_rejects_a_generic_tuple(self) -> None:
        construction = SixRootMestreConstruction(
            tuple(Q(root) for root in (0, 1, 3, 7, 9, 12))
        )
        self.assertEqual(construction.quartic_condition, -14400)
        self.assertFalse(construction.is_quartic_family)
        with self.assertRaisesRegex(ValueError, "degree-five"):
            construction.quartic_coefficients(Q(1))

    def test_quartic_identity_and_twelve_points_are_exact(self) -> None:
        construction = SixRootMestreConstruction(
            tuple(Q(root) for root in (0, 1, 2, 8, 9, 10))
        )
        parameter = Q(2, 3)
        product = construction.product_coefficients(parameter)
        approximant = construction.square_approximant_coefficients(parameter)
        points = construction.visible_points(parameter)
        self.assertEqual(len(points), 12)
        self.assertEqual(len(set(points)), 12)
        for x, y in points:
            self.assertEqual(y**2, construction.quartic_value(parameter, x))
        for x in (Q(-2), Q(0), Q(7, 5)):
            product_value = sum(value * x**index for index, value in enumerate(product))
            approximant_value = sum(
                value * x**index for index, value in enumerate(approximant)
            )
            self.assertEqual(
                approximant_value**2 - product_value,
                parameter**2 * construction.quartic_value(parameter, x),
            )

    def test_general_api_recovers_pinned_fermigier_quartic(self) -> None:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
        self.assertEqual(construction.quartic_condition, 0)
        self.assertFalse(construction.is_reflection_symmetric)
        self.assertEqual(construction.quartic_content, 50616**2)
        self.assertEqual(construction.quartic_square_scale, 50616)
        for parameter in (Q(1), Q(7, 3)):
            general = construction.quartic_coefficients(parameter)
            pinned = tuple(
                reversed(FermigierMestreFamily.quartic_coefficients(parameter))
            )
            self.assertEqual(
                tuple(value / 50616**2 for value in general),
                pinned,
            )
            self.assertEqual(
                construction.primitive_quartic_discriminant(parameter),
                16 * FermigierMestreFamily.discriminant_factor(parameter),
            )
            self.assertEqual(
                construction.primitive_jacobian_coefficients(parameter),
                FermigierMestreFamily.coefficients(parameter),
            )

    def test_discriminant_interpolation_and_generic_singularity(self) -> None:
        fermigier = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
        self.assertEqual(len(fermigier.primitive_discriminant_polynomial) - 1, 20)
        self.assertEqual(
            fermigier.primitive_discriminant_value(Q(0)),
            16 * FermigierMestreFamily.discriminant_factor(Q(0)),
        )
        singular = SixRootMestreConstruction(
            tuple(Q(root) for root in (0, 1, 2, 4, 5, 6))
        )
        self.assertTrue(singular.is_quartic_family)
        self.assertTrue(
            all(value == 0 for value in singular.primitive_discriminant_polynomial)
        )

    def test_nagao_quadratic_base_change_degree_proxy(self) -> None:
        nagao = SixRootMestreConstruction(
            tuple(Q(root) for root in (-17, -16, 10, 11, 14, 17))
        )
        self.assertEqual(
            discriminant_degree_profile(nagao),
            {
                "status": "exact polynomial computation",
                "degree": 20,
                "squarefree_degree": 20,
                "repeated_factor_degree": 0,
            },
        )
        base_change = nagao_base_change_profile(nagao)
        self.assertEqual(
            base_change["squarefree_discriminant_degree_proxy"], 40
        )
        self.assertEqual(
            base_change["rational_collision_preimages"],
            [
                {
                    "collision_parameter_T": "3",
                    "rational_parameters_U": ["1287/239", "infinity"],
                    "collision_loss": 1,
                },
                {
                    "collision_parameter_T": "17",
                    "rational_parameters_U": ["-330/7", "13"],
                    "collision_loss": 1,
                },
            ],
        )

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is optional")
    def test_two_nonreflection_specialization_calibrations(self) -> None:
        cases = (
            ((0, 2, 8, 9, 11, 14), Q(19), "54.9708843220", -1, [9, 9]),
            ((0, 1, 7, 8, 9, 11), Q(8), "19.025832263", 1, [4, 4]),
        )
        for roots, parameter, log_prefix, root_number, bounds in cases:
            construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
            data = minimal_curve_data(
                construction.primitive_jacobian_coefficients(parameter),
                rank_effort=0,
                timeout=10,
                stack_bytes=64_000_000,
            )
            self.assertTrue(data["log_conductor"].startswith(log_prefix))
            self.assertEqual(data["root_number"], root_number)
            self.assertEqual(
                [
                    data["pari_ellrank"]["lower_bound"],
                    data["pari_ellrank"]["upper_bound"],
                ],
                bounds,
            )

    def test_visible_collision_profile_detects_section_coalescence(self) -> None:
        construction = SixRootMestreConstruction(
            tuple(Q(root) for root in (0, 1, 2, 8, 9, 10))
        )
        degeneracy = construction.visible_point_degeneracy(Q(1, 2))
        self.assertEqual(degeneracy.distinct_abscissae, 8)
        self.assertEqual(degeneracy.collision_loss, 4)
        self.assertEqual(degeneracy.zero_ordinates, 0)

    def test_fermigier_seven_adic_branching_is_replayed(self) -> None:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
        profile = p_adic_profile(construction, 7, 2, {})
        self.assertEqual(profile["fixed_discriminant_valuation"], 0)
        self.assertEqual(profile["root_counts_by_exponent"], [1, 7])
        self.assertEqual(
            profile["root_residues_at_requested_exponent"],
            [0, 7, 14, 21, 28, 35, 42],
        )
        self.assertEqual(profile["branching_ratio_last_to_first"], "7")


if __name__ == "__main__":
    unittest.main()

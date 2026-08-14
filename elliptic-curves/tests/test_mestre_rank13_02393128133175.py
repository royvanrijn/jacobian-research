#!/usr/bin/env python3
"""Focused exact tests for the new six-root rank-at-least-13 family."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))

import mestre_rank13_02393128133175 as formulas  # noqa: E402
from verify_mestre_rank13_02393128133175 import (  # noqa: E402
    BASE_CHANGE_CONSTANT,
    CENTERED_ROOTS,
    COMPANION_DATA,
    ROOTS,
    build_verification,
    companion_point,
    reconstruct_symbolically,
)


Q = Fraction


class MestreRank13NewFamilyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.symbolic = reconstruct_symbolically()
        cls.result = build_verification(full_polynomials=True)

    def test_root_geometry_and_symbolic_quartic_reconstruction(self) -> None:
        self.assertEqual(ROOTS, (0, 23, 93, 128, 133, 175))
        self.assertEqual(CENTERED_ROOTS, (-92, -69, 1, 36, 41, 83))
        self.assertEqual(formulas.CONSTRUCTION.quartic_condition, 0)
        self.assertEqual(formulas.CONSTRUCTION.quartic_content, Q(5_760_000))
        self.assertEqual(formulas.CONSTRUCTION.quartic_square_scale, Q(2_400))
        for parameter in (Q(1), Q(2), Q(-5, 3), Q(17, 4)):
            reconstructed = tuple(
                Q(str(value.subs(self.symbolic.t, parameter)))
                for value in self.symbolic.quartic_coefficients
            )
            self.assertEqual(
                reconstructed,
                formulas.primitive_quartic_coefficients(parameter),
            )

    def test_all_six_nonvisible_companion_squares(self) -> None:
        self.assertEqual(len(COMPANION_DATA), 6)
        self.assertEqual(
            tuple((row[0], row[1]) for row in COMPANION_DATA),
            (
                (31, Q(619, 3)),
                (-31, Q(619, 3)),
                (19, Q(304, 3)),
                (-19, Q(304, 3)),
                (1, Q(115, 3)),
                (-1, Q(115, 3)),
            ),
        )
        for parameter in (Q(1), Q(-2), Q(7, 3)):
            quartic = formulas.primitive_quartic_coefficients(parameter)
            for index in range(6):
                x_value, y_value = companion_point(parameter, index)
                value = sum(
                    coefficient * x_value**degree
                    for degree, coefficient in enumerate(quartic)
                )
                self.assertEqual(y_value**2, value)
        self.assertEqual(
            companion_point(Q(13, 5), 0),
            formulas.linear_extra_point(Q(13, 5)),
        )

    def test_split_infinity_is_pivotal_in_exact_mod3_certificate(self) -> None:
        specialization = self.result["specialization_certificate"]
        self.assertTrue(specialization["quartic_discriminant_nonzero"])
        self.assertEqual(specialization["displayed_point_count"], 14)
        self.assertEqual(
            specialization["visible_mod3"]["combined_exact_rank_over_F3"], 11
        )
        self.assertEqual(
            specialization["affine_with_one_companion_mod3"][
                "combined_exact_rank_over_F3"
            ],
            12,
        )
        all_companions = specialization["affine_with_all_six_companions_mod3"]
        self.assertEqual(all_companions["point_count"], 18)
        self.assertEqual(all_companions["combined_exact_rank_over_F3"], 12)
        self.assertEqual(
            all_companions["independent_subset_indices_one_based"],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
        )
        full = specialization["with_split_infinity_mod3"]
        self.assertEqual(full["combined_exact_rank_over_F3"], 13)
        self.assertEqual(
            full["independent_subset_indices_one_based"],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14],
        )
        self.assertEqual(
            full["certificate_primes"],
            [23, 29, 37, 41, 47, 59, 71, 97, 101, 127, 131, 137],
        )
        self.assertEqual(full["rational_3_torsion_exclusion"]["prime"], 7)
        self.assertEqual(full["rational_3_torsion_exclusion"]["group_order"], 10)
        self.assertEqual(
            full["point_sha256"],
            "0616cd332c65b76785cbb39f7c561533e1087b0dec6b1de3568101c1b1c38bfb",
        )
        self.assertEqual(
            full["independent_subset_sha256"],
            "f1c52d7f854ac75d81c5501d219a773250b65cb3660e4916823a66d90dfe262e",
        )

    def test_degree_40_squarefree_frontier(self) -> None:
        geometry = self.result["discriminant_geometry"]
        self.assertEqual(geometry["original_core_degree_in_T"], 20)
        self.assertTrue(geometry["original_core_irreducible_over_Q"])
        self.assertTrue(geometry["original_core_squarefree"])
        self.assertEqual(geometry["base_changed_core_degree_in_u"], 40)
        self.assertTrue(geometry["base_changed_core_irreducible_over_Q"])
        self.assertTrue(geometry["base_changed_core_squarefree"])
        self.assertEqual(
            geometry["original_core_coefficients_ascending_sha256"],
            "1f9ed250329133ae7510c14a6044fd396fecbbd9463426baf347f8e29108691f",
        )
        self.assertEqual(
            geometry["base_changed_core_coefficients_ascending_sha256"],
            "7accef47d6942c6d4cd531bdac5bd772fc9c768e6e61f322772f7ad9c0758751",
        )
        self.assertEqual(
            len(geometry["original_core_coefficients_ascending"]), 21
        )
        self.assertEqual(
            len(geometry["base_changed_core_coefficients_ascending"]), 41
        )
        self.assertEqual(BASE_CHANGE_CONSTANT, 2 * 3 * 7**4)


if __name__ == "__main__":
    unittest.main()

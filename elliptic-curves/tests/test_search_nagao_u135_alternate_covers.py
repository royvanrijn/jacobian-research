from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from alternate_quartic_covers import alternate_cover, mobius_preimage  # noqa: E402
from nagao_skew_height import load_rank17_target  # noqa: E402
from search_nagao_u135_alternate_covers import (  # noqa: E402
    best_cross_ratio_charts,
    enumerate_cover_plans,
    full_coset_identity_frontier,
    mask_indices,
    projective_height,
)


Q = Fraction


class AlternateCoverSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target = load_rank17_target(
            ROOT
            / "artifacts/generated-results/elliptic_nagao_rank17_frontier_certificate.json",
            Q(135, 2),
        )

    def test_projective_height_and_masks(self) -> None:
        self.assertEqual(projective_height(Q(-17, 23)), 23)
        self.assertEqual(mask_indices(0b10101, 5), (0, 2, 4))
        with self.assertRaises(ValueError):
            mask_indices(0, 5)

    def test_cross_ratio_chart_normalizes_three_known_parameters(self) -> None:
        coefficients = (Q(0), Q(0), Q(0), Q(-1), Q(1))
        basis = (
            (Q(0), Q(1)),
            (Q(-1), Q(1)),
            (Q(1), Q(1)),
            (Q(3), Q(5)),
            (Q(5), Q(11)),
        )
        cover = alternate_cover(coefficients, basis[0])
        chart = best_cross_ratio_charts(cover, basis, count=1)[0]
        parameters = {
            index: cover.curve_point_to_cover(point)[0]
            for index, point in enumerate(basis)
            if point != basis[0]
        }
        preimages = tuple(
            mobius_preimage(chart.matrix, parameters[index])
            for index in chart.basis_indices
        )
        self.assertEqual(preimages, (Q(0), Q(1), None))

    def test_small_full_coset_scan_is_exact_and_sorted(self) -> None:
        basis = self.target.saturated_basis[:3]
        frontier = full_coset_identity_frontier(
            self.target.jacobian_coefficients, basis, retain_count=7
        )
        self.assertEqual(len(frontier), 7)
        self.assertEqual([item[0] for item in frontier], sorted(item[0] for item in frontier))
        self.assertEqual(
            {item[0][2] for item in frontier}, set(range(1, 1 << len(basis)))
        )

    def test_low_weight_plan_count(self) -> None:
        basis = self.target.saturated_basis[:4]
        plans = enumerate_cover_plans(
            self.target.jacobian_coefficients,
            basis,
            maximum_subset_weight=1,
            charts_per_cover=1,
        )
        self.assertEqual(len(plans), 4)
        self.assertTrue(all(len(plan.charts) == 1 for plan in plans))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression tests for the structural-search groundwork."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from build_structural_search_groundwork import (  # noqa: E402
    build_manifest,
    canonical_digest,
)
from structural_search import (  # noqa: E402
    BranchDivisor,
    IntegralLattice,
    ProjectivePadicBall,
    ProjectiveRational,
    projective_congruence_lattice,
    ResidualSelmerBudget,
    TwistCharacter,
    V4CoverDecomposition,
    enumerate_isotropic_fibration_candidates,
    enumerate_k3_divisor_candidates,
    two_division_cubic,
)


Q = Fraction


class StructuralSearchGroundworkTests(unittest.TestCase):
    def test_two_division_cubic_for_short_and_generalized_models(self) -> None:
        self.assertEqual(two_division_cubic((0, 0, 0, -1, 1)), (1, -1, 0, 1))
        # y^2 + y = x^3 - x has b2=0, b4=-2, b6=1.
        self.assertEqual(two_division_cubic((0, 0, 1, -1, 0)), (1, -4, 0, 4))
        with self.assertRaises(ValueError):
            two_division_cubic((0, 0, 0, 1))

    def test_residual_selmer_budget_distinguishes_upper_bound_from_rank_gain(self) -> None:
        budget = ResidualSelmerBudget(
            ell=2,
            selmer_dimension=22,
            known_free_rank=20,
            rational_ell_torsion_dimension=0,
        )
        self.assertEqual(budget.rank_upper_bound, 22)
        self.assertEqual(budget.unexplained_selmer_dimension, 2)
        self.assertTrue(budget.target_not_excluded(21))
        self.assertIn("Sha", budget.to_record()["claim_boundary"])
        with self.assertRaises(ValueError):
            ResidualSelmerBudget(2, 19, 20)

    def test_k3_lattice_candidate_and_isotropic_enumeration(self) -> None:
        lattice = IntegralLattice(((0, 1), (1, 0)), ("F", "S"))
        candidates = enumerate_k3_divisor_candidates(
            lattice,
            fiber_vector=(1, 0),
            coefficient_bound=2,
            fiber_degrees=(1,),
            self_intersections=(-2, 0, 2, 4),
        )
        by_square = {candidate.self_intersection: candidate for candidate in candidates}
        self.assertEqual(by_square[-2].arithmetic_genus, 0)
        self.assertEqual(by_square[0].arithmetic_genus, 1)
        self.assertEqual(by_square[2].arithmetic_genus, 2)
        rays = enumerate_isotropic_fibration_candidates(
            lattice, ample_vector=(1, 1), coefficient_bound=2
        )
        self.assertIn((0, 1), rays)
        self.assertIn((1, 0), rays)

    def test_v4_genus_decomposition_and_twist_character(self) -> None:
        first = BranchDivisor.squarefree_polynomial(6, prefix="f")
        second = BranchDivisor.squarefree_polynomial(6, prefix="g")
        cover = V4CoverDecomposition(first, second)
        self.assertEqual(cover.quotient_genera, (2, 2, 5))
        self.assertEqual(cover.cover_genus, 9)
        twist = TwistCharacter("d", first)
        self.assertEqual(twist.base_change_genus, 2)
        self.assertEqual(
            twist.rank_lower_bound_after_base_change(
                base_rank_lower_bound=12, twist_rank_lower_bound=1
            ),
            13,
        )
        with self.assertRaises(ValueError):
            V4CoverDecomposition(first, first)

    def test_projective_padic_affine_and_infinity_charts(self) -> None:
        affine = ProjectivePadicBall(7, 2, "affine", residue=3)
        self.assertTrue(affine.matches(ProjectiveRational.normalized(3, 1)))
        self.assertTrue(affine.matches(ProjectiveRational.normalized(52, 1)))
        self.assertFalse(affine.matches(ProjectiveRational.normalized(1, 7)))

        infinity = ProjectivePadicBall(7, 2, "infinity", residue=0)
        self.assertTrue(infinity.matches(ProjectiveRational.normalized(1, 49)))
        self.assertFalse(infinity.matches(ProjectiveRational.normalized(7, 1)))

        lattice = projective_congruence_lattice(
            (
                ProjectivePadicBall(5, 1, "infinity", residue=0),
                ProjectivePadicBall(7, 1, "affine", residue=3),
            )
        )
        self.assertEqual(lattice.index, 35)
        self.assertTrue(all(lattice.contains(vector) for vector in lattice.basis))
        self.assertTrue(
            lattice.point_matches_charts(ProjectiveRational.normalized(1, 5))
        )
        self.assertFalse(
            lattice.point_matches_charts(ProjectiveRational.normalized(2, 1))
        )

    def test_pinned_manifest_is_deterministic_and_groundwork_only(self) -> None:
        manifest = build_manifest()
        declared = manifest.pop("result_sha256")
        self.assertEqual(canonical_digest(manifest), declared)
        manifest["result_sha256"] = declared
        self.assertEqual(manifest["priority_order"][0], "R20-RESIDUAL-2SELMER")
        self.assertEqual(
            manifest["exact_calibrations"]["R20_quotient_escape"][
                "marginal_dimension"
            ],
            8,
        )
        self.assertEqual(
            manifest["exact_calibrations"]["V4_pair_cover"]["quotient_genera"],
            [2, 2, 5],
        )
        self.assertIn("no new rank", manifest["claim_level"])

        pinned_path = (
            ROOT
            / "artifacts/generated-results/elliptic-curves/elliptic_structural_search_groundwork.json"
        )
        if pinned_path.exists():
            self.assertEqual(json.loads(pinned_path.read_text()), manifest)


if __name__ == "__main__":
    unittest.main()

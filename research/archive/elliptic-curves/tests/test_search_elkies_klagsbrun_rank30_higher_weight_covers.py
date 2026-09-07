#!/usr/bin/env python3

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from alternate_quartic_covers import mobius_preimage  # noqa: E402
from elkies_klagsbrun_rank29 import (  # noqa: E402
    from_short_point,
    published_short_points,
)
from search_elkies_klagsbrun_rank30 import (  # noqa: E402
    exact_linear_combination,
    poly_evaluate,
)
from search_elkies_klagsbrun_rank30_higher_weight_covers import (  # noqa: E402
    DEFAULT_SEED,
    WEIGHT_BANDS,
    build_candidate,
    build_cover_chart_plans,
    canonical_state,
    evaluate_state,
    map_plan_points,
    search_signed_classes,
    select_diverse_evaluations,
    signed_public_short_points,
    signed_short_sum,
)


Q = Fraction
SCRIPT = CAS / "search_elkies_klagsbrun_rank30_higher_weight_covers.py"
ARTIFACT = GENERATED / "elliptic_elkies_klagsbrun_rank30_higher_weight_covers.json"


class HigherWeightRank30CoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.state = canonical_state(
            (1 << 0) | (1 << 2) | (1 << 5) | (1 << 9),
            (1 << 0) | (1 << 5),
        )
        cls.evaluation = evaluate_state(cls.state)
        cls.candidate = build_candidate(cls.evaluation)

    def test_global_sign_is_canonical_but_complement_is_distinct(self) -> None:
        self.assertEqual(
            canonical_state(self.state.mask, self.state.negative_mask ^ self.state.mask),
            self.state,
        )
        complement = canonical_state(
            ((1 << 29) - 1) ^ self.state.mask,
            0,
        )
        self.assertNotEqual(complement.mask, self.state.mask)
        self.assertEqual(self.state.weight, 4)
        self.assertGreater(complement.weight, 3)

    def test_signed_representative_and_all_58_seed_parameters_are_exact(self) -> None:
        self.assertEqual(len(self.candidate.parameters), 58)
        self.assertEqual(signed_short_sum(self.state), self.candidate.base_short_point)
        self.assertEqual(
            exact_linear_combination(self.state.coefficients),
            from_short_point(self.candidate.base_short_point),
        )
        cover = self.candidate.cover
        for point in signed_public_short_points()[::11]:
            cover_point = cover.curve_point_to_cover(point)
            self.assertEqual(cover.value(cover_point[0]), cover_point[1] ** 2)
            self.assertEqual(cover.cover_point_to_curve(cover_point), point)

    def test_affine_and_mobius_chart_maps_replay_exactly(self) -> None:
        plans = build_cover_chart_plans(
            (self.candidate,),
            offset_count=1,
            affine_count=1,
            mobius_count=1,
            affine_pool_size=6,
            mobius_pool_size=6,
            offset_height=10,
            affine_height=10,
            mobius_height=10,
            skew_numerator_bound=100,
            skew_denominator_bound=10,
        )
        self.assertEqual(len(plans), 4)
        cover = self.candidate.cover
        short_points = published_short_points()
        for plan in plans:
            original_point = next(
                cover.curve_point_to_cover(point)
                for point in short_points
                if plan.matrix is None
                or mobius_preimage(plan.matrix, cover.curve_point_to_cover(point)[0])
                is not None
            )
            parameter, ordinate = original_point
            if plan.matrix is None:
                local = (parameter - plan.center) / plan.scale
                raw = local, ordinate
            else:
                local = mobius_preimage(plan.matrix, parameter)
                self.assertIsNotNone(local)
                a_value, b_value, c_value, d_value = plan.matrix
                del a_value, b_value
                raw = local, ordinate * (c_value * local + d_value) ** 2
            self.assertEqual(poly_evaluate(plan.polynomial, raw[0]), raw[1] ** 2)
            images, pole_count = map_plan_points(plan, self.candidate, (raw,))
            self.assertEqual(pole_count, 0)
            self.assertEqual(images, (from_short_point(cover.cover_point_to_curve(original_point)),))

    def test_small_beam_is_deterministic_stratified_and_disjoint(self) -> None:
        evaluations, bands = search_signed_classes(
            seed=DEFAULT_SEED,
            evaluation_budget=100,
            beam_width=4,
            rounds=2,
            mutations_per_state=2,
        )
        self.assertEqual(len(evaluations), 100)
        self.assertEqual(len(bands), len(WEIGHT_BANDS))
        self.assertTrue(all(item.state.weight >= 4 for item in evaluations))
        selected = select_diverse_evaluations(
            evaluations, count=5, minimum_distance=4
        )
        self.assertEqual([item.state.weight for item in selected], [4, 7, 11, 18, 26])
        self.assertEqual(
            selected[0].state.identifier,
            "w04_m00000684_n00000200",
        )

    def test_generated_artifact_if_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the higher-weight rank-30 artifact has not been generated")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertGreaterEqual(data["quotient_semantics"]["minimum_evaluated_weight"], 4)
        self.assertTrue(data["search_budget"]["one_pass_no_retry"])
        self.assertTrue(
            all(not record["retried"] for record in data["search_result"]["run_records"])
        )
        allowed = {
            "exact_seed_or_companion_in_rank29_subgroup",
            "exactly_in_published_rank29_subgroup",
            "unresolved_after_exact_relation_and_mod2_search",
            "exact_independent_30th_point",
        }
        self.assertTrue(
            all(
                record["classification"] in allowed
                for record in data["search_result"]["point_records"]
            )
        )


if __name__ == "__main__":
    unittest.main()

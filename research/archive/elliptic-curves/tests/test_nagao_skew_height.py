from fractions import Fraction
import json
from pathlib import Path
import unittest

from nagao_1994 import rank13_known_quartic_points
from nagao_linear_sections import omitted_companion_sections
from nagao_skew_height import (
    build_mobius_chart_plan,
    checkpoint_reference,
    classify_uniform_checkpoint,
    load_rank17_target,
    target_plan,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts" / "generated-results"


class ParameterizedSkewHeightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target = load_rank17_target(
            GENERATED / "elliptic_nagao_rank17_frontier_certificate.json",
            Q(135, 2),
        )

    def test_u135_target_loads_from_exact_rank17_certificate(self) -> None:
        self.assertEqual(self.target.parameter_t, Q(5065, 36))
        self.assertEqual(self.target.certified_rank_lower_bound, 17)
        self.assertEqual(len(self.target.saturated_basis), 17)
        self.assertEqual(
            self.target.log_conductor,
            "144.92745591457677468873080648944135759260092871198652625659362117960894403350375",
        )

    def test_u135_reference_defines_the_next_no_compute_plan(self) -> None:
        reference = checkpoint_reference(
            GENERATED / "elliptic_nagao_rank13_rank_gain_search.json",
            Q(135, 2),
        )
        self.assertEqual(reference.unexpected_x_count, 66)
        plan = target_plan(self.target, reference)
        self.assertFalse(plan["launches_search"])
        self.assertEqual(
            plan["mobius"]["expected_chart_count_after_checkpoint"], 132
        )
        self.assertEqual(len(plan["skew_boxes"]), 10)

    def test_checkpoint_classification_excludes_all_known_linear_sections(self) -> None:
        displayed = rank13_known_quartic_points(self.target.parameter_t)
        companions = tuple(
            section.point(self.target.parameter_t)
            for section in omitted_companion_sections()
        )
        raw = tuple(
            point
            for chosen in displayed + companions
            for point in (chosen, (chosen[0], -chosen[1]))
        )
        checkpoint = classify_uniform_checkpoint(
            self.target, raw, height_bound=1_000_000
        )
        self.assertEqual(len(checkpoint.unexpected_points), 0)
        self.assertEqual(checkpoint.displayed_x_count, 13)
        self.assertEqual(checkpoint.companion_x_count, 5)

    def test_chart_plan_is_parameter_independent(self) -> None:
        displayed = rank13_known_quartic_points(self.target.parameter_t)
        # Treat two exact displayed points as synthetic extras to exercise the
        # planner without launching PARI.  Classification is intentionally not
        # used here because it would correctly remove them.
        from nagao_skew_height import UniformCheckpoint

        checkpoint = UniformCheckpoint(
            target=self.target,
            height_bound=1,
            raw_signed_points=displayed[:2],
            signless_points=displayed[:2],
            unexpected_points=displayed[:2],
            displayed_x_count=0,
            companion_x_count=0,
            zero_ordinate_count=0,
        )
        charts = build_mobius_chart_plan(checkpoint)
        self.assertEqual(len(charts), 4)
        for _, (a_value, b_value, c_value, d_value) in charts:
            self.assertEqual(a_value * d_value - b_value * c_value, 1)


if __name__ == "__main__":
    unittest.main()

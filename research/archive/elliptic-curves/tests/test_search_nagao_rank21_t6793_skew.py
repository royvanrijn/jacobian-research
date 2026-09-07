from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from alternate_quartic_covers import mobius_preimage, point_on_short_curve  # noqa: E402
from nagao_1994 import (  # noqa: E402
    RANK21_CONSTRUCTION,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
)
from search_nagao_rank21_t6793_skew import (  # noqa: E402
    INPUT_BASIS_DIGEST,
    INPUT_RANK,
    PARAMETER_T,
    load_exact_basis,
    optimized_cross_ratio_charts,
)
from search_nagao_rank21_t956_skew import exact_linear_combination  # noqa: E402
from search_nagao_u42_skew_height import (  # noqa: E402
    SEARCH_BOXES,
    transform_binary_quartic,
)


INPUT = ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_unbiased.json"
ARTIFACT = ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t6793_skew.json"


class NagaoRank21T6793SkewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.short, cls.basis, cls.conductor = load_exact_basis(INPUT)
        cls.quartic = primitive_quartic_coefficients(
            RANK21_CONSTRUCTION, PARAMETER_T
        )
        cls.artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_exact_input_basis_and_conductor_are_pinned(self) -> None:
        self.assertEqual(len(self.basis), INPUT_RANK)
        self.assertEqual(
            self.artifact["input"]["saturated_basis_sha256"],
            INPUT_BASIS_DIGEST,
        )
        self.assertTrue(all(point_on_short_curve(self.short, point) for point in self.basis))
        self.assertEqual(self.conductor["root_number"], -1)
        self.assertEqual(
            self.conductor["log_conductor"],
            "158.572648489303177852801115450776268252943114125935864368897",
        )

    def test_cross_ratio_charts_transform_known_points_exactly(self) -> None:
        visible = primitive_visible_points(RANK21_CONSTRUCTION, PARAMETER_T)
        charts = optimized_cross_ratio_charts(
            [point[0] for point in visible], count=6
        )
        self.assertEqual(len(charts), 6)
        for chart in charts:
            a_value, b_value, c_value, d_value = chart.matrix
            self.assertNotEqual(a_value * d_value - b_value * c_value, 0)
            transformed = transform_binary_quartic(self.quartic, chart.matrix)
            checked = 0
            for x_value, z_value in visible:
                preimage = mobius_preimage(chart.matrix, x_value)
                if preimage is None:
                    continue
                transformed_z = z_value * (c_value * preimage + d_value) ** 2
                self.assertEqual(
                    transformed_z**2,
                    quartic_value(transformed, preimage),
                )
                checked += 1
            self.assertGreaterEqual(checked, len(visible) - 1)

    def test_declared_boxes_and_chart_attempts_are_audited(self) -> None:
        boxes = self.artifact["skew_staircase"]["boxes"]
        self.assertEqual(
            [record["id"] for record in boxes],
            [search_box.identifier for search_box in SEARCH_BOXES],
        )
        self.assertTrue(all(record["status"] == "completed" for record in boxes))
        self.assertTrue(all(not record["retried"] for record in boxes))
        chart_records = self.artifact["mobius_search"]["records"]
        self.assertEqual(len(chart_records), 108 + 16 + 4)
        self.assertEqual(sum(record["status"] == "completed" for record in chart_records), 124)
        self.assertEqual(sum(record["status"] == "timeout" for record in chart_records), 4)
        self.assertTrue(all(not record["retried"] for record in chart_records))
        scope = self.artifact["bounded_scope"]
        self.assertEqual(scope["completed_box_calls"], 10)
        self.assertEqual(scope["completed_chart_calls"], 124)
        self.assertEqual(scope["timed_out_chart_calls"], 4)
        self.assertFalse(scope["all_declared_runs_completed"])

    def test_all_new_images_replay_in_the_exact_rank19_span(self) -> None:
        analysis = self.artifact["new_point_analysis"]
        self.assertEqual(analysis["outside_checkpoint_quartic_abscissa_count"], 27)
        self.assertEqual(analysis["distinct_nonbasis_jacobian_sign_pairs"], 27)
        self.assertEqual(analysis["exact_relations_in_certified_rank19_span"], 27)
        self.assertEqual(analysis["unresolved_by_exact_relation_replay"], 0)
        for record in analysis["records"]:
            quartic_point = Q(record["quartic_x"]), Q(record["quartic_z"])
            jacobian_point = Q(record["jacobian_x"]), Q(record["jacobian_y"])
            self.assertEqual(
                quartic_point[1] ** 2,
                quartic_value(self.quartic, quartic_point[0]),
            )
            self.assertEqual(
                quartic_point_to_short_jacobian(
                    RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
                ),
                jacobian_point,
            )
            self.assertEqual(
                exact_linear_combination(
                    self.short, self.basis, record["basis_relation"]
                ),
                jacobian_point,
            )
        self.assertEqual(self.artifact["height_selection"]["stable_numerical_rank"], 19)
        self.assertEqual(self.artifact["exact_rank_gain_attempt"]["status"], "not_triggered")


if __name__ == "__main__":
    unittest.main()

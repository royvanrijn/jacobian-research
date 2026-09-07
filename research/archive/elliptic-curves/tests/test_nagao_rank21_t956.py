from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS_DIRECTORY = ROOT / "elliptic-curves" / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from alternate_quartic_covers import (  # noqa: E402
    point_on_short_curve,
)
from certify_nagao_rank21_t956 import (  # noqa: E402
    EXPECTED_HEIGHT_SUBSET,
    PARAMETER_T,
    point_digest,
)
from mod2_reduction_independence import (  # noqa: E402
    Mod2ReductionSignature,
    combined_mod2_rank,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)
from nagao_1994 import (  # noqa: E402
    RANK21_CONSTRUCTION,
    quartic_point_to_short_jacobian,
    quartic_value,
    primitive_quartic_coefficients,
    short_jacobian_coefficients,
)
from search_nagao_rank21_t956_skew import (  # noqa: E402
    CHART_CENTER_COUNT,
    CHART_SHIFTS,
    SEARCH_BOXES,
    chart_plan,
    exact_linear_combination,
    load_checkpoint,
)


CERTIFICATE_PATH = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank21_t956_rank17_certificate.json"
)
SKEW_PATH = (
    ROOT / "artifacts/generated-results/elliptic_nagao_rank21_t956_skew.json"
)


def signature(record: dict[str, object]) -> Mod2ReductionSignature:
    return Mod2ReductionSignature(
        prime=int(record["prime"]),
        group_order=int(record["group_order"]),
        doubled_subgroup_order=int(record["doubled_subgroup_order"]),
        quotient_dimension=int(record["quotient_dimension"]),
        rows=tuple(tuple(int(value) for value in row) for row in record["rows"]),
    )


class NagaoRank21T956Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        cls.skew = json.loads(SKEW_PATH.read_text(encoding="utf-8"))
        cls.quartic = primitive_quartic_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
        cls.short = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)

    def test_exact_rank17_certificate_replays(self) -> None:
        exact = self.certificate["exact_rank_certificate"]
        basis = tuple(
            (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
            for record in exact["saturated_basis"]
        )
        self.assertEqual(len(basis), 17)
        self.assertTrue(all(point_on_short_curve(self.short, point) for point in basis))
        signatures = tuple(signature(record) for record in exact["finite_reduction_signatures"])
        self.assertEqual(combined_mod2_rank(signatures, 17), 17)
        self.assertEqual(exact["certified_algebraic_rank_lower_bound"], 17)
        prime = exact["two_torsion_certificate_prime"]
        self.assertTrue(
            short_curve_has_no_rational_2_torsion_modular_certificate(
                self.short, prime
            )
        )
        self.assertEqual(
            self.certificate["height_selection"]["selected_pool_indices_one_based"],
            list(EXPECTED_HEIGHT_SUBSET),
        )
        self.assertEqual(
            exact["small_prime_saturation"]["height_determinant_ratio"].split(".")[0],
            str(2**32),
        )

    def test_uniform_new_points_map_exactly(self) -> None:
        records = self.certificate["uniform_search"]["new_points"]
        self.assertEqual(len(records), 43)
        mapped = []
        for record in records:
            quartic_point = Q(record["quartic_x"]), Q(record["quartic_z"])
            self.assertEqual(
                quartic_point[1] ** 2,
                quartic_value(self.quartic, quartic_point[0]),
            )
            image = quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
            )
            self.assertEqual(
                image,
                (Q(record["jacobian_x"]), Q(record["jacobian_y"])),
            )
            mapped.append(image)
        self.assertEqual(
            point_digest(mapped), self.certificate["uniform_search"]["new_image_sha256"]
        )

    def test_skew_boxes_and_chart_plan_are_exactly_the_declared_scope(self) -> None:
        boxes = self.skew["skew_staircase"]["boxes"]
        self.assertEqual([record["id"] for record in boxes], [box.identifier for box in SEARCH_BOXES])
        self.assertTrue(all(record["status"] == "completed" for record in boxes))
        self.assertTrue(all(not record["retried"] for record in boxes))
        self.assertEqual(SEARCH_BOXES[0].denominator_lower, 1)
        self.assertEqual(SEARCH_BOXES[-1].denominator_upper, 128_000)
        for left, right in zip(SEARCH_BOXES, SEARCH_BOXES[1:]):
            self.assertEqual(left.denominator_upper + 1, right.denominator_lower)

        basis, checkpoint_points, _ = load_checkpoint(CERTIFICATE_PATH)
        del basis
        plan = chart_plan(checkpoint_points)
        self.assertEqual(len(plan), CHART_CENTER_COUNT * len(CHART_SHIFTS))
        stored = self.skew["unimodular_charts"]["records"]
        self.assertEqual([record["id"] for record in stored], [item[0] for item in plan])
        for record, (_, center, matrix) in zip(stored, plan):
            self.assertEqual(Q(record["center"]), center)
            self.assertEqual(tuple(record["matrix_a_b_c_d"]), matrix)
            self.assertEqual(record["determinant"], 1)
            self.assertEqual(record["status"], "completed")
            self.assertFalse(record["retried"])

    def test_all_42_new_images_have_exact_rank17_basis_relations(self) -> None:
        basis, _, _ = load_checkpoint(CERTIFICATE_PATH)
        records = self.skew["new_point_analysis"]["records"]
        new_records = [
            record
            for record in records
            if not record["duplicate_checkpoint_or_prior_jacobian_sign_pair"]
        ]
        self.assertEqual(len(new_records), 42)
        self.assertTrue(
            self.skew["new_point_analysis"][
                "all_new_images_exactly_in_certified_rank17_span"
            ]
        )
        for record in new_records:
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
                    self.short, basis, record["certified_basis_relation"]
                ),
                jacobian_point,
            )
        self.assertEqual(self.skew["height_selection"]["stable_numerical_rank"], 17)
        self.assertEqual(self.skew["height_selection"]["stable_numerical_rank_gain"], 0)
        self.assertEqual(self.skew["exact_rank_gain_attempt"]["status"], "not_triggered")

    def test_certificate_hash_is_pinned_by_skew_artifact(self) -> None:
        digest = hashlib.sha256(CERTIFICATE_PATH.read_bytes()).hexdigest()
        self.assertEqual(self.skew["certificate_input"]["sha256"], digest)


if __name__ == "__main__":
    unittest.main()


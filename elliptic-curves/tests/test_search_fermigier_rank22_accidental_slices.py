from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CAS))

from alternate_quartic_covers import short_add  # noqa: E402
from fermigier_mestre import (  # noqa: E402
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
)
from search_fermigier_rank22_accidental_slices import (  # noqa: E402
    MODEL_CHANGE_SHORT_TO_MINIMAL,
    CONJUGATE_EXTRA_SECTION_X,
    ESCALATION_CHECKPOINT_SHA256,
    RECORD_SEARCH_RELATIVE_ACCIDENTAL_X,
    Slice,
    T0,
    build_slices,
    canonical_signless_points,
    candidate_point_pools,
    generic_group_seed_points,
    minimum_to_short,
    poly_evaluate,
    published_accidental_points,
    quartic_group_pullback,
    record_search_relative_accidentals,
    select_reconstruction_convention,
    short_negate,
    short_sum,
    short_to_minimum,
    slice_polynomial,
    slice_parameter_incidences,
)
from ek_k3 import rational_square_root  # noqa: E402
from verify_fermigier_rank22_points import PUBLISHED_POINTS  # noqa: E402


class FermigierRank22AccidentalSliceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        (
            cls.base_index,
            cls.sign,
            cls.reconstruction,
            cls.trials,
        ) = select_reconstruction_convention()
        cls.published_slices = build_slices(
            published_accidental_points(cls.reconstruction)
        )
        replay_points = []
        for x_value in RECORD_SEARCH_RELATIVE_ACCIDENTAL_X:
            ordinate = rational_square_root(
                FermigierMestreFamily.quartic_value(T0, x_value)
            )
            if ordinate is None:
                raise AssertionError("a pinned record source ceased to be rational")
            replay_points.extend(((x_value, ordinate), (x_value, -ordinate)))
        cls.record_accidentals = record_search_relative_accidentals(
            tuple(FermigierMestreFamily.known_quartic_points(T0))
            + tuple(replay_points)
        )
        cls.slices = build_slices(cls.record_accidentals)

    def test_pinned_minimal_model_change_round_trips_all_points(self) -> None:
        self.assertEqual(MODEL_CHANGE_SHORT_TO_MINIMAL[0], Q(14, 507))
        for point in PUBLISHED_POINTS:
            self.assertEqual(short_to_minimum(minimum_to_short(point)), point)

    def test_unique_convention_recovers_all_twenty_two_preimages(self) -> None:
        self.assertEqual((self.base_index, self.sign), (0, 1))
        maximum = max(trial["generic_abscissa_matches"] for trial in self.trials)
        winners = [
            trial
            for trial in self.trials
            if trial["generic_abscissa_matches"] == maximum
        ]
        self.assertEqual(len(winners), 1)
        self.assertEqual(maximum, 11)
        self.assertEqual(len(self.reconstruction), 22)

        generic = FermigierMestreFamily.known_quartic_points(T0)
        generic_x = {point[0] for point in generic}
        accidental_labels = [
            f"P{index}"
            for index, point in enumerate(self.reconstruction, start=1)
            if point[0] not in generic_x
        ]
        self.assertEqual(
            accidental_labels,
            ["P6", *[f"P{index}" for index in range(13, 23)]],
        )

        coefficients = FermigierMestreFamily.coefficients(T0)
        offset = FermigierMestreFamily.quartic_point_to_jacobian(T0, generic[0])
        for published, preimage in zip(
            PUBLISHED_POINTS, self.reconstruction, strict=True
        ):
            self.assertEqual(
                preimage[1] ** 2,
                FermigierMestreFamily.quartic_value(T0, preimage[0]),
            )
            short_point = minimum_to_short(published)
            doubled = short_add(coefficients, short_point, short_point)
            self.assertIsNotNone(doubled)
            target = short_add(coefficients, offset, doubled)
            self.assertEqual(
                FermigierMestreFamily.quartic_point_to_jacobian(T0, preimage),
                target,
            )
            self.assertEqual(quartic_group_pullback(T0, preimage), short_point)

    def test_all_priority_slices_are_exact_quartics(self) -> None:
        self.assertEqual(T0, NORMALIZED_RECORD_PARAMETER)
        self.assertEqual(len(self.slices), 28)
        self.assertEqual(len(self.record_accidentals), 14)
        self.assertEqual(
            sum(
                point[0] == CONJUGATE_EXTRA_SECTION_X
                for _, point in self.record_accidentals
            ),
            1,
        )
        self.assertEqual({slice_data.slope for slice_data in self.slices}, {-1, 1})
        for slice_data in self.slices:
            self.assertEqual(len(slice_data.coefficients) - 1, 4)
            self.assertEqual(
                poly_evaluate(slice_data.coefficients, T0),
                slice_data.source_point[1] ** 2,
            )
            self.assertEqual(slice_data.x_value(T0), slice_data.source_point[0])

    def test_exact_slice_point_generic_collision_is_filtered(self) -> None:
        slice_data = next(
            item for item in self.published_slices if item.identifier == "p6_p1"
        )
        parameter = Q(4522, 39)
        ordinate = Q(306326456219, 1521)
        self.assertEqual(
            ordinate**2,
            poly_evaluate(slice_data.coefficients, parameter),
        )
        raw = ((parameter, ordinate), (parameter, -ordinate))
        pools, provenance = candidate_point_pools(
            (slice_data,), ((raw, {"status": "completed"}),)
        )
        self.assertEqual(pools, {})
        self.assertEqual(provenance, {})

    def test_negative_parameter_incidence_keeps_raw_x(self) -> None:
        point = canonical_signless_points(
            (FermigierMestreFamily.known_quartic_points(T0)[0],)
        )[0]
        tau = -T0
        slope = 1
        intercept = point[0] - slope * tau
        slice_data = Slice(
            "TEST",
            point,
            slope,
            intercept,
            slice_polynomial(slope, intercept),
        )
        self.assertEqual(poly_evaluate(slice_data.coefficients, tau), point[1] ** 2)
        incidences = slice_parameter_incidences(slice_data, ((tau, point[1]),))
        self.assertEqual(incidences[0]["quartic_x"], str(point[0]))
        self.assertEqual(
            incidences[0]["canonical_even_fiber_slope"], -slice_data.slope
        )

    def test_generic_group_seed_coordinate_has_thirteen_classes(self) -> None:
        seeds = generic_group_seed_points(T0)
        self.assertEqual(len(seeds), 13)
        self.assertEqual(len({point[0] for point in seeds}), 13)

    def test_conjugate_extra_source_has_exact_record_fiber_dependency(self) -> None:
        generic = FermigierMestreFamily.known_quartic_points(T0)
        pullbacks = tuple(
            quartic_group_pullback(T0, point) for point in generic[1:]
        )
        self.assertNotIn(None, pullbacks)
        conjugate = next(
            point
            for _, point in self.record_accidentals
            if point[0] == CONJUGATE_EXTRA_SECTION_X
        )
        conjugate_pullback = quartic_group_pullback(T0, conjugate)
        self.assertIsNotNone(conjugate_pullback)
        relation = (
            pullbacks[5],
            short_negate(pullbacks[6]),
            short_negate(pullbacks[7]),
            pullbacks[8],
            short_negate(pullbacks[11]),
            conjugate_pullback,
        )
        self.assertIsNone(
            short_sum(FermigierMestreFamily.coefficients(T0), relation)
        )

    def test_pinned_bounded_artifact_frontier(self) -> None:
        artifact_path = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic-curves"
            / "elliptic_fermigier_rank22_accidental_slices.json"
        )
        artifact = json.loads(artifact_path.read_text())
        summary = artifact["summary"]
        self.assertEqual(summary["published_points_reconstructed_exactly"], 22)
        self.assertEqual(summary["published_accidental_preimages"], 11)
        self.assertEqual(summary["record_replay_signed_points"], 54)
        self.assertEqual(summary["record_replay_signless_abscissas"], 27)
        self.assertEqual(summary["record_replay_extra_vs_positive_catalog"], 14)
        self.assertEqual(
            summary["record_replay_specialization_accidental_sources"], 13
        )
        self.assertEqual(summary["record_replay_exact_certified_rank_lower_bound"], 22)
        self.assertFalse(summary["record_replay_meets_strict_conductor_target"])
        self.assertEqual(summary["genus_one_slices"], 28)
        self.assertEqual(summary["slice_searches_completed"], 28)
        self.assertEqual(summary["slice_search_timeouts"], 0)
        self.assertEqual(summary["distinct_new_parameters"], 14)
        self.assertEqual(summary["multi_slice_parameter_count"], 3)
        self.assertEqual(summary["maximum_slice_incidence_count"], 2)
        self.assertEqual(
            summary[
                "maximum_forced_pool_numerical_rank_before_specialized_search"
            ],
            12,
        )
        self.assertEqual(summary["maximum_exact_pool_point_count"], 15)
        self.assertEqual(summary["specialized_quartic_height_50000_screens"], 5)
        self.assertEqual(summary["specialized_quartic_height_50000_completed"], 5)
        self.assertEqual(
            summary["specialized_quartic_height_50000_ranks"],
            {
                "88893/26": 12,
                "3115/3": 15,
                "9191/30": 12,
                "11305/6": 13,
                "121919/260": 12,
            },
        )
        self.assertEqual(
            summary["specialized_quartic_frontier_parameter_t"], "3115/3"
        )
        self.assertEqual(
            summary["specialized_quartic_frontier_H1000000_rank"], 15
        )
        self.assertEqual(summary["finite_reduction_certificates_triggered"], 0)
        self.assertEqual(summary["target_hits"], [])
        self.assertEqual(
            [record["parameter_t"] for record in artifact["candidate_conductor_screen"]],
            [
                "59262/13",
                "79016/351",
                "88893/26",
                "3115/3",
                "5614/201",
                "9191/135",
                "9191/30",
                "11305/6",
                "46687/260",
                "96663/26",
                "121919/1170",
                "121919/260",
                "154357/162",
                "154357/48",
            ],
        )
        certificate = artifact["record_fiber_height_1000000_replay"][
            "exact_rank_certificate"
        ]["canonical_published_point_finite_reduction_certificate"]
        self.assertEqual(certificate["status"], "certified")
        self.assertEqual(certificate["combined_exact_rank_over_F2"], 22)
        self.assertEqual(certificate["two_torsion_certificate_prime"], 31)
        self.assertEqual(
            certificate["certificate_primes"],
            [29, 43, 67, 73, 79, 83, 89, 101, 103, 107, 109, 127, 131, 137, 149, 191, 223],
        )
        frontier = next(
            record
            for record in artifact["candidate_conductor_screen"]
            if record["parameter_t"] == "3115/3"
        )["specialized_quartic_height_50000_screen"]
        self.assertEqual(frontier["signed_point_count"], 46)
        self.assertEqual(frontier["new_distinct_direct_images_modulo_sign"], 10)
        self.assertEqual(frontier["height_rank"]["stable_numerical_rank"], 15)
        escalation = frontier["height_1000000_escalation"]
        self.assertEqual(escalation["signed_quartic_points_found"], 134)
        self.assertEqual(escalation["new_x_values_beyond_visible_sections"], 54)
        self.assertEqual(escalation["stable_numerical_rank"], 15)
        checkpoint_path = (
            ROOT
            / "archive"
            / "elliptic-curves"
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_3115_3_h1000000_checkpoint.json"
        )
        self.assertEqual(
            hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
            ESCALATION_CHECKPOINT_SHA256,
        )
        script = CAS / "search_fermigier_rank22_accidental_slices.py"
        self.assertEqual(
            hashlib.sha256(script.read_bytes()).hexdigest(),
            artifact["script_sha256"],
        )


if __name__ == "__main__":
    unittest.main()

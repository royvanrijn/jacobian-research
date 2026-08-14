from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    RANK13_BASE_CHANGE_CONSTANT,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
)
from search_nagao_rank13_rank_gain import (  # noqa: E402
    DEFAULT_ANCHORS,
    DEFAULT_CENTERS,
    canonical_positive_u,
    companion_section_x_values,
    generate_population,
)
from triage_nagao_rank13_finalists import point_on_short_curve  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts" / "generated-results"


class NagaoRankGainSearchTests(unittest.TestCase):
    def test_u_symmetry_is_exact_and_canonical(self) -> None:
        constant = Q(RANK13_BASE_CHANGE_CONSTANT)
        for parameter_u in (Q(42), Q(471, 11), Q(200)):
            with self.subTest(parameter_u=parameter_u):
                partner = constant / parameter_u
                self.assertEqual(
                    rank13_base_changed_short_jacobian_coefficients(parameter_u),
                    rank13_base_changed_short_jacobian_coefficients(partner),
                )
                self.assertEqual(
                    canonical_positive_u(parameter_u),
                    canonical_positive_u(partner),
                )
        with self.assertRaises(ValueError):
            canonical_positive_u(Q(0))

    def test_population_is_reduced_deduplicated_and_keeps_anchors(self) -> None:
        population = generate_population(
            farey_denominator=3,
            mutation_denominator=4,
            mutation_numerator_radius=2,
            centers=DEFAULT_CENTERS,
            anchors=DEFAULT_ANCHORS,
        )
        parameters = [candidate.parameter_u for candidate in population]
        self.assertEqual(parameters, sorted(set(parameters)))
        self.assertTrue(all(value**2 < RANK13_BASE_CHANGE_CONSTANT for value in parameters))
        by_u = {candidate.parameter_u: candidate for candidate in population}
        for anchor in DEFAULT_ANCHORS:
            self.assertIn(canonical_positive_u(anchor), by_u)
            self.assertIn(
                "forced-anchor", by_u[canonical_positive_u(anchor)].origins
            )

    def test_five_companions_are_present_in_the_old_u42_raw_label(self) -> None:
        data = json.loads(
            (GENERATED / "elliptic_nagao_rank13_finalist_triage.json").read_text()
        )
        record = next(
            candidate for candidate in data["candidates"]
            if int(candidate["parameter_u"]) == 42
        )
        old_new_x = {
            Q(point["quartic_x"])
            for point in record["bounded_search"]["new_point_records"]
        }
        companions = set(companion_section_x_values(rank13_base_parameter(Q(42))))
        self.assertEqual(len(companions), 5)
        self.assertTrue(companions <= old_new_x)

    def test_artifact_pins_two_new_rank17_numerical_frontiers(self) -> None:
        artifact_path = GENERATED / "elliptic_nagao_rank13_rank_gain_search.json"
        data = json.loads(artifact_path.read_text())
        self.assertEqual(data["population"]["count"], 9196)
        self.assertTrue(data["method"]["conductor_not_used_for_selection"])
        self.assertEqual(data["summary"]["maximum_stable_numerical_rank"], 17)
        records = {
            Q(record["parameter_u"]): record
            for record in data["escalation_box"]["records"]
        }
        for parameter_u in (Q(471, 11), Q(135, 2)):
            with self.subTest(parameter_u=parameter_u):
                record = records[parameter_u]
                self.assertEqual(record["quartic_naive_height_bound"], 1_000_000)
                self.assertTrue(record["height_rank_stable_across_precisions"])
                self.assertEqual(record["stable_baseline_numerical_rank"], 13)
                self.assertEqual(record["stable_pool_numerical_rank"], 17)
                self.assertEqual(record["stable_numerical_rank_gain"], 4)
                self.assertTrue(
                    record["conductor_probe"][
                        "below_strict_log_conductor_target"
                    ]
                )
                subset = record["explicit_numerically_independent_subset"]
                self.assertEqual(len(subset), 17)
                coefficients = rank13_base_changed_short_jacobian_coefficients(
                    parameter_u
                )
                points = tuple(
                    (Q(point["jacobian_x"]), Q(point["jacobian_y"]))
                    for point in subset
                )
                self.assertTrue(
                    all(point_on_short_curve(coefficients, point) for point in points)
                )

        script_path = CAS / "search_nagao_rank13_rank_gain.py"
        self.assertEqual(
            data["script_sha256"], hashlib.sha256(script_path.read_bytes()).hexdigest()
        )

    def test_mutation_artifact_pins_u74_rank17_numerical_frontier(self) -> None:
        artifact_path = GENERATED / "elliptic_nagao_rank13_rank_gain_mutations.json"
        data = json.loads(artifact_path.read_text())
        self.assertEqual(data["population"]["count"], 5133)
        self.assertEqual(data["population"]["mutation_denominator_bound"], 64)
        record = next(
            item
            for item in data["escalation_box"]["records"]
            if Q(item["parameter_u"]) == Q(74)
        )
        self.assertEqual(record["quartic_naive_height_bound"], 1_000_000)
        self.assertTrue(record["height_rank_stable_across_precisions"])
        self.assertEqual(record["stable_baseline_numerical_rank"], 13)
        self.assertEqual(record["stable_pool_numerical_rank"], 17)
        self.assertTrue(
            record["conductor_probe"]["below_strict_log_conductor_target"]
        )
        subset = record["explicit_numerically_independent_subset"]
        self.assertEqual(len(subset), 17)
        coefficients = rank13_base_changed_short_jacobian_coefficients(Q(74))
        self.assertTrue(
            all(
                point_on_short_curve(
                    coefficients,
                    (Q(point["jacobian_x"]), Q(point["jacobian_y"])),
                )
                for point in subset
            )
        )
        script_path = CAS / "search_nagao_rank13_rank_gain.py"
        self.assertEqual(
            data["script_sha256"], hashlib.sha256(script_path.read_bytes()).hexdigest()
        )


if __name__ == "__main__":
    unittest.main()

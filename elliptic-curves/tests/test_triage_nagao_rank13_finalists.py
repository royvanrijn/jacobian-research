from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    RANK13_CONSTRUCTION,
    quartic_point_to_short_jacobian,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
)
from triage_nagao_rank13_finalists import (  # noqa: E402
    load_finalists,
    parse_u_values,
    point_digest,
    point_on_short_curve,
    split_infinity_jacobian_point,
    stable_height_rank,
)
from triage_nagao_rank13_local_candidates import (  # noqa: E402
    load_local_integer_candidates,
)
from extend_nagao_u42_frontier import checkpoint_points  # noqa: E402
from extend_nagao_u118_height import validate_checkpoint  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]


class NagaoFinalistTriageTests(unittest.TestCase):
    def test_split_infinity_limit_is_an_exact_point(self) -> None:
        for parameter_u in (Q(1), Q(50), Q(84), Q(1256)):
            with self.subTest(parameter_u=parameter_u):
                coefficients = rank13_base_changed_short_jacobian_coefficients(
                    parameter_u
                )
                positive = split_infinity_jacobian_point(parameter_u)
                negative = split_infinity_jacobian_point(parameter_u, sign=-1)
                self.assertEqual(positive[0], negative[0])
                self.assertEqual(positive[1], -negative[1])
                self.assertTrue(point_on_short_curve(coefficients, positive))
        with self.assertRaises(ValueError):
            split_infinity_jacobian_point(Q(84), sign=0)

    def test_thirteen_affine_sections_map_exactly(self) -> None:
        parameter_u = Q(84)
        parameter_t = rank13_base_parameter(parameter_u)
        sections = rank13_known_quartic_points(parameter_t)
        images = tuple(
            quartic_point_to_short_jacobian(
                RANK13_CONSTRUCTION, parameter_t, point
            )
            for point in sections
        )
        coefficients = rank13_base_changed_short_jacobian_coefficients(parameter_u)
        self.assertEqual(len(sections), 13)
        self.assertEqual(len(set(sections)), 13)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in images))
        self.assertEqual(len(point_digest(sections)), 64)
        self.assertEqual(point_digest(sections), point_digest(tuple(sections)))

    def test_primary_and_scaled_frontiers_merge_without_duplicates(self) -> None:
        primary = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_integer_u.json"
        )
        scaled = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_integer_u2000.json"
        )
        finalists = load_finalists(primary, scaled, (1256, 42))
        parameters = [candidate.parameter_u for candidate in finalists]
        self.assertEqual(len(parameters), 13)
        self.assertEqual(len(set(parameters)), 13)
        self.assertEqual(parameters[-1], 1256)
        candidate = finalists[-1]
        self.assertEqual(candidate.parameter_t, Q(-4949, 8))
        self.assertTrue(candidate.below_target)
        self.assertEqual(candidate.root_number, -1)
        self.assertEqual(candidate.last_numerical_prime, 19_997)

    def test_parser_and_stability_guards(self) -> None:
        self.assertEqual(parse_u_values("84,2,189"), (84, 2, 189))
        self.assertEqual(parse_u_values(""), ())
        for value in ("0", "84,84", "x"):
            with self.assertRaises(argparse.ArgumentTypeError):
                parse_u_values(value)
        runs = (
            {"numerical_rank": 16, "subset_indices_one_based": [1, 3]},
            {"numerical_rank": 16, "subset_indices_one_based": [1, 3]},
        )
        self.assertEqual(stable_height_rank(runs), 16)
        with self.assertRaises(AssertionError):
            stable_height_rank(
                (runs[0], {"numerical_rank": 15, "subset_indices_one_based": [1, 3]})
            )

    def test_local_crt_integer_candidate_loader(self) -> None:
        path = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_local_crt.json"
        )
        candidates = load_local_integer_candidates(path, (118, 316))
        self.assertEqual([candidate.parameter_u for candidate in candidates], [118, 316])
        self.assertEqual(candidates[0].parameter_t, Q(4813, 118))
        self.assertEqual(candidates[0].root_number, -1)
        self.assertTrue(candidates[0].below_target)
        self.assertEqual(candidates[1].parameter_t, Q(-38153, 316))
        with self.assertRaises(ValueError):
            load_local_integer_candidates(path, (999_999,))

    def test_u42_checkpoint_and_saturated_basis_are_exact(self) -> None:
        checkpoint_path = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_finalist_triage.json"
        )
        checkpoint = checkpoint_points(checkpoint_path)
        self.assertEqual(len(checkpoint), 17)
        extension_path = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_u42_height_10000000.json"
        )
        extension = json.loads(extension_path.read_text())
        saturated = tuple(
            (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
            for record in extension["small_prime_saturation"]["saturated_basis"]
        )
        coefficients = rank13_base_changed_short_jacobian_coefficients(Q(42))
        self.assertEqual(len(saturated), 17)
        self.assertTrue(
            all(point_on_short_curve(coefficients, point) for point in saturated)
        )
        self.assertEqual(
            extension["small_prime_saturation"]["height_determinant_ratio"].split(
                ".", 1
            )[0],
            str(2**32),
        )
        self.assertEqual(extension["pari_ellrank_effort_zero"]["status"], "timeout")

    def test_u118_checkpoint_is_pinned_before_extension(self) -> None:
        checkpoint_path = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank13_local_candidate_triage.json"
        )
        validate_checkpoint(checkpoint_path, expected_parameter_t="4813/118")
        with self.assertRaises(ValueError):
            validate_checkpoint(checkpoint_path, expected_parameter_t="0")


if __name__ == "__main__":
    unittest.main()

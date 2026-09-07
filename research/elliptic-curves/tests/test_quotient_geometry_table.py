from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
TABLE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/quotient_geometry_table_v1.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class QuotientGeometryTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.table = json.loads(TABLE.read_text())

    def test_inventory_and_claim_boundary(self) -> None:
        table = self.table
        self.assertEqual(table["status"], "PASS_COMPLETE_QUOTIENT_GEOMETRY_TABLE")
        self.assertEqual(
            table["observation_unit"],
            {
                "usable_known_R17_controls": 5,
                "new_R17_ladder_fibres": 16,
                "MW16_A1_presentations": 9,
                "total_presentations": 30,
            },
        )
        self.assertEqual(len(table["cases"]), 30)
        self.assertEqual(
            sum(
                len(case["blind_recovery"]["direction_witnesses"])
                for case in table["cases"]
            ),
            230,
        )
        self.assertEqual(
            table["answer"]["successive_minimum_cutoff_hypothesis"], "NO"
        )
        self.assertEqual(
            table["answer"]["strictly_partial_nonempty_initial_stage_count"], 13
        )
        self.assertEqual(
            table["answer"][
                "strictly_partial_nonempty_initial_stage_counterexample_count"
            ],
            13,
        )
        self.assertEqual(
            set(table["answer"]["counterexample_case_ids"]),
            {
                "r17-control-rank28",
                "r17-ladder-curve535",
                "mw16-a1-curve398-p16875",
                "mw16-a1-curve398-p63669",
                "mw16-a1-curve400-p53042",
                "mw16-a1-curve400-p62992",
                "mw16-a1-curve401-p57487",
            },
        )

    def test_every_case_has_complete_quotient_geometry(self) -> None:
        for case in self.table["cases"]:
            geometry = case["quotient_geometry"]
            quotient_rank = int(geometry["quotient_rank"])
            self.assertEqual(case["j_displayed"], quotient_rank, case["case_id"])
            self.assertEqual(len(geometry["quotient_gram"]), quotient_rank)
            self.assertTrue(
                all(len(row) == quotient_rank for row in geometry["quotient_gram"])
            )
            self.assertGreater(float(geometry["regulator"]), 0.0)
            minima = list(map(float, geometry["successive_minima"]))
            self.assertEqual(len(minima), quotient_rank)
            self.assertEqual(minima, sorted(minima))
            self.assertEqual(
                len(geometry["successive_minimum_witnesses"]), quotient_rank
            )
            self.assertTrue(geometry["enumeration"]["complete_for_successive_minima"])

    def test_every_recovered_direction_has_the_requested_decomposition(self) -> None:
        for case in self.table["cases"]:
            for direction in case["blind_recovery"]["direction_witnesses"]:
                for field in (
                    "intrinsic_q_t_P",
                    "projection_away_from_inherited_span",
                    "optimal_half_lattice_position",
                    "half_lattice_phase",
                    "coordinate_distortion_term",
                    "blind_recovery_stage",
                    "predicted_search_height_window",
                ):
                    self.assertIn(field, direction)
                self.assertGreaterEqual(float(direction["intrinsic_q_t_P"]), 0.0)
                projection = direction["projection_away_from_inherited_span"]
                self.assertAlmostEqual(
                    float(projection["schur_complement_lambda"]),
                    float(direction["intrinsic_q_t_P"]),
                    places=7,
                )
                optimum = direction["optimal_half_lattice_position"]
                self.assertTrue(optimum["stable_center_across_rounding_scales"])
                self.assertEqual(len(optimum["rounded_cvp_candidates"]), 2)
                self.assertGreaterEqual(float(optimum["minimum_phase_energy"]), 0.0)
                self.assertLess(
                    abs(float(optimum["completing_square_residual"])), 1e-10
                )
                self.assertGreaterEqual(
                    float(direction["half_lattice_phase"]["energy"]), 0.0
                )
                window = direction["predicted_search_height_window"]
                self.assertTrue(window["centered_quotient_energy_lies_under_window"])
                self.assertLess(abs(float(window["decomposition_residual"])), 1e-10)
                self.assertLessEqual(
                    direction["coordinate_geometry"]["reduced_coordinate_height"],
                    window["height_bound_B"],
                )

    def test_curve398_separates_intrinsic_height_from_chart_visibility(self) -> None:
        diagnostic = self.table["answer"]["curve398_diagnostic"]
        self.assertEqual(diagnostic["displayed_quotient_rank"], 14)
        self.assertEqual(diagnostic["initial_recovered_rank"], 5)
        self.assertEqual(
            diagnostic[
                "successive_minimum_witness_indices_in_recovered_subspace_one_based"
            ],
            [5, 8],
        )
        self.assertGreater(
            float(diagnostic["recovered_direction_lambda_range"][1]),
            float(diagnostic["successive_minima_range"][1]),
        )
        self.assertEqual(
            self.table["summary"]["stable_projection_half_lattice_cvp_count"],
            230,
        )

    def test_fail_closed_subspace_comparisons(self) -> None:
        by_id = {case["case_id"]: case for case in self.table["cases"]}
        unresolved = {
            case_id
            for case_id, case in by_id.items()
            if not case["successive_minimum_recovery_test"][
                "displayed_containment_complete"
            ]
        }
        self.assertEqual(
            unresolved,
            {
                "r17-ladder-curve478",
                "r17-ladder-curve539",
                "mw16-a1-curve542-p30486",
            },
        )
        for case_id in unresolved:
            comparison = by_id[case_id]["successive_minimum_recovery_test"]
            self.assertIsNone(comparison["initial_equals_successive_minimum_prefix"])
            self.assertIsNone(comparison["final_equals_successive_minimum_prefix"])

    def test_input_hash_chain(self) -> None:
        for name, expected in self.table["input_hashes"].items():
            self.assertEqual(digest(ROOT / name), expected, name)


if __name__ == "__main__":
    unittest.main()

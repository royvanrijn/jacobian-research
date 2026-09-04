from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import quotient_rank_escape_detector_v2 as detector  # noqa: E402


def load_sample_builder():
    path = (
        ROOT
        / "elkies-k3/scripts"
        / "build_r17_quotient_rank_escape_detector_v2_sample.py"
    )
    spec = importlib.util.spec_from_file_location("detector_v2_sample", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class QuotientRankEscapeDetectorV2Tests(unittest.TestCase):
    def test_canonical_rref_and_kernel(self) -> None:
        rows = [[1, 1, 0, 0], [0, 1, 1, 0], [1, 0, 1, 0]]
        reduced = detector.canonical_rref(rows, 4)
        self.assertEqual(reduced, [[1, 0, 1, 0], [0, 1, 1, 0]])
        kernel = detector.nullspace_basis(rows, 4)
        self.assertEqual(kernel, [[1, 1, 1, 0], [0, 0, 0, 1]])
        self.assertTrue(
            all(detector.dot(row, vector) == 0 for row in rows for vector in kernel)
        )

    def test_exact_record_calibration_and_leave_one_out(self) -> None:
        width = 32
        unit = lambda index: [int(column == index) for column in range(width)]
        local = {
            "2": [unit(0)],
            "3": [unit(1)],
            "5": [unit(2)],
            "infinity": [],
        }
        mw17 = [unit(index) for index in range(3, 20)]
        exceptional = [unit(index) for index in range(20, 32)]
        result = detector.analyze_complete_descent(
            ambient_dimension=width,
            local_condition_rows=local,
            required_places=list(local),
            mw17_rows=mw17,
            exceptional_rows=exceptional,
            local_metadata={
                place: {"local_kummer_image_dimension": 1}
                for place in local
            },
        )
        self.assertEqual(result["two_selmer_dimension"], 29)
        self.assertEqual(result["actual_mw17_image_dimension"], 17)
        self.assertEqual(result["s_res"], 12)
        calibration = result["held_out_record_calibration"]
        self.assertEqual(calibration["exceptional_image_dimension_modulo_mw17"], 12)
        self.assertEqual(calibration["actual_mw29_image_dimension"], 29)
        self.assertEqual(calibration["selmer_modulo_mw29_dimension"], 0)
        self.assertTrue(calibration["exact_rank_29_if_selmer_dimension_29"])
        suppression = {
            row["place"]: row["single_place_suppression_of_residual_intersection"]
            for row in result["places"]
        }
        self.assertEqual(suppression, {"2": 1, "3": 1, "5": 1, "infinity": 0})
        self.assertIsNone(result["pairing_claim"])

    def test_rejects_known_point_outside_selmer(self) -> None:
        width = 18
        unit = lambda index: [int(column == index) for column in range(width)]
        with self.assertRaises(detector.DetectorInputError):
            detector.analyze_complete_descent(
                ambient_dimension=width,
                local_condition_rows={"2": [unit(0)], "infinity": []},
                required_places=["2", "infinity"],
                mw17_rows=[unit(index) for index in range(17)],
            )

    def test_rejects_omitted_declared_bad_place(self) -> None:
        width = 18
        unit = lambda index: [int(column == index) for column in range(width)]
        with self.assertRaises(detector.DetectorInputError):
            detector.analyze_complete_descent(
                ambient_dimension=width,
                local_condition_rows={"2": [unit(0)], "infinity": []},
                required_places=["2", "37", "infinity"],
                mw17_rows=[unit(index) for index in range(1, 18)],
            )

    def test_checkpointed_condition_conversion(self) -> None:
        identity = [[int(row == column) for row in range(4)] for column in range(4)]
        selmer = {
            "two_selmer_dimension": 2,
            "finite_local_condition_primes": ["2"],
            "local_condition_matrix": {
                "global_s_squareclass_dimension": 4,
                "global_norm_square_subspace_dimension": 4,
                "global_norm_square_subspace_basis_columns_in_s_squareclasses": identity,
                "selmer_basis_columns_in_global_norm_square_subspace": [identity[2], identity[3]],
                "places": [
                    {
                        "place": "2",
                        "allowed_subspace_basis_columns_in_global_s_squareclasses": identity[1:],
                        "norm_subspace_intersection_dimension_for_this_place_alone": 3,
                        "ambient_local_kummer_dimension": 1,
                        "computed_local_kummer_image_dimension": 1,
                        "localized_global_s_squareclass_image_dimension": 3,
                    },
                    {
                        "place": "infinity",
                        "allowed_subspace_basis_columns_in_global_s_squareclasses": identity,
                        "norm_subspace_intersection_dimension_for_this_place_alone": 4,
                        "ambient_local_kummer_dimension": 0,
                        "computed_local_kummer_image_dimension": 0,
                        "localized_global_s_squareclass_image_dimension": 0,
                    },
                ],
            },
        }
        rows = detector.checkpointed_simon_condition_rows(selmer)
        self.assertEqual(rows["2"], [[1, 0, 0, 0]])
        self.assertEqual(rows["infinity"], [])
        point_rows = detector.checkpointed_point_rows_in_normspace(
            selmer, {"point_selmer_rows": [[1, 1]]}
        )
        self.assertEqual(point_rows, [[0, 0, 1, 1]])

    def test_sample_is_balanced_and_blinded(self) -> None:
        blinded, key = load_sample_builder().build()
        stage1 = [row for row in key["rows"] if row["stage_1_included"]]
        stage2 = [row for row in key["rows"] if row["stage_2_included"]]
        self.assertEqual((len(stage1), len(stage2)), (10, 30))
        self.assertEqual(
            {name: sum(row["aggregate_cohort"] == name for row in stage1) for name in {
                row["aggregate_cohort"] for row in stage1
            }},
            {
                "full_cylinders": 2,
                "matched_ordinary": 2,
                "two_only": 2,
                "odd_only": 2,
                "random_equal_codimension": 2,
            },
        )
        self.assertTrue(
            all(
                set(row)
                == {
                    "blind_id",
                    "chart",
                    "parameter",
                    "projective_pair",
                    "stage_1_included",
                    "stage_2_included",
                }
                for row in blinded["rows"]
            )
        )

    def test_outcome_d_control_certificate_is_fail_closed(self) -> None:
        path = (
            ROOT
            / "artifacts/generated-results"
            / "elkies-k3-r17-quotient-rank-escape-detector-v2-controls-v1.json"
        )
        document = json.loads(path.read_text())
        self.assertEqual(
            document["status"],
            "OUTCOME_D_COMPLETE_DESCENTS_BLOCKED_RECORD_INPUTS_AND_LOCAL_DATA_CERTIFIED",
        )
        self.assertFalse(document["detector_v2"]["control_gate_passed"])
        self.assertFalse(document["detector_v2"]["stage_1_application_authorized"])
        for control, bad_count in zip(document["controls"], (13, 16)):
            self.assertEqual(
                control["specialized_mw17"]["actual_global_mod2_image_dimension"],
                17,
            )
            self.assertTrue(
                control["specialized_mw17"]["two_saturated_inside_full_E_Q"]
            )
            self.assertIsNone(
                control["specialized_mw17"]["all_prime_saturation_inside_full_E_Q"]
            )
            self.assertEqual(
                control["known_mw29"]["exceptional_image_dimension_modulo_mw17"],
                12,
            )
            self.assertEqual(
                len(control["specialized_mw17"]["global_kummer_half_ideals"]),
                17,
            )
            self.assertEqual(
                len(control["known_mw29"]["exceptional_global_kummer_representatives"]),
                12,
            )
            self.assertEqual(
                len(control["known_mw29"]["exceptional_kummer_half_ideals"]),
                12,
            )
            finite_places = control["local_conditions"]["all_bad_finite_places"]
            self.assertEqual(len(finite_places), bad_count)
            self.assertTrue(
                all(
                    len(place["exceptional_generator_local_data"]) == 12
                    and place["local_factor_descriptors"]
                    for place in finite_places
                )
            )
            self.assertIsNone(control["complete_two_selmer"]["dimension"])
            self.assertIsNone(control["complete_two_selmer"]["s_res"])


if __name__ == "__main__":
    unittest.main()

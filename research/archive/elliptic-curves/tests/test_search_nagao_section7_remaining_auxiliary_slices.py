from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_section7_remaining_auxiliary_slices.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_section7_remaining_auxiliary_slices.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_section7_remaining_auxiliary_slices", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class RemainingAuxiliarySliceTests(unittest.TestCase):
    def test_scope_is_the_exact_eleven_slice_complement(self) -> None:
        self.assertEqual(len(MODULE.OPEN_SLICE_IDS), 11)
        self.assertEqual(len(MODULE.COORDINATION_EXCLUDED_SLICE_IDS), 5)
        self.assertEqual(
            set(MODULE.OPEN_SLICE_IDS)
            | set(MODULE.COORDINATION_EXCLUDED_SLICE_IDS),
            {
                "a01_sp01",
                "a02_sp01",
                "a03_sp01",
                "a04_sp01",
                "a05_sp01",
                "a06_sp01",
                "a07_sp01",
                "a08_sm01",
                "a09_sm01",
                "a10_sp01",
                "a11_sm01",
                "a12_sp01",
                "a13_sp01",
                "a14_sp01",
                "a15_sp01",
                "a16_sp01",
            },
        )
        prepared, metadata = MODULE.load_open_slices(MODULE.DEFAULT_INPUT)
        self.assertEqual(
            tuple(item.source.identifier for item in prepared),
            MODULE.OPEN_SLICE_IDS,
        )
        self.assertTrue(metadata["all_associations_replayed_by_exact_square_tests"])
        self.assertEqual(metadata["pinned_identity_height"], 200_000)
        for item in prepared:
            self.assertEqual(
                item.source.slope * MODULE.T0 + item.source.intercept,
                item.source.source_point[0],
            )
            self.assertEqual(
                item.auxiliary.quartic_value(MODULE.T0),
                item.source.source_point[1] ** 2,
            )

    def test_ternary_cap_and_fast_generic_filter(self) -> None:
        self.assertEqual(
            [len(MODULE.ternary_vectors(dimension)) for dimension in range(7)],
            [0, 2, 8, 26, 80, 242, 728],
        )
        self.assertEqual(3**MODULE.MAX_BASIS_DIMENSION, 59_049)

        from search_nagao_rank21_accidental_slices import generic_quartic_points

        for label, point in generic_quartic_points(MODULE.T0)[:4]:
            labels = MODULE.generic_labels_for_point(MODULE.T0, point)
            self.assertIn(f"quartic-x:{label}", labels)
            self.assertIn(f"jacobian-sign-pair:{label}", labels)

    def test_generated_artifact_closes_every_declared_box(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"], MODULE.file_sha256(SCRIPT)
        )
        self.assertEqual(data["scope"]["open_slice_ids"], list(MODULE.OPEN_SLICE_IDS))
        ranks = {
            record["id"]: record["stable_numerical_height_rank"]
            for record in data["slices"]
        }
        self.assertEqual(
            ranks,
            {
                "a02_sp01": 6,
                "a03_sp01": 0,
                "a05_sp01": 3,
                "a06_sp01": 4,
                "a07_sp01": 5,
                "a11_sm01": 3,
                "a12_sp01": 2,
                "a13_sp01": 2,
                "a14_sp01": 5,
                "a15_sp01": 2,
                "a16_sp01": 6,
            },
        )
        self.assertEqual(
            sum(record["ternary_vector_count_excluding_zero"] for record in data["slices"]),
            2_096,
        )
        self.assertTrue(
            all(
                not any(record["generation_counts"][key] for key in (
                    "identity_or_torsion_relations",
                    "exceptional_inverse_images",
                    "zero_parameters",
                    "prior_H200000_or_T0_parameters",
                    "generic_section_images",
                    "singular_section7_parameters",
                    "new_parameters_at_height_at_most_200000",
                ))
                for record in data["slices"]
            )
        )
        generation = data["generation"]
        self.assertEqual(generation["unique_accepted_parameters"], 2_096)
        self.assertEqual(
            generation["candidate_parameter_sha256"],
            "66344a5e5be263613548ccff9c0f286e543f266b93fdf8683d2ecc149de3dadf",
        )
        self.assertEqual(len(generation["all_candidate_records"]), 2_096)
        proxy = data["proxy_filter"]
        self.assertGreater(proxy["minimum_log_radical_upper_proxy"], 1_000)
        self.assertEqual(proxy["below_threshold_count"], 0)
        self.assertEqual(data["exact_conductors"]["attempted"], 0)
        self.assertEqual(data["exact_conductors"]["sub_182_72_count"], 0)
        self.assertFalse(data["outcome"]["rank21_certified"])
        self.assertFalse(data["outcome"]["breakthrough_curve_found"])


if __name__ == "__main__":
    unittest.main()

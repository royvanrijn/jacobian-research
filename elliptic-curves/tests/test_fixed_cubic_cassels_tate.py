#!/usr/bin/env sage
"""Arithmetic and fail-closed regressions; run with sage -python."""

import copy
import gzip
import importlib.machinery
import importlib.util
import json
from pathlib import Path
import unittest

from sage.all import GF, QQ, matrix

ROOT = Path(__file__).resolve().parents[2]
LOADER = importlib.machinery.SourceFileLoader(
    "fixed_cubic_ct_verifier",
    str(ROOT / "elliptic-curves/cas/verify_fixed_cubic_cassels_tate.sage"),
)
SPEC = importlib.util.spec_from_loader(LOADER.name, LOADER)
VERIFIER = importlib.util.module_from_spec(SPEC)
LOADER.exec_module(VERIFIER)


class FixedCubicCasselsTateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = json.loads(gzip.decompress(VERIFIER.EVIDENCE.read_bytes()))
        cls.result = VERIFIER.verify(cls.evidence)
        cls.source, cls.local_run, cls.A, cls.B, cls.I, cls.J = VERIFIER.source_data()
        cls.covers = {row["anchor_mask"]: VERIFIER.R(list(map(QQ, row["quartic"])))
                      for row in cls.evidence["covers"]}

    def verify_pair(self, row):
        return VERIFIER.verify_pair(row, self.covers, self.local_run["complete_finite_place_support"],
                                    self.I, self.J, set())

    def test_complete_arithmetic_replay_has_rank_sixteen_and_radical_two(self):
        self.assertEqual(self.result["verified_pairing_entry_count"], 153)
        self.assertEqual(self.result["pairing_rank"], 16)
        self.assertEqual(self.result["restricted_radical_dimension"], 2)
        self.assertEqual(self.result["obstructed_class_count"], 262140)
        self.assertEqual(self.result["nonzero_compatible_class_count"], 3)
        self.assertIsNone(self.result["full_curve_rank_upper"])
        self.assertFalse(self.result["full_selmer_radical_computed"])

    def test_published_nonzero_example(self):
        VERIFIER.verify_published_control(self.evidence["published_control"])
        self.assertTrue(self.result["published_nonzero_control_verified"])

    def test_combinations_can_survive_when_no_basis_vector_does(self):
        M = matrix(GF(2), [[0, 1, 1], [1, 0, 0], [1, 0, 0]])
        result = VERIFIER.decompose(M, [{"mask": 1 << i} for i in range(3)])
        self.assertEqual(result["restricted_radical_basis"][0]["W_coordinates"], [0, 1, 1])
        self.assertEqual(result["pairing_rank"], 2)

    def test_rational_squareclasses_keep_negative_valuations_exact(self):
        self.assertTrue(VERIFIER.local_square(QQ(1)/4, 2))
        self.assertFalse(VERIFIER.local_square(QQ(5)/4, 2))
        self.assertEqual(VERIFIER.hilbert_symbol(QQ(-1)/4, -1, 2), -1)

    def test_changing_one_entry_cannot_certify_a_different_matrix(self):
        row = copy.deepcopy(self.evidence["pairings"][0])
        row["value"] ^= 1
        with self.assertRaisesRegex(ArithmeticError, "pairing entry"):
            self.verify_pair(row)

    def test_missing_place_does_not_mean_zero_contribution(self):
        row = copy.deepcopy(self.evidence["pairings"][0])
        row["local_terms"].pop(0)
        with self.assertRaisesRegex(ArithmeticError, "support"):
            self.verify_pair(row)

    def test_fabricated_square_root_is_rejected(self):
        row = copy.deepcopy(self.evidence["pairings"][0])
        row["square_root_phi_coefficients"] = ["1", "0", "0"]
        with self.assertRaisesRegex(ArithmeticError, "square-root"):
            self.verify_pair(row)

    def test_quartic_with_wrong_map_is_rejected(self):
        row = copy.deepcopy(self.evidence["covers"][0])
        row["d_over_quartic_y"] = str(2*QQ(row["d_over_quartic_y"]))
        with self.assertRaisesRegex(ArithmeticError, "two-quadric cover"):
            VERIFIER.verify_cover(row, self.source, self.A, self.B, self.I, self.J)

    def test_bounded_radical_misses_remain_unknown(self):
        rows = self.result["radical_search_candidates_by_coefficient_height"]
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["point_or_sha_status"] == "UNKNOWN" for row in rows))
        self.assertEqual([row["affine_point_count"] for row in rows], [0, 0, 0])


if __name__ == "__main__":
    unittest.main()

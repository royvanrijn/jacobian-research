#!/usr/bin/env python3
"""Pinned tests for the exact Fermigier arithmetic generic-rank theorem."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
PROGRAM = ROOT / "elliptic-curves"
sys.path.insert(0, str(CAS))
sys.path.insert(0, str(PROGRAM))

import verify_fermigier_generic_rank_exact as verifier  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_fermigier_generic_rank_exact.json"
)
EXPECTED_ARTIFACT_SHA256 = (
    "61bf11ae14db1aedcf7809697c96e78f2c2978c22f4e7b8fd894de74628b3de7"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierGenericRankExactTests(unittest.TestCase):
    def test_pinned_theorem_artifact(self) -> None:
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(
            data["result_sha256"],
            "0f9b048f7d56a68b7f3343d2bb49153bd5989da903040a7c75118e55c15d3a61",
        )
        self.assertEqual(
            data["provenance"]["script_sha256"],
            "74de70ce1c362a8e73bafbded149abad576e363fcbc97b086213bc6d5af719d5",
        )
        theorem = data["theorem"]
        self.assertEqual(
            theorem["arithmetic_generic_Mordell_Weil_rank_over_Q_of_u"], 12
        )
        self.assertEqual(
            theorem[
                "geometric_generic_Mordell_Weil_rank_interval_over_Qbar_of_u"
            ],
            [12, 13],
        )
        self.assertFalse(theorem["Tate_conjecture_assumed"])
        self.assertEqual(data["fiber_inventory"]["trivial_lattice"]["rank"], 5)
        self.assertEqual(
            data["exact_generic_sections"]["combined_exact_matrix_rank"], 12
        )
        self.assertEqual(
            data["H2_frobenius"]["full_positive_p_eigenvalue_multiplicity"],
            17,
        )
        self.assertEqual(
            data["H2_frobenius"][
                "residual_characteristic_polynomial_at_positive_p"
            ],
            2_136_275_316,
        )

    def test_independent_exact_replay(self) -> None:
        result = verifier.verify_generic_rank()
        self.assertEqual(
            result["model_bridge"]["literal_symmetric_shift"], "s=2*u"
        )
        self.assertEqual(
            result["fiber_inventory"]["infinity_fiber"]["Kodaira_type"],
            "I4",
        )
        self.assertTrue(
            result["fiber_inventory"]["infinity_fiber"]["split_over_Q"]
        )
        self.assertEqual(
            result["point_counts"]["1"]["surface_point_count"], 2_244
        )
        self.assertEqual(
            result["point_counts"]["2"]["surface_point_count"], 2_856_000
        )
        self.assertEqual(
            result["point_counts"]["1"]["h2_frobenius_trace"], 562
        )
        self.assertEqual(
            result["point_counts"]["2"]["h2_frobenius_trace"], 30_238
        )
        self.assertEqual(
            result["H2_frobenius"]["residual_traces"],
            {"degree_1": -135, "degree_2": 1_661},
        )
        self.assertEqual(
            result["theorem"][
                "arithmetic_generic_Mordell_Weil_rank_over_Q_of_u"
            ],
            12,
        )


if __name__ == "__main__":
    unittest.main()

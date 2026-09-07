#!/usr/bin/env python3
"""Focused replay checks for the exact rank-15 T=490/9 certificate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/certify_mestre_rank15_490_9.py"
FRONTIER_SCRIPT = (
    ROOT / "elliptic-curves/cas/search_mestre_rank14_pair_rational_frontier.py"
)
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_mestre_rank15_490_9.json"
FRONTIER = (
    ROOT
    / "artifacts/generated-results"
    / "elliptic_mestre_rank14_pair_rational_frontier.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "622f0d563f7b34d9a06635d79da992066b86b3797dff3f496c7b7959c8f7bd12"
)
EXPECTED_ARTIFACT_SHA256 = (
    "50b2b9c8bd24bcb5533534446af6404f3a9a761b5f33e0e28e04dc572227f950"
)
EXPECTED_FRONTIER_SHA256 = (
    "87e2d278cc1ee0653d1a4f871c1e34ed3d03babe1c1cd2ffe6712b7608efaee7"
)
EXPECTED_FRONTIER_SCRIPT_SHA256 = (
    "2f6251c67e2eb3cee2eca37d7e866913e9d5de73d30e3bfcb253641454d40d5f"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRank15CertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_inputs_and_provenance(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(sha256(FRONTIER), EXPECTED_FRONTIER_SHA256)
        self.assertEqual(sha256(FRONTIER_SCRIPT), EXPECTED_FRONTIER_SCRIPT_SHA256)
        provenance = self.data["provenance"]
        self.assertEqual(provenance["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(provenance["frontier_sha256"], EXPECTED_FRONTIER_SHA256)
        self.assertEqual(
            provenance["frontier_script_sha256"],
            EXPECTED_FRONTIER_SCRIPT_SHA256,
        )
        self.assertEqual(provenance["external_process_calls"], 0)
        self.assertEqual(
            self.data["result_sha256"],
            "31ff386e802c368efdcdd027168883e8edd459921eff431de71bfcfa5e517648",
        )

    def test_exact_curve_and_conductor(self) -> None:
        curve = self.data["curve"]
        self.assertEqual(curve["roots"], [0, 7, 121, 128, 183, 194])
        self.assertEqual(Fraction(curve["parameter"]), Fraction(490, 9))
        self.assertEqual(Fraction(curve["sign_equivalent_parameter"]), Fraction(-490, 9))
        self.assertEqual(
            curve["conductor"],
            "1468617013201344525723305189723172651230013579564220710",
        )
        self.assertLess(float(curve["log_conductor"]), 182.72)
        self.assertTrue(curve["below_strict_log_conductor_target_numerically"])
        a1, a2, a3, coefficient_a, coefficient_b = map(
            Fraction, curve["weierstrass_coefficients"]
        )
        self.assertEqual((a1, a2, a3), (0, 0, 0))
        for point in self.data["points"]:
            x_value, y_value = Fraction(point["x"]), Fraction(point["y"])
            self.assertEqual(
                y_value**2,
                x_value**3 + coefficient_a * x_value + coefficient_b,
            )

    def test_exact_finite_reduction_certificate(self) -> None:
        certificate = self.data["finite_reduction_certificate"]
        self.assertEqual(certificate["descent_modulus"], 3)
        self.assertEqual(certificate["point_count"], 15)
        self.assertEqual(certificate["combined_exact_rank_over_F3"], 15)
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 15)
        self.assertEqual(
            certificate["independent_subset_indices_one_based"], list(range(1, 16))
        )
        self.assertEqual(
            certificate["certificate_primes"],
            [37, 47, 59, 61, 67, 71, 83, 89, 97, 101, 113, 131, 149, 181, 191],
        )
        self.assertEqual(
            certificate["rational_3_torsion_exclusion"]["prime"], 19
        )
        self.assertEqual(
            certificate["point_sha256"],
            "fcbc9e86472490b3e5db4980d8d579c42657c41e3685be61953e1b95f5f9ed8e",
        )
        claim = self.data["claim"]
        self.assertEqual(claim["certified_algebraic_rank_lower_bound"], 15)
        self.assertTrue(claim["does_not_claim_exact_mordell_weil_rank"])
        self.assertTrue(claim["does_not_hit_rank21_or_rank30_target"])


if __name__ == "__main__":
    unittest.main()

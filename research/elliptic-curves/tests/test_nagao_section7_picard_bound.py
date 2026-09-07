from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results/elliptic-curves"
sys.path.insert(0, str(CAS))

from verify_nagao_section7_picard_bound import (  # noqa: E402
    EXPECTED_H2_TRACES,
    EXPECTED_SURFACE_POINT_COUNTS,
    reconstruct_residual_frobenius,
    surface_point_count,
    verify_good_reduction,
)


SCRIPT = CAS / "verify_nagao_section7_picard_bound.py"
LINEAR_VERIFIER = CAS / "verify_nagao_section7_linear_sections.py"
ARTIFACT = GENERATED / "elliptic_nagao_section7_picard_bound.json"


class NagaoSection7PicardBoundTests(unittest.TestCase):
    def test_good_reduction_and_exact_point_counts(self) -> None:
        reduction = verify_good_reduction()
        self.assertTrue(reduction["good_reduction_of_resolved_k3"])
        self.assertEqual(reduction["infinity_fiber"], "split I4")
        for degree in (1, 2):
            result = surface_point_count(29, degree)
            self.assertEqual(
                result["surface_point_count"],
                EXPECTED_SURFACE_POINT_COUNTS[degree],
            )
            self.assertEqual(
                result["h2_frobenius_trace"], EXPECTED_H2_TRACES[degree]
            )

    def test_residual_frobenius_has_only_one_cyclotomic_eigenvalue(self) -> None:
        result = reconstruct_residual_frobenius(29, 370, 16318)
        self.assertEqual(result["unique_real_eigenvalue"], "-29")
        self.assertEqual(
            result["residual_characteristic_polynomial"],
            "X**5 + 123*X**4 + 6554*X**3 + 190066*X**2 + 2999847*X + 20511149",
        )
        self.assertTrue(
            all(
                record["gcd_degree"] == 0
                for record in result["cyclotomic_gcd_checks"]
            )
        )
        self.assertEqual(result["residual_positive_p_eigenvalue_multiplicity"], 0)
        self.assertEqual(result["rational_Neron_Severi_rank"], 17)
        self.assertEqual(
            result["arithmetic_generic_Mordell_Weil_rank_over_Q_of_T"], 12
        )
        self.assertEqual(result["geometric_Picard_rank_upper_bound"], 18)
        self.assertEqual(
            result["geometric_generic_Mordell_Weil_rank_upper_bound"], 13
        )

    def test_generated_artifact_and_hashes(self) -> None:
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["proved_arithmetic_generic_rank_over_Q_of_T"], 12)
        self.assertEqual(data["proved_generic_rank_interval"], [12, 13])
        self.assertEqual(
            data["residual_frobenius"]["geometric_Picard_rank_upper_bound"],
            18,
        )
        self.assertFalse(data["target_hit"])
        self.assertEqual(
            data["script_sha256"], hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(
            data["rank12_certificate_verifier_sha256"],
            hashlib.sha256(LINEAR_VERIFIER.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()

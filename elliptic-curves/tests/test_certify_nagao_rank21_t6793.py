from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
sys.path.insert(0, str(CAS))

from certify_nagao_rank21_t6793 import (  # noqa: E402
    CERTIFIED_RANK_LOWER_BOUND,
    EXPECTED_CERTIFICATE_PRIMES,
    EXPECTED_CONDUCTOR,
    EXPECTED_INPUT_POINT_SHA256,
    EXPECTED_SATURATED_BASIS_SHA256,
    EXPECTED_SOURCE_SHA256,
    PARAMETER_T,
    load_pinned_input,
    sha256_file,
)
from nagao_1994 import RANK21_CONSTRUCTION, short_jacobian_coefficients  # noqa: E402
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve  # noqa: E402


SCRIPT = CAS / "certify_nagao_rank21_t6793.py"
SOURCE = GENERATED / "elliptic_nagao_rank21_unbiased.json"
ARTIFACT = GENERATED / "elliptic_nagao_rank21_t6793_rank19_certificate.json"


class NagaoRank21T6793CertificateTests(unittest.TestCase):
    def test_pinned_source_subset(self) -> None:
        self.assertEqual(str(PARAMETER_T), "6793/64")
        self.assertEqual(sha256_file(SOURCE), EXPECTED_SOURCE_SHA256)
        points, provenance = load_pinned_input(SOURCE)
        self.assertEqual(len(points), CERTIFIED_RANK_LOWER_BOUND)
        self.assertEqual(point_digest(points), EXPECTED_INPUT_POINT_SHA256)
        self.assertEqual(provenance["stable_numerical_rank"], 19)
        coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in points))

    def test_certificate_constants_are_pinned(self) -> None:
        self.assertEqual(CERTIFIED_RANK_LOWER_BOUND, 19)
        self.assertEqual(len(EXPECTED_CERTIFICATE_PRIMES), 16)
        self.assertEqual(EXPECTED_CERTIFICATE_PRIMES[0], 17)
        self.assertEqual(EXPECTED_CERTIFICATE_PRIMES[-1], 131)
        self.assertEqual(len(str(EXPECTED_CONDUCTOR)), 69)
        self.assertEqual(len(EXPECTED_SATURATED_BASIS_SHA256), 64)

    def test_generated_certificate_is_exact_and_pinned_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the dedicated rank-19 certificate has not been generated")
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["candidate"]["parameter_t"], "6793/64")
        self.assertEqual(data["candidate"]["conductor"], str(EXPECTED_CONDUCTOR))
        self.assertEqual(data["candidate"]["root_number"], -1)
        self.assertTrue(data["candidate"]["below_strict_log_conductor_target"])
        certificate = data["exact_rank_certificate"]
        self.assertEqual(
            certificate["certified_algebraic_rank_lower_bound"],
            CERTIFIED_RANK_LOWER_BOUND,
        )
        self.assertEqual(certificate["combined_exact_rank_over_F2"], 19)
        self.assertEqual(len(certificate["saturated_basis"]), 19)
        self.assertEqual(
            certificate["saturated_basis_sha256"],
            EXPECTED_SATURATED_BASIS_SHA256,
        )
        self.assertFalse(data["interpretation"]["target_reached"])


if __name__ == "__main__":
    unittest.main()

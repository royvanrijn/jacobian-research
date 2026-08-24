from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results/elliptic-curves"
sys.path.insert(0, str(CAS))

from certify_nagao_rank20_t5081 import (  # noqa: E402
    CERTIFIED_RANK_LOWER_BOUND,
    CONSTRUCTION,
    EXPECTED_CERTIFICATE_PRIMES,
    EXPECTED_CONDUCTOR,
    EXPECTED_POOL_SHA256,
    EXPECTED_SATURATED_BASIS_SHA256,
    EXPECTED_SELECTED_INDICES,
    PAPER_PARAMETER,
    PARAMETER_T,
    ROOTS,
    exact_curve_data,
    sha256_file,
)
from triage_nagao_rank13_finalists import point_on_short_curve  # noqa: E402


SCRIPT = CAS / "certify_nagao_rank20_t5081.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank20_t5081_rank20_certificate.json"


class NagaoRank20T5081CertificateTests(unittest.TestCase):
    def test_exact_constructor_convention_and_visible_points(self) -> None:
        self.assertEqual(ROOTS, (346, 260, 255, 146, 55, 0))
        self.assertEqual(str(PAPER_PARAMETER), "5081/94")
        self.assertEqual(str(PARAMETER_T), "5081/47")
        self.assertEqual(PARAMETER_T, 2 * PAPER_PARAMETER)
        self.assertEqual(CONSTRUCTION.quartic_square_scale, 7800)
        _, visible, coefficients = exact_curve_data()
        self.assertEqual(len(visible), 12)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in visible))

    def test_pinned_exact_certificate_constants(self) -> None:
        self.assertEqual(CERTIFIED_RANK_LOWER_BOUND, 20)
        self.assertEqual(len(EXPECTED_SELECTED_INDICES), 20)
        self.assertEqual(len(EXPECTED_CERTIFICATE_PRIMES), 15)
        self.assertEqual(EXPECTED_CERTIFICATE_PRIMES[-1], 173)
        self.assertEqual(len(EXPECTED_POOL_SHA256), 64)
        self.assertEqual(len(EXPECTED_SATURATED_BASIS_SHA256), 64)
        self.assertEqual(len(str(EXPECTED_CONDUCTOR)), 76)

    def test_generated_certificate_is_exact_and_pinned_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the rank-20 certificate has not been generated")
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["candidate"]["paper_parameter_t"], "5081/94")
        self.assertEqual(data["candidate"]["constructor_parameter_T"], "5081/47")
        self.assertEqual(data["candidate"]["conductor"], str(EXPECTED_CONDUCTOR))
        self.assertEqual(data["candidate"]["root_number"], 1)
        self.assertTrue(data["candidate"]["below_strict_log_conductor_target"])
        certificate = data["exact_rank_certificate"]
        self.assertEqual(certificate["combined_exact_rank_over_F2"], 20)
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 20)
        self.assertEqual(len(certificate["saturated_basis"]), 20)
        self.assertEqual(
            certificate["saturated_basis_sha256"], EXPECTED_SATURATED_BASIS_SHA256
        )
        self.assertFalse(data["interpretation"]["target_reached"])


if __name__ == "__main__":
    unittest.main()

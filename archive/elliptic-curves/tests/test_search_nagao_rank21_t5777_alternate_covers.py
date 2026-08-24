from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from search_nagao_rank21_t5777_alternate_covers import (  # noqa: E402
    CERTIFICATE_SHA256,
    load_exact_basis,
    sha256_file,
)
from search_nagao_u135_alternate_covers import (  # noqa: E402
    full_coset_identity_frontier,
)


CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank21_t5777_rank17_certificate.json"
)


class NagaoRank21T5777AlternateCoverTests(unittest.TestCase):
    def test_pinned_exact_basis_loads(self) -> None:
        coefficients, basis, data = load_exact_basis(CERTIFICATE)
        self.assertEqual(len(coefficients), 5)
        self.assertEqual(len(basis), 17)
        self.assertEqual(
            data["exact_rank_certificate"]["certified_algebraic_rank_lower_bound"],
            17,
        )
        self.assertEqual(sha256_file(CERTIFICATE), CERTIFICATE_SHA256)

    def test_small_full_coset_scan_has_all_masks(self) -> None:
        coefficients, basis, _ = load_exact_basis(CERTIFICATE)
        frontier = full_coset_identity_frontier(
            coefficients, basis[:3], retain_count=7
        )
        self.assertEqual({score[2] for score, _ in frontier}, set(range(1, 8)))


if __name__ == "__main__":
    unittest.main()

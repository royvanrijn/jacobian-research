from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from alternate_quartic_covers import alternate_cover, short_subset_sum  # noqa: E402
from search_nagao_rank21_t6793_alternate_covers import (  # noqa: E402
    CERTIFICATE_SHA256,
    INPUT_CERTIFIED_RANK,
    PARAMETER_T,
    load_exact_basis,
    sha256_file,
)
from search_nagao_u135_alternate_covers import (  # noqa: E402
    full_coset_identity_frontier,
)


CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank21_t6793_rank19_certificate.json"
)


class NagaoRank21T6793AlternateCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.coefficients, cls.basis, cls.certificate, cls.conductor = (
            load_exact_basis(CERTIFICATE)
        )

    def test_pinned_exact_basis_loads(self) -> None:
        self.assertEqual(PARAMETER_T.numerator, 6793)
        self.assertEqual(PARAMETER_T.denominator, 64)
        self.assertEqual(len(self.coefficients), 5)
        self.assertEqual(len(self.basis), INPUT_CERTIFIED_RANK)
        self.assertEqual(
            self.certificate["certified_algebraic_rank_lower_bound"], 19
        )
        self.assertEqual(self.conductor["root_number"], -1)
        self.assertEqual(sha256_file(CERTIFICATE), CERTIFICATE_SHA256)

    def test_small_full_coset_scan_has_every_nonzero_mask(self) -> None:
        frontier = full_coset_identity_frontier(
            self.coefficients, self.basis[:4], retain_count=15
        )
        self.assertEqual({score[2] for score, _ in frontier}, set(range(1, 16)))

    def test_alternate_cover_round_trip_is_exact(self) -> None:
        base_point = short_subset_sum(
            self.coefficients, self.basis, (0, 1, 2)
        )
        self.assertIsNotNone(base_point)
        cover = alternate_cover(self.coefficients, base_point)
        for point in self.basis[:6]:
            if point == base_point:
                continue
            cover_point = cover.curve_point_to_cover(point)
            self.assertEqual(cover.cover_point_to_curve(cover_point), point)


if __name__ == "__main__":
    unittest.main()

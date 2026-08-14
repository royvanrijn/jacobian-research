from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from search_nagao_u135_alternate_covers import (  # noqa: E402
    full_coset_identity_frontier,
)
from search_nagao_u75_alternate_covers import (  # noqa: E402
    EXPECTED_HEIGHT_SUBSET,
    PARAMETER_T,
    PARAMETER_U,
    SIEVE_ROWS,
    exact_seed_pool,
)
from nagao_1994 import rank13_base_parameter  # noqa: E402
from triage_nagao_rank13_finalists import point_on_short_curve  # noqa: E402


Q = Fraction


class U75AlternateCoverSearchTests(unittest.TestCase):
    def test_pinned_base_change_and_sieve_rows(self) -> None:
        self.assertEqual(PARAMETER_U, Q(75, 2))
        self.assertEqual(PARAMETER_T, Q(1181, 4))
        self.assertEqual(rank13_base_parameter(PARAMETER_U), PARAMETER_T)
        self.assertEqual(len(SIEVE_ROWS), 4)
        self.assertEqual(len(EXPECTED_HEIGHT_SUBSET), 15)

    def test_exact_eighteen_point_seed_pool(self) -> None:
        coefficients, pool, records = exact_seed_pool()
        self.assertEqual(len(pool), 18)
        self.assertEqual(len(records), 4)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in pool))
        self.assertTrue(
            all(record["exact_quartic_and_jacobian_membership_checked"] for record in records)
        )

    def test_small_full_coset_frontier_has_every_nonzero_mask(self) -> None:
        coefficients, pool, _ = exact_seed_pool()
        basis = tuple(pool[index - 1] for index in EXPECTED_HEIGHT_SUBSET[:3])
        frontier = full_coset_identity_frontier(
            coefficients, basis, retain_count=7
        )
        self.assertEqual({score[2] for score, _ in frontier}, set(range(1, 8)))


if __name__ == "__main__":
    unittest.main()

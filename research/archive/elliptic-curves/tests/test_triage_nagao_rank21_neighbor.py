from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import shutil
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from nagao_1994 import (  # noqa: E402
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_POINTS,
    point_on_extended_weierstrass,
)
from triage_nagao_rank21_neighbor import (  # noqa: E402
    CANDIDATE_PARAMETER,
    bounded_quartic_points,
    exact_visible_seeds,
    load_candidate,
    map_and_deduplicate,
)


ROOT = Path(__file__).resolve().parents[2]
SOURCE_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_neighborhood.json"
)


class NagaoRank21NeighborTriageTests(unittest.TestCase):
    def test_source_candidate_is_pinned_without_a_rank_claim(self) -> None:
        candidate = load_candidate(SOURCE_ARTIFACT)
        self.assertEqual(candidate.parameter, Q(6041, 198))
        self.assertEqual(candidate.parameter, CANDIDATE_PARAMETER)
        self.assertTrue(candidate.below_target)
        self.assertEqual(candidate.root_number, 1)
        self.assertTrue(candidate.log_conductor.startswith("170.765123121845207336"))
        self.assertEqual(
            [check["observed_discriminant_valuation"] for check in candidate.source_local_checks],
            [4, 5, 4, 3, 3],
        )

    def test_candidate_and_published_visible_seeds_map_exactly(self) -> None:
        for parameter in (CANDIDATE_PARAMETER, RANK21_CONSTRUCTOR_PARAMETER):
            quartic, jacobian, coefficients = exact_visible_seeds(parameter)
            self.assertEqual(len(quartic), 12)
            self.assertEqual(len({point[0] for point in quartic}), 12)
            self.assertEqual(len(jacobian), 12)
            self.assertEqual(len({point[0] for point in jacobian}), 12)
            self.assertEqual(len(coefficients), 5)

    def test_sign_pair_deduplication_recognizes_a_visible_seed(self) -> None:
        quartic, jacobian, _ = exact_visible_seeds(CANDIDATE_PARAMETER)
        x_value, y_value = quartic[0]
        records, new_images, zero_ordinates = map_and_deduplicate(
            CANDIDATE_PARAMETER,
            ((x_value, y_value), (x_value, -y_value)),
            quartic,
            jacobian,
        )
        self.assertEqual(len(records), 1)
        self.assertTrue(records[0]["visible_section_abscissa"])
        self.assertTrue(records[0]["duplicate_seed_or_prior_sign_pair"])
        self.assertEqual(new_images, ())
        self.assertEqual(zero_ordinates, ())

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is optional")
    def test_height_50000_search_returns_only_visible_sign_pairs(self) -> None:
        quartic, jacobian, _ = exact_visible_seeds(CANDIDATE_PARAMETER)
        raw, _, _ = bounded_quartic_points(
            CANDIDATE_PARAMETER,
            height_bound=50_000,
            timeout=10,
            stack_bytes=128_000_000,
        )
        records, new_images, zero_ordinates = map_and_deduplicate(
            CANDIDATE_PARAMETER, raw, quartic, jacobian
        )
        self.assertEqual(len(raw), 12)
        self.assertEqual(len(records), 6)
        self.assertTrue(all(record["visible_section_abscissa"] for record in records))
        self.assertEqual(new_images, ())
        self.assertEqual(zero_ordinates, ())

    def test_all_published_points_are_exact_calibration_points(self) -> None:
        self.assertEqual(len(RANK21_PUBLISHED_POINTS), 21)
        self.assertTrue(
            all(
                point_on_extended_weierstrass(RANK21_PUBLISHED_MODEL, point)
                for point in RANK21_PUBLISHED_POINTS
            )
        )


if __name__ == "__main__":
    unittest.main()

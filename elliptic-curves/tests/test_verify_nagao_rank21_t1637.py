from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    RANK21_CONSTRUCTION,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from triage_nagao_rank13_finalists import point_on_short_curve  # noqa: E402
from verify_nagao_rank21_t1637 import (  # noqa: E402
    EXPECTED_HEIGHT_SUBSET,
    EXPECTED_NEW_IMAGES,
    EXPECTED_POOL_SIZE,
    EXPECTED_SIGNED_POINTS,
    PARAMETER_T,
    sha256_file,
)


SCRIPT = CAS / "verify_nagao_rank21_t1637.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank21_t1637_rank16_certificate.json"


class NagaoRank21T1637VerifierTests(unittest.TestCase):
    def test_pinned_search_counts_and_subset(self) -> None:
        self.assertEqual(str(PARAMETER_T), "1637/12")
        self.assertEqual(EXPECTED_SIGNED_POINTS, 252)
        self.assertEqual(EXPECTED_NEW_IMAGES, 114)
        self.assertEqual(EXPECTED_POOL_SIZE, 126)
        self.assertEqual(len(EXPECTED_HEIGHT_SUBSET), 16)
        self.assertEqual(EXPECTED_HEIGHT_SUBSET, tuple(range(1, 12)) + (13, 14, 15, 16, 17))

    def test_visible_points_map_exactly(self) -> None:
        quartic = primitive_quartic_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
        short = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
        points = primitive_visible_points(RANK21_CONSTRUCTION, PARAMETER_T)
        self.assertEqual(len(points), 12)
        for point in points:
            self.assertEqual(point[1] ** 2, quartic_value(quartic, point[0]))
            image = quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION, PARAMETER_T, point
            )
            self.assertTrue(point_on_short_curve(short, image))

    def test_generated_certificate_is_pinned_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the exact rank certificate has not been generated")
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["candidate"]["parameter_t"], "1637/12")
        self.assertEqual(
            data["exact_rank_certificate"]["certified_algebraic_rank_lower_bound"],
            16,
        )
        self.assertEqual(data["candidate"]["root_number"], 1)


if __name__ == "__main__":
    unittest.main()

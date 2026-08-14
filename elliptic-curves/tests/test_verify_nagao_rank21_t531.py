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
from verify_nagao_rank21_t531 import (  # noqa: E402
    EXPECTED_CONDUCTOR,
    EXPECTED_HEIGHT_SUBSET,
    EXPECTED_NEW_IMAGES,
    EXPECTED_POOL_SIZE,
    EXPECTED_SIGNED_POINTS,
    PARAMETER_T,
    sha256_file,
)


SCRIPT = CAS / "verify_nagao_rank21_t531.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank21_t531_rank17_certificate.json"


class NagaoRank21T531VerifierTests(unittest.TestCase):
    def test_pinned_search_counts_subset_and_conductor(self) -> None:
        self.assertEqual(str(PARAMETER_T), "531/2")
        self.assertEqual(EXPECTED_SIGNED_POINTS, 228)
        self.assertEqual(EXPECTED_NEW_IMAGES, 102)
        self.assertEqual(EXPECTED_POOL_SIZE, 114)
        self.assertEqual(len(EXPECTED_HEIGHT_SUBSET), 17)
        self.assertEqual(
            EXPECTED_HEIGHT_SUBSET,
            tuple(range(1, 12)) + (13, 14, 15, 16, 17, 19),
        )
        self.assertEqual(
            EXPECTED_CONDUCTOR,
            17363106716312881727260568367360064041745970511465978238830,
        )

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
        self.assertEqual(data["candidate"]["parameter_t"], "531/2")
        self.assertEqual(data["candidate"]["conductor"], str(EXPECTED_CONDUCTOR))
        self.assertEqual(data["candidate"]["root_number"], -1)
        self.assertEqual(
            data["exact_rank_certificate"]["certified_algebraic_rank_lower_bound"],
            17,
        )


if __name__ == "__main__":
    unittest.main()

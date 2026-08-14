from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.kihara import (  # noqa: E402
    fifteen_visible_points,
    kihara_quartic,
    kihara_rank14_replay,
    verify_kihara_rank14_manifest,
)
from ecsearch.rank_certification import is_on_weierstrass_curve  # noqa: E402


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "kihara_rank14_t2_v1.json"
)


class KiharaRank14Tests(unittest.TestCase):
    def test_t2_quartic_and_visible_points(self) -> None:
        model = kihara_quartic(2)
        self.assertEqual((model.p, model.q, model.u), (80, -288, 18050187264))
        self.assertEqual(len(model.roots), 12)
        self.assertEqual(len(model.product), 13)
        self.assertEqual(len(model.square_part), 7)
        self.assertEqual(len(model.quartic), 5)
        points = fifteen_visible_points(model)
        self.assertEqual(len(points), 15)
        self.assertEqual(points[12][0], Fraction(29183481217024, 421))
        self.assertEqual(points[13][0], Fraction(4755487195136, 269))
        self.assertEqual(points[14][0], -517996544)

    def test_weierstrass_map(self) -> None:
        replay = kihara_rank14_replay()
        self.assertEqual(len(replay.weierstrass_points), 14)
        self.assertTrue(
            all(
                is_on_weierstrass_curve(
                    replay.weierstrass_coefficients,
                    point,
                )
                for point in replay.weierstrass_points
            )
        )

    def test_pinned_exact_independence_certificate(self) -> None:
        manifest = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        verify_kihara_rank14_manifest(manifest)
        certificate = manifest["independence_certificate"]
        self.assertEqual(certificate["relation_prime"], 5)
        self.assertEqual(len(certificate["rows"]), 14)


if __name__ == "__main__":
    unittest.main()

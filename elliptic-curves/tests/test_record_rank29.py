from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.rank_certification import is_on_weierstrass_curve  # noqa: E402
from ecsearch.record_rank29 import (  # noqa: E402
    load_rank29_baseline,
    verify_rank29_manifest,
)


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "elkies_klagsbrun_e29_independence_v1.json"
)


class RecordRank29Tests(unittest.TestCase):
    def test_all_published_points_are_on_curve(self) -> None:
        coefficients, points = load_rank29_baseline()
        self.assertEqual(len(points), 29)
        self.assertTrue(
            all(is_on_weierstrass_curve(coefficients, point) for point in points)
        )

    def test_pinned_exact_independence_certificate(self) -> None:
        manifest = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        verify_rank29_manifest(manifest)
        certificate = manifest["independence_certificate"]
        self.assertEqual(certificate["relation_prime"], 2)
        self.assertEqual(len(certificate["rows"]), 29)
        self.assertIn("no exact-rank claim", manifest["claim_status"])


if __name__ == "__main__":
    unittest.main()

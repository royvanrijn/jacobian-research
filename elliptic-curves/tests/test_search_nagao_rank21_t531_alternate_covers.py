from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
sys.path.insert(0, str(CAS))

from search_nagao_rank21_t531_alternate_covers import (  # noqa: E402
    CERTIFICATE_SHA256,
    PARAMETER_T,
    load_exact_basis,
    sha256_file,
)
from search_nagao_u135_alternate_covers import (  # noqa: E402
    full_coset_identity_frontier,
)


CERTIFICATE = GENERATED / "elliptic_nagao_rank21_t531_rank17_certificate.json"
ARTIFACT = GENERATED / "elliptic_nagao_rank21_t531_alternate_covers.json"
SCRIPT = CAS / "search_nagao_rank21_t531_alternate_covers.py"


class NagaoRank21T531AlternateCoverTests(unittest.TestCase):
    def test_pinned_exact_basis_loads(self) -> None:
        coefficients, basis, data = load_exact_basis(CERTIFICATE)
        self.assertEqual(str(PARAMETER_T), "531/2")
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

    def test_generated_artifact_is_pinned_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the bounded alternate-cover artifact has not been run")
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["candidate"]["parameter_t"], "531/2")
        self.assertEqual(
            data["declared_budget"][
                "all_nonzero_certified_mod2_classes_identity_scored"
            ],
            131071,
        )
        self.assertEqual(data["declared_budget"]["pilot_chart_count"], 60)
        self.assertEqual(data["declared_budget"]["escalation_chart_count"], 8)
        self.assertLessEqual(data["declared_budget"]["deep_chart_count"], 2)
        self.assertGreaterEqual(
            data["results"]["certified_rank_lower_bound_after_search"], 17
        )


if __name__ == "__main__":
    unittest.main()

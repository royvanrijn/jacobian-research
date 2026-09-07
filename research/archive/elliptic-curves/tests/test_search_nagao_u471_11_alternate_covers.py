from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
sys.path.insert(0, str(CAS))

from search_nagao_u135_alternate_covers import (  # noqa: E402
    full_coset_identity_frontier,
)
from search_nagao_u471_11_alternate_covers import (  # noqa: E402
    CERTIFICATE_SHA256,
    PARAMETER_T,
    PARAMETER_U,
    RANK_GAIN_SHA256,
    load_exact_inputs,
    sha256_file,
)


Q = Fraction
CERTIFICATE = GENERATED / "elliptic_nagao_rank17_frontier_certificate.json"
RANK_GAIN = GENERATED / "elliptic_nagao_rank13_rank_gain_search.json"
ARTIFACT = GENERATED / "elliptic_nagao_u471_11_alternate_covers.json"
SCRIPT = CAS / "search_nagao_u471_11_alternate_covers.py"


class U471AlternateCoverSearchTests(unittest.TestCase):
    def test_pinned_exact_certificate_and_subset_lineage(self) -> None:
        target, checkpoint, selected_digest = load_exact_inputs(
            CERTIFICATE, RANK_GAIN
        )
        self.assertEqual(PARAMETER_U, Q(471, 11))
        self.assertEqual(PARAMETER_T, Q(5579, 22))
        self.assertEqual(target.parameter_t, PARAMETER_T)
        self.assertEqual(target.certified_rank_lower_bound, 17)
        self.assertEqual(len(target.saturated_basis), 17)
        self.assertEqual(checkpoint.stable_pool_numerical_rank, 17)
        self.assertEqual(
            selected_digest,
            "81d4881b4da5ea3b35d9e8b79419acc9a7193ffe6365fc58abba9268f5459c11",
        )
        self.assertEqual(sha256_file(CERTIFICATE), CERTIFICATE_SHA256)
        self.assertEqual(sha256_file(RANK_GAIN), RANK_GAIN_SHA256)

    def test_small_full_coset_scan_has_all_masks(self) -> None:
        target, _, _ = load_exact_inputs(CERTIFICATE, RANK_GAIN)
        frontier = full_coset_identity_frontier(
            target.jacobian_coefficients,
            target.saturated_basis[:3],
            retain_count=7,
        )
        self.assertEqual({score[2] for score, _ in frontier}, set(range(1, 8)))

    def test_generated_artifact_is_pinned_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the bounded alternate-cover artifact has not been run")
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["candidate"]["parameter_u"], "471/11")
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

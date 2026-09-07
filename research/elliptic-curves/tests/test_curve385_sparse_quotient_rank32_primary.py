#!/usr/bin/env python3
"""Regression tests for the promoted curve-385 primary sparse campaign."""

from __future__ import annotations

import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
MANIFEST = ART / "curve385_sparse_quotient_rank32_primary_v1.json"
LEDGER_GZ = ART / "curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz"
PROMOTER = ROOT / "elliptic-curves/cas/promote_curve385_sparse_quotient_rank32_primary.py"

MANIFEST_SHA256 = "ff1b0a2e8dd29b9a34a3b81cf8db3bed5350d12ee0dcfd9dd3936f226d255d61"
LEDGER_GZIP_SHA256 = "08a2e416255910f733ef98283332e3a60a947350646329e4e2045cbc08d802c0"
LEDGER_RAW_SHA256 = "17600ab552c8c4c5184d8ec02c6743c475424998e2672a8eeacc3ee75df5b77d"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class Curve385SparseQuotientRank32PrimaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST.read_text())

    def test_frozen_artifact_hashes(self) -> None:
        self.assertEqual(digest(MANIFEST), MANIFEST_SHA256)
        self.assertEqual(digest(LEDGER_GZ), LEDGER_GZIP_SHA256)
        self.assertEqual(
            sha256(gzip.decompress(LEDGER_GZ.read_bytes())).hexdigest(),
            LEDGER_RAW_SHA256,
        )
        self.assertEqual(self.manifest["program"]["sha256"], digest(PROMOTER))

    def test_primary_stage_accounting(self) -> None:
        campaign = self.manifest["campaign"]
        self.assertEqual(campaign["starting_certified_rank_lower_bound"], 29)
        self.assertEqual(campaign["final_certified_rank_lower_bound"], 29)
        self.assertEqual(campaign["target_rank_lower_bound"], 32)
        self.assertFalse(campaign["exact_group_growth"])
        self.assertFalse(campaign["target_reached"])
        self.assertEqual(campaign["total_planned_chart_count"], 3_354)
        self.assertEqual(campaign["total_fresh_search_count"], 3_116)
        self.assertEqual(campaign["total_previously_searched_exact_chart_count"], 238)
        self.assertEqual(3_116 + 238, 3_354)

        stages = campaign["stage_summaries"]
        self.assertEqual([row["id"] for row in stages], ["natural-weight-1", "natural-weight-2"])
        self.assertEqual([row["planned_chart_count"] for row in stages], [516, 2_838])
        self.assertEqual([row["fresh_search_count"] for row in stages], [394, 2_722])
        self.assertEqual(
            [row["previously_searched_exact_chart_count"] for row in stages],
            [122, 116],
        )
        self.assertTrue(all(row["timeout_count"] == 0 for row in stages))
        self.assertTrue(all(row["pari_failure_count"] == 0 for row in stages))
        self.assertTrue(all(row["new_independent_direction_count"] == 0 for row in stages))
        self.assertTrue(all(row["finite_index_saturation_event_count"] == 0 for row in stages))

    def test_claim_boundary_and_next_stage(self) -> None:
        self.assertEqual(
            self.manifest["status"],
            "PASS_COMPLETE_PRIMARY_SPARSE_CAMPAIGN_BOUNDED_NO_GROWTH",
        )
        boundary = " ".join(self.manifest["claim_boundary"])
        self.assertIn("not a rank upper bound", boundary)
        self.assertIn("Rank at least 32 remains open", boundary)
        self.assertEqual(
            self.manifest["next_precommitted_stage"],
            {
                "id": "alternate-a-weight-at-most-2",
                "index": 3,
                "new_chart_count": 3_268,
                "requires_explicit_stage_limit_escalation": True,
            },
        )

    def test_full_auditor_check(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(PROMOTER), "--check"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("C385SPARSEPRIMARY|status=PASS", completed.stdout)


if __name__ == "__main__":
    unittest.main()

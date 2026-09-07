"""Regression for the preserved curve-302 17-to-31 recovery chain."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "elliptic-curves/cas/audit_curve302_residual_strict_bootstrap_chain.py"
RESULT = ROOT / "artifacts/generated-results/elliptic-curves/curve302_residual_strict_bootstrap_chain_v1.json"


class Curve302ResidualStrictBootstrapChainTest(unittest.TestCase):
    def test_replays(self) -> None:
        completed = subprocess.run([sys.executable, str(AUDIT), "--check"], cwd=ROOT, check=True, capture_output=True, text=True)
        self.assertIn("rank_at_least=31|status=PASS", completed.stdout)
        payload = json.loads(RESULT.read_text())
        self.assertEqual(payload["status"], "PASS_EXACT_CHAIN_PRESERVED_RANK_AT_LEAST_31")
        self.assertEqual([(row["rank_before"], row["rank_after"]) for row in payload["recovery_arms"]], [
            (24, 25), (25, 26), (26, 27), (27, 28), (28, 29), (29, 30), (30, 31),
        ])
        self.assertEqual([row["certificate"]["mod2_rank_lower_bound"] for row in payload["recovery_arms"]], list(range(25, 32)))

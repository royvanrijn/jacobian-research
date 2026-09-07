from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from pointed_regression_sources import historical_digest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/curve385_sparse_restart_budget_v2.json"
)
BUILDER = CAS / "build_curve385_sparse_restart_budget.py"
sys.path.insert(0, str(CAS))

from curve385_sparse_restart_policy import (  # noqa: E402
    RANK_CHANGING,
    SATURATION_ONLY,
    classify_group_change,
    simulate_unit_rank_path,
    validate_accounting,
)


class Curve385SparseRestartBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        spec = importlib.util.spec_from_file_location("restart_builder", BUILDER)
        assert spec is not None and spec.loader is not None
        cls.builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.builder)
        cls.artifact = json.loads(ARTIFACT.read_text())

    def test_artifact_replays(self) -> None:
        with patch.object(self.builder, "digest", side_effect=historical_digest):
            self.assertEqual(self.builder.build(), self.artifact)
        self.assertEqual(
            self.artifact["status"],
            "FROZEN_INDEPENDENT_RESTART_BUDGETS_BEFORE_FUTURE_SEARCH",
        )

    def test_rank_and_saturation_changes_are_disjoint(self) -> None:
        rank_change = classify_group_change(
            rank_before=29,
            rank_after=30,
            basis_before_sha256="a",
            basis_after_sha256="b",
            finite_index_saturation_event_count=1,
        )
        saturation = classify_group_change(
            rank_before=29,
            rank_after=29,
            basis_before_sha256="a",
            basis_after_sha256="b",
            finite_index_saturation_event_count=1,
        )
        self.assertEqual(rank_change, RANK_CHANGING)
        self.assertEqual(saturation, SATURATION_ONLY)

    def test_two_saturations_do_not_consume_rank_budget(self) -> None:
        result = simulate_unit_rank_path(
            [
                SATURATION_ONLY,
                SATURATION_ONLY,
                RANK_CHANGING,
                RANK_CHANGING,
                RANK_CHANGING,
            ]
        )
        self.assertEqual(result["status"], "TARGET_REACHED")
        self.assertEqual(result["rank"], 32)
        self.assertEqual(
            result["restart_accounting"],
            {
                "rank_changing_group_change_count": 3,
                "saturation_only_group_change_count": 2,
            },
        )

    def test_same_rank_basis_change_requires_saturation_evidence(self) -> None:
        with self.assertRaisesRegex(ValueError, "finite-index saturation"):
            classify_group_change(
                rank_before=29,
                rank_after=29,
                basis_before_sha256="a",
                basis_after_sha256="b",
                finite_index_saturation_event_count=0,
            )

    def test_saturation_budget_is_independent_and_finite(self) -> None:
        result = simulate_unit_rank_path([SATURATION_ONLY] * 5)
        self.assertEqual(result["status"], "RESTART_BUDGET_EXCEEDED")
        self.assertEqual(result["exceeded"], SATURATION_ONLY)
        self.assertEqual(result["rank"], 29)
        self.assertEqual(
            result["restart_accounting"]["rank_changing_group_change_count"], 0
        )

    def test_checkpoint_accounting_requires_exact_nonnegative_v2_counters(self) -> None:
        with self.assertRaisesRegex(ValueError, "exactly the v2 counters"):
            validate_accounting({"combined_lattice_state_count": 4})
        with self.assertRaisesRegex(ValueError, "nonnegative integers"):
            validate_accounting(
                {
                    "rank_changing_group_change_count": 0,
                    "saturation_only_group_change_count": -1,
                }
            )


if __name__ == "__main__":
    unittest.main()

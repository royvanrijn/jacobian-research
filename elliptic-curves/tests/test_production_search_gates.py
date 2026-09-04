from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "cas/production_search_gates.py"
SPEC = importlib.util.spec_from_file_location("production_search_gates", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
GATES = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = GATES
SPEC.loader.exec_module(GATES)


class ProductionSearchGatesTests(unittest.TestCase):
    def test_incomplete_descent_schedules_but_does_not_exclude(self) -> None:
        record = GATES.production_gate_record(
            target_rank=32,
            search_limits={"chart_count": 344, "height": 100_000},
            scheduling_information=[
                {"kind": "incomplete_descent", "status": "BNF_PENDING"}
            ],
        )
        self.assertEqual(record["proof_gate"]["status"], GATES.OPEN_STATUS)
        self.assertTrue(record["search_budget_gate"]["bounded_search_authorized"])
        self.assertFalse(
            record["scheduling_information"][0]["used_as_mathematical_exclusion"]
        )

    def test_unconditional_upper_bound_below_target_is_the_only_exclusion(self) -> None:
        record = GATES.production_gate_record(
            target_rank=32,
            search_limits={"wall_seconds": 60},
            certified_rank_upper_bound=31,
            upper_bound_kind="complete unconditional 2-Selmer upper bound",
            upper_bound_evidence="certificate.json",
        )
        self.assertEqual(record["proof_gate"]["status"], GATES.EXCLUDED_STATUS)
        self.assertFalse(record["search_budget_gate"]["bounded_search_authorized"])

    def test_upper_bound_at_target_does_not_exclude_lower_bound_search(self) -> None:
        record = GATES.production_gate_record(
            target_rank=32,
            search_limits={"wall_seconds": 60},
            certified_rank_upper_bound=32,
            upper_bound_kind="complete unconditional 2-Selmer upper bound",
            upper_bound_evidence="certificate.json",
        )
        self.assertEqual(record["proof_gate"]["status"], GATES.OPEN_STATUS)
        self.assertTrue(record["search_budget_gate"]["bounded_search_authorized"])

    def test_certified_points_meet_target_without_descent(self) -> None:
        record = GATES.certified_point_lower_bound_record(
            certified_independent_rank=32,
            target_rank=32,
            curve_equations_verified=True,
            independence_evidence="finite-reduction-rank32.json",
        )
        self.assertEqual(record["status"], GATES.LOWER_BOUND_SUCCESS_STATUS)
        self.assertTrue(record["target_lower_bound_met"])
        self.assertFalse(record["descent_completion_required"])
        self.assertFalse(record["exact_rank_claimed"])

    def test_search_budget_must_be_explicit_and_finite(self) -> None:
        for limit in (float('inf'), float('-inf'), float('nan')):
            with self.subTest(limit=limit):
                with self.assertRaises(GATES.ProductionSearchGateError):
                    GATES.production_gate_record(
                        target_rank=32, search_limits={"wall_seconds": limit}
                    )
        with self.assertRaisesRegex(GATES.ProductionSearchGateError, "explicit limits"):
            GATES.production_gate_record(target_rank=32, search_limits={})
        with self.assertRaisesRegex(GATES.ProductionSearchGateError, "positive"):
            GATES.production_gate_record(
                target_rank=32, search_limits={"wall_seconds": 0}
            )


if __name__ == "__main__":
    unittest.main()

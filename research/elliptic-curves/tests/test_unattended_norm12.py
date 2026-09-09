import os
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from run_unattended_norm12 import next_action, process_token


class UnattendedGates(unittest.TestCase):
    def terminal(self, reason, rank=28):
        return dict(stop_reason=reason, rank_lower_bound=rank)

    def test_only_verified_budget_semantics_continue(self):
        self.assertEqual(next_action(28, self.terminal('CHART_BUDGET_EXHAUSTED')), 'CONTINUE')
        self.assertEqual(next_action(28, self.terminal('FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN')), 'EXHAUSTED')

    def test_gain_stops_and_cloud_reconciliation_has_priority(self):
        self.assertEqual(next_action(28, self.terminal('CHART_BUDGET_EXHAUSTED', 29)), 'STOP_CERTIFIED_GAIN')
        self.assertEqual(next_action(28, self.terminal('ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION', 29)), 'RECONCILE')

    def test_unknown_reason_and_decreasing_rank_fail_closed(self):
        with self.assertRaises(ArithmeticError):
            next_action(28, self.terminal('unexpected'))
        with self.assertRaises(ArithmeticError):
            next_action(28, self.terminal('CHART_BUDGET_EXHAUSTED', 27))

    def test_process_identity_does_not_treat_missing_pid_as_live(self):
        self.assertIsNotNone(process_token(os.getpid()))
        self.assertIsNone(process_token(-1))


if __name__ == '__main__':
    unittest.main()

import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from audit_strata60_mw16_accounting import same,totals


class IndependentAccountingTests(unittest.TestCase):
    def test_unresolved_curve_is_allocated_but_not_zero_filled(self):
        rows=[dict(gain=2,boxes=43,discovery=100,verification=10,worker_status='COMPLETE_DECLARED_POINT_ATTEMPT'),
              dict(gain=None,boxes=None,discovery=600,verification=300,worker_status='POINT_FAILED_OR_CENSORED')]
        r=totals(rows)
        self.assertEqual(r['allocated_boxes'],86)
        self.assertEqual(r['certified_gain_sum_known'],2)
        self.assertEqual(r['including_verification_seconds'],1010)
        self.assertEqual(r['unresolved_gain_curves'],1)
        self.assertIsNone(r['gain_per_discovery_second'])
        self.assertIsNone(r['completion_fraction'])

    def test_partial_certified_gain_and_timeout_cost_survive(self):
        rows=[dict(gain=1,boxes=4,discovery=600,verification=50,worker_status='POINT_FAILED_OR_CENSORED')]
        r=totals(rows)
        self.assertEqual(r['gain_per_discovery_second'],1/600)
        self.assertEqual(r['gain_per_completed_box'],1/4)
        self.assertEqual(r['completion_fraction'],4/43)

    def test_corrupt_counts_and_nonfinite_cost_are_rejected(self):
        for actual,expected in [(42,43),(0,None),(float('nan'),1.0),(float('inf'),1.0),(1.000001,1.0)]:
            with self.subTest(actual=actual):
                with self.assertRaises(ArithmeticError):same(actual,expected)
        same(1.0+1e-14,1.0)


if __name__=='__main__':unittest.main()

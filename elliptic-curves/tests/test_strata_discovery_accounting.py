import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from report_strata60_mw16_experiment import summarize,comparison


def population():
    return [dict(family=f,arm=a,certified_gain=g,completed_boxes=43,
                 discovery_worker_seconds=100,verification_seconds=20,
                 worker_status='COMPLETE_DECLARED_POINT_ATTEMPT')
            for f in range(5) for _ in range(4)
            for a,g in [('top',1),('moderate',2),('lower',3)]]


class DiscoveryAccountingTests(unittest.TestCase):
    def test_missing_certificate_is_not_zero_gain(self):
        rows=population();rows[0]['certified_gain']=None
        data=summarize(rows)
        self.assertEqual(data['unresolved_gain_curves'],1)
        self.assertIsNone(data['gain_per_discovery_second'])
        self.assertIsNone(data['gain_per_completed_box'])
        self.assertEqual(comparison(rows)['recommendation'],'INCONCLUSIVE')

    def test_partial_exposure_counts_completed_boxes_and_all_cost(self):
        row=population()[0];row.update(completed_boxes=3,certified_gain=2,
            discovery_worker_seconds=600,worker_status='POINT_FAILED_OR_CENSORED')
        data=summarize([row])
        self.assertEqual(data['completion_fraction'],3/43)
        self.assertEqual(data['gain_per_discovery_second'],2/600)
        self.assertEqual(data['gain_per_completed_box'],2/3)
        self.assertEqual(data['censored_or_failed_workers'],1)

    def test_more_gains_do_not_imply_better_efficiency(self):
        rows=population()
        for r in rows:
            if r['arm']!='top':r['discovery_worker_seconds']=1000
        result=comparison(rows)
        self.assertGreater(result['top_only']['gain_per_discovery_second'],result['equal_three_arm_portfolio']['gain_per_discovery_second'])
        self.assertEqual(result['recommendation'],'INCONCLUSIVE')

    def test_diversification_requires_family_robustness(self):
        rows=population()
        self.assertTrue(comparison(rows)['diversification_criterion_met'])
        for r in rows:
            if r['arm']!='top':r['certified_gain']=20 if r['family']==0 else 0
        result=comparison(rows)
        self.assertGreater(result['equal_three_arm_portfolio']['gain_per_discovery_second'],result['top_only']['gain_per_discovery_second'])
        self.assertFalse(result['diversification_criterion_met'])

    def test_missing_exposure_and_lower_completion_block_recommendation(self):
        rows=population();rows[1]['completed_boxes']=None
        self.assertIsNone(summarize(rows)['completion_fraction'])
        self.assertEqual(comparison(rows)['recommendation'],'INCONCLUSIVE')
        rows[1]['completed_boxes']=42
        self.assertFalse(comparison(rows)['diversification_criterion_met'])


if __name__=='__main__':unittest.main()

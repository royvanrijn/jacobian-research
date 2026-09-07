import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
import report_nearcut60_mw16_experiment as report
import audit_nearcut60_mw16_accounting as independent


class NearcutAccountingTests(unittest.TestCase):
    def test_independent_totals_agree_for_complete_and_partial_evidence(self):
        for gain,boxes,status in [(0,43,'COMPLETE_DECLARED_POINT_ATTEMPT'),(3,7,'POINT_FAILED_OR_CENSORED'),(None,None,'POINT_FAILED_OR_CENSORED')]:
            with self.subTest(gain=gain):
                actual=report.summarize([{'certified_gain':gain,'completed_boxes':boxes,'discovery_worker_seconds':600.25,'verification_seconds':30.5,'worker_status':status}])
                expected=independent.totals([{'gain':gain,'boxes':boxes,'discovery':600.25,'verification':30.5,'worker_status':status}])
                for k,v in expected.items():independent.same(actual[k],v)
                if gain is None:
                    self.assertIsNone(actual['gain_per_discovery_second'])
                    self.assertEqual(actual['unresolved_gain_curves'],1)


if __name__=='__main__':unittest.main()

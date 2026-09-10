"""Pure policy regressions for the X1092 Elkies-pair v2 scheduler."""
from pathlib import Path
import sys
import unittest

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import elkies_pair_policy as p


class SchedulerV2Tests(unittest.TestCase):
    def test_controls_preserve_ranked_candidate(self):
        self.assertEqual(p.effective_control_count(0,2),0)
        self.assertEqual(p.effective_control_count(1,2),0)
        self.assertEqual(p.effective_control_count(2,2),1)
        self.assertEqual(p.effective_control_count(3,2),2)
        self.assertEqual(p.effective_control_count(10,0),0)

    def test_control_reduction_depends_only_on_population_size(self):
        rows=[{'id':'a'},{'id':'b'}]
        n=p.effective_control_count(len(rows),2)
        controls,ranked=p.split_controls(rows,n)
        self.assertEqual(len(controls),1)
        self.assertEqual(len(ranked),1)

    def test_carrier_quality_prefers_primary_yield(self):
        a=p.carrier_quality_summary([100,200,5000],1024)
        b=p.carrier_quality_summary([100,200,300,4000],1024)
        self.assertGreater(p.carrier_quality_key(b),p.carrier_quality_key(a))

    def test_extended_yield_breaks_primary_tie(self):
        a=p.carrier_quality_summary([100,200,8000],1024)
        b=p.carrier_quality_summary([100,200,3000],1024)
        self.assertGreater(p.carrier_quality_key(b),p.carrier_quality_key(a))

    def test_same_cover_is_only_a_tie_break(self):
        weak=p.carrier_quality_summary([100,9000],1024)
        strong=p.carrier_quality_summary([100,200,300],1024)
        self.assertGreater(p.carrier_quality_key(strong,False),p.carrier_quality_key(weak,True))
        equal=p.carrier_quality_summary([100,200],1024)
        self.assertGreater(p.carrier_quality_key(equal,True),p.carrier_quality_key(equal,False))

    def test_empty_quality_loses(self):
        empty=p.carrier_quality_summary([],1024)
        usable=p.carrier_quality_summary([1024],1024)
        self.assertGreater(p.carrier_quality_key(usable),p.carrier_quality_key(empty))


if __name__=='__main__':unittest.main()

"""Scheduling, evidence boundaries and supervisor renewal regressions."""
import copy
from fractions import Fraction
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

CAS=Path(__file__).resolve().parents[1]/'cas';sys.path.insert(0,str(CAS))
import parent_foundry_policy as policy
import run_parent_foundry as controller


class PolicyTests(unittest.TestCase):
    def test_guardian_really_starts_a_second_clean_run(self):
        with tempfile.TemporaryDirectory() as name:
            folder=Path(name);calls=[]
            class Process:
                pid=1
                def wait(self):
                    if len(calls)==2:(folder/'STOP').touch()
                    return 0
            def start(*args,**kwargs):calls.append(args);return Process()
            with patch.object(controller,'guard',return_value=({},folder)), \
                 patch.object(controller.subprocess,'Popen',side_effect=start), \
                 patch.object(controller,'process_info',return_value={'start_token':'test'}), \
                 patch.object(controller.time,'sleep'):
                controller.guardian(folder)
            self.assertEqual(len(calls),2)
            self.assertEqual(json.loads((folder/'renewal.json').read_text())['action'],
                             'START_NEXT_CONTROLLER_RUN')

    def test_clean_run_and_crash_renew_without_daily_budget(self):
        for rc in (0,1,75,-9):
            self.assertTrue(policy.guardian_should_restart(False,rc))
            self.assertFalse(policy.guardian_should_restart(True,rc))

    def test_parameters_are_reproducible_distinct_and_not_endpoints(self):
        for family in ('test-x948','test-x1092'):
            values=[policy.parameter(family,i) for i in range(36)]
            self.assertEqual(len(values),len(set(values)))
            self.assertEqual(values,[policy.parameter(family,i) for i in range(36)])
            self.assertNotIn('0',values)
            for i,q in enumerate(values[:12]):
                t=Fraction(q);height=max(abs(t.numerator),t.denominator)
                self.assertTrue(2<=height<=8 if i<8 else 17<=height<=64)

    def test_censored_slots_are_not_zero_rank_observations(self):
        rows=[{'rank_lower_bound':24,'certified_jump_lower_bound':8,'exposure_complete':True},
              {'rank_lower_bound':None,'certified_jump_lower_bound':None,'exposure_complete':False},
              {'rank_lower_bound':17,'certified_jump_lower_bound':1,'exposure_complete':False}]
        tails=policy.tails(rows)
        self.assertEqual(tails['certified_slots'],2)
        self.assertEqual(tails['8']['fraction'],'1/3')
        self.assertEqual(tails['8']['incomplete_slot_upper_fraction'],'3/3')
        self.assertEqual(tails['incomplete_exposures'],2)

    def test_parent_promotion_requires_its_panel(self):
        rows=[{'certified_jump_lower_bound':None} for _ in range(12)]
        rows[0]['certified_jump_lower_bound']=8
        self.assertFalse(policy.promote(rows[:11]));self.assertTrue(policy.promote(rows))
        rows[0]['certified_jump_lower_bound']=0
        self.assertFalse(policy.promote(rows))

    def test_momentum_can_outrank_a_stalled_high_rank_curve(self):
        lively={'rank':25,'generic_rank':16,'calls':64,'stale_calls':0,'last_gain':9}
        stalled={'rank':28,'generic_rank':16,'calls':2000,'stale_calls':255,'last_gain':0}
        self.assertGreater(policy.utility(lively),policy.utility(stalled))
        stalled['stale_calls']=256;self.assertFalse(policy.continuation(stalled))

    def test_incremental_exploitation_is_bounded_per_job(self):
        self.assertEqual([policy.exploit_allowance({'batches':k}) for k in (0,1,3,100)],
                         [64,96,160,256])

    def test_scheduler_constructs_real_inputs_not_masked_cores(self):
        state={'parents':{},'jobs':{},'cursor':0}
        generic={'proposals':[{'source_family':'det1092','priority':12,'trace_word':[1]*17}],
                 'sources':[{'family':'det1092','generic_rank_lower_bound':17}]}
        request=controller.schedule(Path('/tmp'),state,{},Path('/tmp'),generic)
        self.assertEqual(request['kind'],'construct')
        self.assertEqual(request['source']['generic_rank_lower_bound'],17)
        self.assertEqual(state['cursor'],1)

    def test_baseline_panel_does_not_receive_result_driven_extra_slots(self):
        state={'parents':{'baseline':{'baseline':True,'dispatched':0,'target_slots':12,
            'path':'parent.json','sha256':'hash','generic_rank':17}},'jobs':{},'cursor':0,
            'fibres':{},'seconds':{'panel':1,'exploit':0},'evaluation_dispatches':0}
        request=controller.schedule(Path('/tmp'),state,{'exploit_share':.3,'baseline_every':5,'height':125000},
            Path('/tmp'),{'proposals':[]})
        self.assertEqual(request['allowance'],64)
        self.assertEqual(state['parents']['baseline']['target_slots'],12)


if __name__=='__main__':unittest.main()

"""Metric/stage separation and blind calibration protocol regressions."""
from dataclasses import dataclass
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.chart_policy import ChartPolicy, RepresentationPipeline, calibration_protocol, rank_calibration


class PolicyTests(unittest.TestCase):
    def test_weight_changes_metric_without_changing_enumerator(self):
        generic, quotient = [[1,0],[0,0]], [[0,0],[0,1]]
        first, second = ChartPolicy(), ChartPolicy(quotient_weight='16')
        self.assertEqual(first.metric(generic,quotient),((1,0),(0,1)))
        self.assertEqual(second.metric(generic,quotient),((1,0),(0,16)))
        self.assertEqual(first.enumeration_backend,second.enumeration_backend)
        self.assertNotEqual(first.key,second.key)

    def test_representation_stages_are_separate_and_update_state(self):
        @dataclass(frozen=True)
        class State:
            key: str = 'old'
            arithmetic: str = 'curve'
            basis: tuple = ((1,2),)
        calls = []
        def normalize(state,limits):
            calls.append('raw');return state.model if hasattr(state,'model') else 'raw-model'
        def parameterize(state,model,centre,policy,limits):
            calls.append((model,centre));return 'chart'
        def enumerate_points(state,chart,limits):
            calls.append(chart);return State('new',basis=(*state.basis,(3,4))),{'exact_points':True}
        pipeline=RepresentationPipeline(normalizers={'raw':normalize},parameterizations={'pointed-quartic':parameterize},
                                        enumerators={'gmp-pointed-sieve':enumerate_points})
        state,record=pipeline.run(State(),ChartPolicy(),centre=(1,0),limits={'height':10})
        self.assertEqual(calls,['raw',('raw-model',(1,0)),'chart'])
        self.assertEqual(state.key,'new');self.assertEqual(len(record['measurements']),3)
        self.assertFalse(record['mathematical_exclusion'])

    def test_calibration_requires_full_blind_cells_and_fixed_protocol(self):
        policies=[ChartPolicy(),ChartPolicy(quotient_weight='16')]
        protocol=calibration_protocol(panel=['blind-a','blind-b'],policies=policies,limits={'height':10},
                                      outcome_commitment='frozen-labels',controls=['held-out'])
        cells=[{'case':case,'policy':policy.key,'protocol_hash':protocol['protocol_hash'],
                'certified_independent_recoveries':int(policy.quotient_weight=='16'),'wall_seconds':'1/10'}
               for case in protocol['panel'] for policy in policies]
        self.assertEqual(rank_calibration(protocol,cells)[0]['policy_key'],policies[1].key)
        with self.assertRaises(ValueError):rank_calibration(protocol,cells[:-1])
        with self.assertRaises(ValueError):rank_calibration(protocol,cells+[cells[0]])
        protocol['limits']['height']=20
        with self.assertRaises(ValueError):rank_calibration(protocol,cells)


if __name__=='__main__':unittest.main()

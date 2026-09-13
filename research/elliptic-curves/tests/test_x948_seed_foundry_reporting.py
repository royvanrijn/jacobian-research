"""Regression checks for the experiment's scientific endpoint accounting."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from run_x948_seed_foundry import seed_summary
from report_x948_seed_foundry import acquired_after_23, amplification_summary


def row(rank, *, complete=True, initial=16, phase='seed', address=0):
    return {'family': 'F', 'phase': phase, 'valid': True, 'address_index': address,
            'cost': {'process_tree_cpu_seconds': 2},
            'result': {'rank_lower_bound': rank, 'initial_rank': initial,
                       'binary_outcome_complete': complete, 'gain_timeline': []}}


class EndpointTests(unittest.TestCase):
    def test_missing_preparation_is_not_a_negative(self):
        result = seed_summary([row(16), row(18), row(None, complete=False)], [{'family':'F'}])[0]
        self.assertEqual((result['seeded'], result['complete_binary_outcomes'], result['unresolved']), (1,2,1))
        self.assertEqual(result['seed_fraction'], 0.5)
        self.assertEqual(result['upper_all_slot_fraction'], 2/3)

    def test_controls_do_not_enter_prospective_denominators(self):
        result = seed_summary([row(18,address=None), row(16)], [{'family':'F'}])[0]
        self.assertEqual((result['slots'], result['seeded']), (1,0))

    def test_one_large_cloud_is_not_a_later_deep_gain(self):
        result = {'initial_rank':17, 'gain_timeline':[
            {'call':1,'after':23}, {'call':1,'after':24}]}
        self.assertFalse(acquired_after_23(result))
        result['gain_timeline'].append({'call':2,'after':25})
        self.assertTrue(acquired_after_23(result))

    def test_reaching_threshold_in_seed_is_separate_from_amplifying(self):
        result = amplification_summary([row(24,initial=24,phase='amplify')], [{'family':'F'}])[0]
        self.assertEqual(result['p_amp'], 0)
        self.assertEqual(result['p_reach_23'], 1)
        self.assertEqual(result['already_in_seed_23'], 1)
        self.assertEqual(result['acquiring_after_23'], 0)


if __name__ == '__main__':
    unittest.main()

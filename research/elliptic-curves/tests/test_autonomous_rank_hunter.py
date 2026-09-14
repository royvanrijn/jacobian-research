from pathlib import Path
import sys

CAS = Path(__file__).resolve().parents[1] / 'cas'
sys.path.insert(0, str(CAS))

import autonomous_rank_hunter_policy as p


def test_productive_m27_continues():
    s = {'rank_lower_bound': 27, 'complement_calls': 100,
         'gain_calls': [9, 20, 21, 44, 67], 'complement_directions': 8,
         'bank_calls': 100, 'bank_last_gain_call': 67, 'bank_generation': 0,
         'has_suffix': True, 'seed_cloud_directions': 1}
    assert p.continuation_decision(s)['action'] == 'continue_bank'


def test_stalled_m26_does_not_grind_forever():
    s = {'rank_lower_bound': 26, 'complement_calls': 500,
         'gain_calls': [9, 10, 31, 54, 97], 'complement_directions': 6,
         'bank_calls': 500, 'bank_last_gain_call': 97, 'bank_generation': 2,
         'has_suffix': True}
    assert p.continuation_decision(s)['action'] == 'retire'


def test_stale_productive_curve_pivots_representation():
    s = {'rank_lower_bound': 27, 'complement_calls': 300,
         'gain_calls': [10, 20, 35, 60, 90], 'complement_directions': 5,
         'bank_calls': 300, 'bank_last_gain_call': 90, 'bank_generation': 0,
         'has_suffix': True}
    assert p.continuation_decision(s)['action'] == 'new_bank'


def test_breadth_is_reserved():
    ordinary = {'rank_lower_bound': 27, 'complement_calls': 100,
                'gain_calls': [67], 'complement_directions': 3,
                'bank_calls': 100, 'bank_last_gain_call': 67, 'bank_generation': 0,
                'has_suffix': True}
    assert p.exploration_slots(4, [ordinary]) == 2


def test_only_hot_near_record_can_use_three_exploit_slots():
    hot = {'rank_lower_bound': 29, 'complement_calls': 500,
           'gain_calls': [350, 430, 490], 'complement_directions': 10,
           'bank_calls': 100, 'bank_last_gain_call': 90, 'bank_generation': 1,
           'has_suffix': True}
    assert p.exploration_slots(4, [hot]) == 1


def test_selection_schedule_never_excludes_high_height():
    lanes = {p.lane_for_dispatch(i) for i in range(12)}
    assert any(x.startswith('low_') for x in lanes)
    assert any(x.startswith('high_') for x in lanes)
    assert {'low_score', 'low_spread', 'low_hash', 'high_score', 'high_spread', 'high_hash'} <= lanes

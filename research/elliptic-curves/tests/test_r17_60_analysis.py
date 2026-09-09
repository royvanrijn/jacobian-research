from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from analyze_r17_60_panel import summarize


def row(seed, calls, cloud, final, last):
    return dict(first_seed=seed, seed_calls=calls, seed_cloud_rank=cloud, rank_lower_bound=final,
                complement_calls=100 if seed else 0, complement_added_directions=final-cloud,
                last_complement_gain_call=last, point_timeouts=0, map_timeouts=0, stop_reason='bounded')


def test_seed_cloud_is_not_mistaken_for_amplification():
    r = summarize([row(True, 1, 22, 22, None), row(True, 3, 18, 20, 100), row(False, 86, 17, 17, None)])
    assert r['cases'] == 3 and r['first_seeds'] == 2 and r['first_call_seeds'] == 1
    assert r['seed_cloud_added_directions'] == 6 and r['complement_added_directions'] == 2
    assert r['seeded_without_complement_gain'] == 1 and r['late_gain_cases'] == 1
    assert r['median_calls_to_seed_conditional'] == '2'
    assert r['mean_rank_lower_bound'] == '59/3'


def test_empty_or_all_miss_strata_have_no_conditional_seed_median():
    assert summarize([])['mean_rank_lower_bound'] is None
    assert summarize([row(False, 86, 17, 17, None)])['median_calls_to_seed_conditional'] is None

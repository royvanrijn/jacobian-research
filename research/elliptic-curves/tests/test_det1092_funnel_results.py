from pathlib import Path
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from det1092_funnel_results import compare_bound, benchmarks


def test_large_upper_bound_is_not_an_exclusion():
    assert compare_bound(10,10**100,100)=='UNRESOLVED_COMPARISON'
    assert compare_bound(101,10**100,100)=='PROVED_ABOVE_RECORDED_MINIMUM'
    assert compare_bound(10,99,100)=='PROVED_BELOW_RECORDED_MINIMUM'
    assert compare_bound(100,100,100)=='UNRESOLVED_COMPARISON'


def test_missing_benchmarks_remain_missing():
    rows=[dict(id=1,rank_lower_bound=23,conductor='99'),
          dict(id=2,rank_lower_bound=24,conductor=None),
          dict(id=3,rank_lower_bound=22,conductor='1')]
    assert benchmarks(rows,23)['reported_minimum']=='99'
    assert benchmarks(rows,23)['missing_conductor_ids']==[2]
    assert benchmarks(rows,24)['reported_minimum'] is None
    assert compare_bound(10,100,None)=='NO_RECORDED_BENCHMARK'


def test_invalid_conductor_intervals_rejected():
    for values in [(0,10,5),(10,5,7),(1,3,0)]:
        with pytest.raises(ValueError):
            compare_bound(*values)

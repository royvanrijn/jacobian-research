import copy
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
SCRIPT=ROOT/'elliptic-curves/cas/analyze_curve302_unlock_seed_closure.py'
spec=importlib.util.spec_from_file_location('unlock_seed_closure',SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def toy_states(dim=2):
    # 00: add0 cost10, add1 cost7; 01: add1 cost3; 10: add0 cost4.
    return [
        {'state_mask':0,'retained_numerators':[10,7]},
        {'state_mask':1,'retained_numerators':[None,3]},
        {'state_mask':2,'retained_numerators':[4,None]},
        {'state_mask':3,'retained_numerators':[None,None]},
    ]


def test_minimax_unseeded_prefers_direction1_then0(monkeypatch):
    monkeypatch.setattr(m,'DIM',2)
    value,order=m.minimax_from(toy_states(),0)
    assert value==7
    assert order==[1,0]


def test_single_seed_changes_bottleneck(monkeypatch):
    monkeypatch.setattr(m,'DIM',2)
    value,order=m.minimax_from(toy_states(),1)
    assert value==3
    assert order==[1]


def test_threshold_closure_is_fixed_point(monkeypatch):
    monkeypatch.setattr(m,'DIM',2)
    mask,waves=m.closure_at(toy_states(),0,7)
    assert mask==3
    assert waves==[[1],[0]]


def test_seed_row_exposes_both_order_schema_names(monkeypatch):
    monkeypatch.setattr(m,'DIM',2)
    monkeypatch.setattr(m,'GENERIC_RANK',17)
    row=m.row_for_seed(toy_states(),['a','b'],1,{'toy':7})
    assert row['optimal_followup_order']==['a']
    assert row['optimal_completion_order']==row['optimal_followup_order']


def sample_payload():
    seed={
        'seed_direction':'a','start_rank':18,
        'minimax_bottleneck_numerator':7,
        'minimax_bottleneck_scaled_height':7/4_000_000.0,
        'optimal_followup_order':['b'],
        'optimal_completion_order':['b'],
        'calibrated_closures':{},
    }
    return {
        'schema':'elliptic-curves.curve302-unlock-seed-closure.v2',
        'status':'PASS_RETROSPECTIVE_SINGLE_SEED_CLOSURE_CENSUS',
        'input':{'x':'y'},'direction_ids':['a','b'],'thresholds':{},
        'unseeded':copy.deepcopy(seed),'single_seed_ranking':[copy.deepcopy(seed)],
        'xi_direction':None,'xi_result':None,'interpretation_boundary':'x',
        'reproducing_command':'new',
    }


def test_legacy_missing_order_alias_is_semantically_equal():
    expected=sample_payload()
    legacy=copy.deepcopy(expected)
    legacy['schema']='elliptic-curves.curve302-unlock-seed-closure.v1'
    legacy['reproducing_command']='old'
    legacy['unseeded'].pop('optimal_completion_order')
    legacy['single_seed_ranking'][0].pop('optimal_completion_order')
    assert m.normalized_semantics(legacy)==m.normalized_semantics(expected)


def test_changed_math_is_not_accepted_as_legacy_compatibility():
    expected=sample_payload()
    changed=copy.deepcopy(expected)
    changed['single_seed_ranking'][0]['minimax_bottleneck_numerator']=8
    with pytest.raises(ArithmeticError):
        m.validate_existing(changed,expected)

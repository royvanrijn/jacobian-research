import importlib.util
from pathlib import Path

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

import importlib.util
from pathlib import Path
import sys

CAS = Path(__file__).resolve().parents[1]/'cas'
spec = importlib.util.spec_from_file_location('broad', CAS/'run_broad_mw17_search.py')
broad = importlib.util.module_from_spec(spec); spec.loader.exec_module(broad)


def rows(n=128):
    out=[]
    for i in range(n):
        out.append({'index':i,'control':i%8==0,'control_order':i//8 if i%8==0 else None,
                    'nonsingular':True,'score_units':1000-i,'model_bits':10+(i%7),
                    'discriminant_bits':20+(i%11),'smooth_prime_count':30-(i%5)})
    return out


def test_disjoint_fixed_arm_sizes():
    selected=broad.choose_arms(rows(),20,8,10)
    assert len(selected)==38
    assert len({r['index'] for r in selected})==38
    assert sum(r['broad_arm']=='ranked' for r in selected)==20
    assert sum(r['broad_arm']=='control' for r in selected)==8
    assert sum(r['broad_arm']=='diversity' for r in selected)==10


def test_controls_do_not_depend_on_score():
    a=broad.choose_arms(rows(),20,8,10)
    changed=rows()
    for r in changed:r['score_units']=-r['score_units']
    b=broad.choose_arms(changed,20,8,10)
    ca=[r['index'] for r in a if r['broad_arm']=='control']
    cb=[r['index'] for r in b if r['broad_arm']=='control']
    assert ca==cb


def test_singular_rows_never_dispatched():
    data=rows(); data[0]['nonsingular']=False; data[1]['nonsingular']=False
    selected=broad.choose_arms(data,20,8,10)
    assert 0 not in {r['index'] for r in selected}
    assert 1 not in {r['index'] for r in selected}

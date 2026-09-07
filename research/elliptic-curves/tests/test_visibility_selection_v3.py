import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from visibility_selection_v3 import cheap_shortlist, final_shortlist, tied_take

def test_all_boundary_ties_survive():
    assert len(cheap_shortlist([7]*100,list(range(100))))==100
    assert tied_take(list(range(5)),lambda i:[1,2,2,2,3][i],2)=={0,1,2,3}

def test_mask_enumeration_is_not_a_tiebreak():
    for n in (1,2,16,64,8192):
        norms=[(i*i+17*i)%100003 for i in range(n)]
        keys=[i*123456789+17 for i in range(n)]
        left={keys[i] for i in cheap_shortlist(norms,keys)}
        right={keys[::-1][i] for i in cheap_shortlist(norms[::-1],keys[::-1])}
        assert left==right

def test_finalists_preserve_equal_profiles_and_reverse_order():
    rows=[{'metric_norm':i%3,'quartic_bits':100+i%5,'quartic_max_bits':30+i%2,
           'multiplicity':1,'fingerprint':i+1} for i in range(100)]
    assert {r['fingerprint'] for r in final_shortlist(rows)}=={r['fingerprint'] for r in final_shortlist(rows[::-1])}
    equal=[dict(r,metric_norm=7,quartic_bits=100,quartic_max_bits=20) for r in rows]
    assert len(final_shortlist(equal))==100

"""Mathematical partition and policy isolation regressions; no point search."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))

from finite_cancellation_validation_features import q_tree, gcd_tree, model_features, Q_NAMES, REAL
from finite_cancellation_validation_replay import verify_tree


def test_two_adic_square_units_and_infinity_are_certified():
    for q in ([1,0,0,0,1],[3,0,0,0,0],[4,0,0,0,0],[0,0,0,0,1]):
        tree=q_tree(q,2)
        assert any(t['kind']=='infinity' for t in tree['leaves'])
        assert verify_tree(tree,None,None,q)==len(tree['leaves'])


def test_q_partition_is_separate_from_cancellation_refinement():
    q=[1,0,0,0,1];n=[0,0,0,0,1];d=[0,0,1,0,0]
    for p in (2,3,5,7):
        qt=q_tree(q,p);snapshot=repr(qt)
        gt=gcd_tree(n,d,qt)
        assert repr(qt)==snapshot
        verify_tree(qt,None,None,q)
        verify_tree(gt,n,d,q,True)
        assert gt['unknown_mass']>=qt['unknown_mass']


def test_fitted_arms_share_identical_q_and_real_inputs():
    n=[1,2,3,4,5];d=[7,-3,1,0,1];q=[1,0,0,0,1]
    small,sl=model_features(n,d,q,False);full,fl=model_features(n,d,q,True)
    assert small=={k:full[k] for k in small} and sl['q']==fl['q']
    other,_=model_features([p*3 for p in n],d,q,False)
    assert {k:small[k] for k in Q_NAMES}=={k:other[k] for k in Q_NAMES}


def test_unresolved_mass_is_not_silently_soluble():
    tree=q_tree([0,0,0,0,0],3,depth=2,cap=4)
    assert tree['soluble_mass']==0 and abs(tree['unknown_mass']-1)<1e-12
    verify_tree(tree,None,None,[0,0,0,0,0])

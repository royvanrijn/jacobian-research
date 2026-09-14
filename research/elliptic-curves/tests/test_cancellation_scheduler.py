"""Mathematical/state-boundary regressions, not timing assertions."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from cancellation_scheduler import AnchorExposure, CostModel, neighbour_law, root_distribution
from finite_cancellation_validation_replay import verify_tree


def test_q_partition_keeps_same_prime_balls_disjoint_and_unknown_explicit():
    q=[9,0,6,0,1]  # (t^2+3)^2; fixture is for the local polynomial partition.
    distribution,tree=root_distribution(q,3,[0,1])
    assert set(distribution)=={0,1,None}
    assert sum(distribution.values())==pytest.approx(1)
    assert distribution[0]>0
    verify_tree(tree,[1,0,0,0,0],[1,0,0,0,0],q)


def test_prime_neighbour_exact_height_and_gcd_law_in_both_balls():
    n=[1,2,0,0,1];d=[2,1,0,1,1];p=3;T=[3,1,0,1]
    e=neighbour_law(n,d,T,p)
    from finite_cancellation_features import ev
    from search_observability import transform
    nn=[int(x)//p**e for x in transform(n,tuple(map(F,T)))]
    dd=[int(x)//p**e for x in transform(d,tuple(map(F,T)))]
    deltas=set()
    for m in range(-8,9):
        for v in range(1,9):
            if gcd(m,v)!=1:continue
            M,V=3*m+v,v;g=gcd(M,V);delta=int(g==p);deltas.add(delta)
            assert g in (1,p)
            M//=g;V//=g
            assert F(gcd(ev(nn,m,v),ev(dd,m,v)),gcd(ev(n,M,V),ev(d,M,V)))==F(p)**(4*delta-e)
            adj=(M-V,3*V);ag=gcd(*adj)
            assert ag==(p if (M-V)%p==0 else 1)
            assert max(abs(x//ag) for x in adj)==max(abs(m),v)
    assert deltas=={0,1}


def test_timeout_does_not_erase_overlap_and_complete_search_does():
    fit={'anchor_priors':[.2]*3,'radial_heights':[1000,10000,100000],
         'backend_cpu_at_125000':.5,'per_call_overhead_cpu':.01,'independent_certificate_cpu_prior':1.}
    prepared={'weights':[.5,.5],'ratios':[[1.,.1],[1.,10.]],'models':[{},{}]}
    anchor=AnchorExposure(prepared,fit,0)
    anchor.observe(0,False)
    assert np.all(anchor.covered==0)
    assert 0 not in [o[2] for o in anchor.options(CostModel(fit))]
    anchor.observe(2,True)
    assert np.all(anchor.covered==1)
    assert not list(anchor.options(CostModel(fit)))


def test_bad_prime_neighbour_fails_closed():
    with pytest.raises(ArithmeticError):neighbour_law([1]*5,[2]*5,[4,0,0,1],3)

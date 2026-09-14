"""Mathematical regressions for local probability and model covariance."""
import sys
from fractions import Fraction as F
from math import gcd
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from finite_cancellation_features import local_tree, leaf_for, ev, integral_pair, affine, vp


def test_soluble_projective_measure_and_infinity():
    # N=x^4+1, D=x^4: common cancellation is impossible, including
    # infinity. Verify every asserted square class for q=x^4+1.
    n=[1,0,0,0,1];d=[0,0,0,0,1];q=[1,0,0,0,1]
    for p in (2,3,5,7):
        tree=local_tree(n,d,q,p)
        assert sum(F(z['mass']) for z in tree['leaves'])==1
        for a in range(-25,26):
            for b in range(26):
                if gcd(a,b)!=1:continue
                z=leaf_for(tree,a,b)
                assert z['vg_exact'] and z['vg']==0
                actual=ev(q,a,b)
                if z['status']=='square':
                    v=vp(actual,p)
                    assert actual==0 or (v%2==0 and (actual//p**v%8==1 if p==2 else pow(actual//p**v%p,(p-1)//2,p)==1))


def test_neighbour_cancellation_has_no_other_prime_support():
    n=[4,0,-4,4,0];d=[4,-8,0,0,1]  # E:y^2=x^3-x+1, anchor(0,1)
    for p in (3,5,7):
        for r in range(p):
            nn,dd=integral_pair(affine(n,p,r),affine(d,p,r))
            for a in range(-10,11):
                for b in range(1,11):
                    if gcd(a,b)!=1:continue
                    x,y=p*a+r*b,b;content=gcd(x,y);x//=content;y//=content
                    old=gcd(ev(n,x,y),ev(d,x,y));new=gcd(ev(nn,a,b),ev(dd,a,b))
                    while old%p==0:old//=p
                    while new%p==0:new//=p
                    assert old==new


def test_unresolved_mass_is_not_soluble_mass():
    # A capped tree around a multiple root must retain uncertainty.
    q=[0,0,0,0,1];n=[0,0,0,0,1];d=[0,0,0,1,0]
    tree=local_tree(n,d,q,3,depth=1)
    assert tree['unknown_mass']>0
    assert tree['soluble_mass']+tree['unknown_mass']<=1
    assert leaf_for(tree,0,1)['status']=='unknown'


def test_joint_primitive_normalization_removes_artificial_gcd():
    n=[4,0,-4,4,0];d=[4,-8,0,0,1]
    assert integral_pair(n,d)==integral_pair([F(z)*F(25,49) for z in n],[F(z)*F(25,49) for z in d])


def test_drop_in_selector_matches_every_retained_cpu_choice():
    import json
    from importlib.machinery import SourceFileLoader
    from finite_cancellation_corpus import OUT
    module=SourceFileLoader('tested_finite_selector',str(Path(__file__).resolve().parents[1]/'cas/lean_finite_cancellation_pari_mapping.sage')).load_module()
    plan=json.loads((OUT/'cpu/protocol.json').read_text());count=0
    for p in (OUT/'cpu/arms').glob('*/finite_selector/chart-*.json'):
        row=json.loads(p.read_text());selection=row['selection']
        index,scores=module.select(selection['features'],plan['fit'])
        assert scores==selection['scores'] and index==min(range(len(scores)),key=lambda i:scores[i])
        count+=1
    assert count==237

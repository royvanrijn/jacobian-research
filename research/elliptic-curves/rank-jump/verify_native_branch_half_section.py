#!/usr/bin/env python3
"""Pure rational quadratic-algebra replay of all native branch halves."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
import retrospective as r
import native_branch_half_section as run
import branch_divisibility_capacity as branch
import verify_branch_divisibility_capacity as prior
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_native_branch_half_section_verification_v1.json'


def compute():
    d=r.read(run.INPUT);out=r.read(run.OUTPUT);old={x['label']:x for x in r.read(prior.OUTPUT)['rows']}
    for obj in (d,out):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    rows=[];operations=0
    for c,result in zip(d['rows'],out['rows']):
        q=list(map(Q,c['q']));K=Algebra([x/q[-1] for x in q]);b=K.elt([0,1])
        ev=lambda coeff:K.evaluate(list(map(Q,coeff)),b)
        A,B=ev(d['A']),ev(d['B']);zero=K.elt([0]);one=K.elt([1])
        scale=lambda a,n:tuple(n*x for x in a)
        sub=lambda a,b:K.add(a,K.neg(b))
        def check(P):
            if P is not None:
                x,y=P;assert K.mul(y,y)==K.add(K.add(K.power(x,3),K.mul(A,x)),B)
        def add(P,Q):
            nonlocal operations
            if P is None:return Q
            if Q is None:return P
            x,y=P;xx,yy=Q
            if x==xx:
                if K.add(y,yy)==zero:return None
                assert y==yy
                slope=K.mul(K.add(scale(K.mul(x,x),3),A),K.inverse(scale(y,2)))
            else:slope=K.mul(sub(yy,y),K.inverse(sub(xx,x)))
            xxx=sub(sub(K.mul(slope,slope),x),xx);yyy=sub(K.mul(slope,sub(x,xxx)),y)
            ans=(xxx,yyy);check(ans);operations+=1;return ans
        def multiply(P,n):
            if n<0:P=(P[0],K.neg(P[1]));n=-n
            ans=None
            while n:
                if n&1:ans=add(ans,P)
                n//=2
                if n:P=add(P,P)
            return ans
        points=[(ev(s['x']),ev(s['y'])) for s in d['sections']]
        half=(ev(c['half_x']),ev(c['half_y']));check(half)
        trace=None
        for n,P in zip(c['word'],points):check(P);trace=add(trace,multiply(P,n))
        assert trace is not None and add(half,half)==trace
        enc=lambda P:[list(map(str,x)) for x in P]
        assert enc(half)==result['half_point'] and enc(trace)==result['trace_point']
        mask=sum((n%2)<<i for i,n in enumerate(c['word']))
        assert mask==result['kernel_mask']==old[c['label']]['finite_kernel_mask'] and mask
        rows.append({'label':c['label'],'kernel_mask':mask,'exact_branch_kernel_dimension':1})
    assert len(rows)==37
    return {'schema':'rank-jump.native-branch-half-section-verification.v1','status':'PASS','rows':rows,
        'rational_group_operations':operations,'branch_kernels_exact':37,
        'method':'Rational multiplication and inversion modulo the monic quadratic, explicit chord/tangent group law; no Sage elliptic or number-field arithmetic.',
        'bindings':branch.bindings([Path(__file__),run.INPUT,run.OUTPUT,prior.OUTPUT,
            Path(__file__).with_name('verify_unpointed_governing_norm.py')])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['branch_kernels_exact'],'exact branch kernels;',result['rational_group_operations'],'rational group operations')

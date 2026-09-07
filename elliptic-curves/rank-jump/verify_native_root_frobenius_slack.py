#!/usr/bin/env python3
"""Independent finite-field factor counts and reciprocal mod-two algebra."""
import argparse
from pathlib import Path
from math import isqrt
import retrospective as r
import native_root_frobenius_slack as run
import branch_divisibility_capacity as branch

OUTPUT=r.OUT/'rank_jump_native_root_frobenius_slack_verification_v1.json'


def compute():
    from sage.all import QQ,GF,PolynomialRing
    d=r.read(branch.INPUT);out=r.read(run.OUTPUT)
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    R=PolynomialRing(QQ,'t');A=R(d['A']);B=R(d['B']);D=-4*A**3-27*B**2
    expected=[];counts=[];fibres=0
    for p in range(5,1010):
        if any(p%i==0 for i in range(2,isqrt(p)+1)):continue
        P=PolynomialRing(GF(p),'x');x=P.gen()
        try:a,b,delta=P(A),P(B),P(D)
        except (ValueError,ZeroDivisionError):continue
        if delta.degree()!=24 or delta.gcd(delta.derivative())!=1 or delta.gcd(a)!=1:continue
        expected.append(p);row=out['rows'][len(expected)-1];assert row['p']==p
        finite=[len((x**3+a(t)*x+b(t)).roots(multiplicities=False)) for t in GF(p)]
        infinity=sorted(int(z) for z in (x**3+a[8]*x+b[12]).roots(multiplicities=False))
        assert finite==row['finite_fibre_counts'] and infinity==row['infinity_roots']
        total=sum(finite)+len(infinity);assert total==row['point_count'] and p+1-total==row['frobenius_trace']
        assert bool(total%2)==row['odd_point_count'];fibres+=p+1
        counts.append({'p':p,'point_count':total,'trace_mod2':int((p+1-total)%2)})
        if len(expected)==12:break
    F=PolynomialRing(GF(2),'T');T=F.gen();possibilities=[]
    for a in (0,1):
        residual=T**3+a*T**2+a*T+1;full=(T+1)**17*residual
        assert full.reverse()==full
        q=full;multiplicity=0
        while q(1)==0:q=q//(T+1);multiplicity+=1
        trace=int(full[19]);assert multiplicity==(18 if trace else 20)
        possibilities.append({'trace_mod2':trace,'residual_coefficients':list(map(int,residual.list())),
            'eigenvalue_one_algebraic_multiplicity':multiplicity})
    odd=[x['p'] for x in counts if x['trace_mod2']]
    assert odd==out['odd_count_primes']==[211,229] and out['global_pool_upper_bound']==18
    assert out['single_support_twist_upper_bound']==2 and out['multiple_support_twist_upper_bound']==1
    return {'schema':'rank-jump.native-root-frobenius-slack-verification.v1','status':'PASS',
        'counts':counts,'finite_and_infinity_fibres_verified':fibres,'reciprocal_possibilities':possibilities,
        'global_pool_dimension_interval':[17,18],'odd_count_primes':odd,
        'method':'Sage finite polynomial roots versus exhaustive enumeration; independent good-prime selection; explicit reciprocal mod-two polynomial enumeration. Fixed-space dimension is bounded by algebraic multiplicity, never equated to it.',
        'bindings':branch.bindings([Path(__file__),run.OUTPUT,branch.INPUT,Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['finite_and_infinity_fibres_verified'],'root fibres; global pool dimension in',result['global_pool_dimension_interval'])

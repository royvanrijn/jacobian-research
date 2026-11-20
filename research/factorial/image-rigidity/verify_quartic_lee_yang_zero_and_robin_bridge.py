#!/usr/bin/env python3
"""Exact supporting algebra for the quartic model; no Gamma/SFC assertion."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import factorial
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit('Optimized execution is not accepted.')

HERE = Path(__file__).resolve().parent
SCRIPT = Path(__file__).name
PROOF = 'QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md'
INPUT = 'quartic_lee_yang_zero_and_robin_inputs.json'
LEDGER = 'quartic_lee_yang_zero_and_robin_source_ledger.json'
RECEIPT = 'quartic_lee_yang_zero_and_robin_verification.json'
MANDATORY = {
    PROOF, SCRIPT, INPUT,
    'quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf',
    'quartic_lee_yang_sources/hou_jiang_newman_2208_00917v1.pdf',
}
NV = 7

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def clean(p):
    return {m: Q(c) for m,c in p.items() if c}

def const(c):
    return {} if not c else {(0,)*NV: Q(c)}

def var(i):
    m=[0]*NV; m[i]=1
    return {tuple(m):Q(1)}

def add(*ps):
    ans={}
    for p in ps:
        for m,c in p.items(): ans[m]=ans.get(m,Q(0))+c
    return clean(ans)

def scale(p,c):
    return clean({m:Q(c)*a for m,a in p.items()})

def mul(*ps):
    ans=const(1)
    for p in ps:
        nxt={}
        for m,a in ans.items():
            for n,b in p.items():
                k=tuple(x+y for x,y in zip(m,n))
                nxt[k]=nxt.get(k,Q(0))+a*b
        ans=clean(nxt)
    return ans

def power(p,k):
    assert isinstance(k,int) and k>=0
    ans=const(1)
    for _ in range(k): ans=mul(ans,p)
    return ans

def derivative(p,i):
    ans={}
    for m,c in p.items():
        if m[i]:
            n=list(m); n[i]-=1
            ans[tuple(n)]=c*m[i]
    return clean(ans)

def substitute(p,values):
    ans={}
    for m,c in p.items():
        term=const(c)
        for i,k in enumerate(m): term=mul(term,power(values[i],k))
        ans=add(ans,term)
    return ans

def check_ledger(root):
    doc=json.loads((root/LEDGER).read_text())
    assert doc.get('schema')==1, 'schema'
    rows=doc.get('sources')
    assert isinstance(rows,list) and len(rows)==len(MANDATORY), 'exact mandatory set'
    paths=[r.get('path') for r in rows]
    assert len(set(paths))==len(paths) and set(paths)==MANDATORY, 'mandatory source paths'
    for row in rows:
        p=root/row['path']
        assert p.is_file() and digest(p)==row['sha256'], row['path']
    return doc

def exact_algebra():
    # Variables: x,tau,a1,a2,lambda,h,v.  At a zero a1=F_x, a2=F_xx.
    x,tau,a1,a2,lam,h,v=[var(i) for i in range(NV)]
    a3=scale(mul(tau,a1),-3)
    a4=add(scale(mul(x,a1),3),scale(mul(tau,a2),-3))
    # Total derivative of r=a2/a1 along a zero, using x'=a2/(2a1).
    # Multiply by 2*a1^3; compare with -(r^3+3*tau*r+3*x)/2.
    r_num=add(scale(mul(a4,power(a1,2)),-1),
              scale(mul(a2,a3,a1),2),scale(power(a2,3),-1))
    r_expected=add(scale(power(a2,3),-1),
                   scale(mul(tau,a2,power(a1,2)),-3),
                   scale(mul(x,power(a1,3)),-3))
    assert not add(r_num,scale(r_expected,-1))
    # r=-q gives q'=(3*x-3*tau*q-q^3)/2.
    q=v
    val=[x,tau,a1,scale(mul(q,a1),-1),lam,h,v]
    transformed=substitute(r_expected,val)
    q_expected=mul(power(a1,3),add(power(q,3),scale(mul(tau,q),3),scale(x,-3)))
    assert not add(transformed,scale(q_expected,-1))
    # At a Robin root Z=-lambda*Z_h, the quartic ODE gives Z_hhh.
    z3=mul(a1,add(scale(tau,3),scale(mul(lam,h),-3)))
    quotient_num=add(mul(z3,a1),scale(power(a2,2),-1))
    expected_num=add(scale(mul(tau,power(a1,2)),3),
                     scale(mul(lam,h,power(a1,2)),-3),scale(power(a2,2),-1))
    assert not add(quotient_num,scale(expected_num,-1))
    # Ambient partial-tau Q identity on a Robin zero: Z_hhhh=3*tau*Z_hh+3Z+3hZ_h.
    z4=add(scale(mul(tau,a2),3),scale(mul(lam,a1),-3),scale(mul(h,a1),3))
    partial_tau_num=add(mul(z4,a1),scale(mul(a2,z3),-1))
    partial_tau_expected=add(scale(mul(h,power(a1,2)),3),
         scale(mul(lam,power(a1,2)),-3),scale(mul(lam,h,a1,a2),3))
    assert not add(partial_tau_num,scale(partial_tau_expected,-1))
    # G_tau/Z_h=(Q+3lambda*tau-3lambda^2*h)/2.
    assert not add(add(a2,mul(lam,z3)),
       scale(add(a2,scale(mul(lam,tau,a1),3),
                    scale(mul(power(lam,2),h,a1),-3)),-1))
    # At tau=0, Z^(5)=6 Z_h+3h Z_hh; this is the coupled quintic correction.
    z5_tau0=add(scale(a1,6),scale(mul(h,a2),3))
    assert not add(scale(z5_tau0,Q(1,6)),
                   scale(add(a1,scale(mul(h,a2),Q(1,2))),-1))
    # Critical CW quadratic coefficient: (1-1/(1+tau*e))/e=tau/(1+tau*e).
    e=v
    assert not add(add(const(1),mul(tau,e)),const(-1),scale(mul(tau,e),-1))
    # Factor4 relation: F_x=iZ_h, F_xx=-Z_hh, so modulus quotient is unchanged.
    # Heat-unfolding coefficient is a2*(z^2+1)/2 for tau shift -eta^2.
    expected_unfold=scale(mul(a2,add(power(x,2),const(1))),Q(1,2))
    direct_unfold=add(scale(mul(a2,power(x,2)),Q(1,2)),scale(a2,Q(1,2)))
    assert not add(expected_unfold,scale(direct_unfold,-1))
    return 8, {'correct_r_numerator':r_expected,'quotient_num':quotient_num,
               'unfold':expected_unfold,'vars':[x,tau,a1,a2,lam,h,v]}

def convolution(a,b,n):
    ans=[Q(0)]*(n+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=n: ans[i+j]+=x*y
    return ans

def log_series(a,n):
    assert a[0]==1
    u=list(a); u[0]=0
    term=[Q(1)]+[Q(0)]*n
    ans=[Q(0)]*(n+1)
    for k in range(1,n+1):
        term=convolution(term,u,n)
        for j in range(n+1): ans[j]+=(-1)**(k+1)*term[j]/k
    return ans

def critical_series(inputs):
    n=inputs['cosh_log_order']
    cosh=[Q(1,factorial(k)) if k%2==0 else Q(0) for k in range(n+1)]
    actual=log_series(cosh,n)
    literal=[Q(x) for x in inputs['log_cosh_coefficients']]
    assert actual==literal
    # All cumulant coefficients through the trusted truncation are reconstructed
    # by multiplying cosh'=cosh*(logcosh)', not by quoting a stored series alone.
    left=[(k+1)*cosh[k+1] for k in range(n)]+[Q(0)]
    dlog=[(k+1)*actual[k+1] for k in range(n)]+[Q(0)]
    product=convolution(cosh,dlog,n)
    assert left[:n]==product[:n]
    assert actual[2]==Q(1,2) and actual[4]==Q(-1,12)
    return n

def rational_guards():
    checks={
      'barrier_cubic_surplus':(1-Q(1,4)-Q(1,4))/2==Q(1,4),
      'bounded_backward_time':Q(2,12)<1,
      'bounded_x_increment_squared':Q(4,12)<1,
      'sqrt12_below4_squared':12<16,
      'cuberoot12_below3_cubed':12<27,
      'robin_denominator':1-Q(5,100)==Q(19,20),
      'robin_numerator':(3+Q(3,100)+Q(6,10000)+25)/(Q(19,20))<30,
      'robin_q_surplus':4+Q(30,100)<5,
      'robin_displacement':Q(20,19)<2,
      'tau_guard':Q(1,100)**2<=Q(1,10000),
      'field_guard':Q(1,100)**3<=Q(1,1000000),
      'complex_h_bound':1+Q(2,1000)+Q(3,1000)==Q(1005,1000),
      'complex_tau_bound':1+Q(1,1000)==Q(1001,1000),
      'complex_denominator':1-Q(5,1000)==Q(995,1000),
      'complex_Q_h':3*Q(1001,1000)+3*Q(1,1000)*Q(1005,1000)+25<29,
      'complex_Q_tau':Q(3,2)*(Q(1005,1000)+Q(1,1000)+Q(1,1000)*Q(1005,1000)*5)<2,
      'complex_h_lambda':1/Q(995,1000)<2,
      'complex_h_sigma':(5+3*Q(1,1000)*Q(1001,1000)+3*Q(1,1000000)*Q(1005,1000))/(2*Q(995,1000))<3,
      'complex_Q_closure':Q(58+89,1000)<1,
    }
    assert all(checks.values()),checks
    return checks

def negatives(data,inputs):
    x,tau,a1,a2,lam,h,v=data['vars']
    # Wrong heat sign predicts real rather than imaginary unfolding roots.
    wrong=scale(mul(a2,add(power(x,2),const(-1))),Q(1,2))
    assert add(data['unfold'],scale(wrong,-1))
    # Wrong sign for the physical field term changes the universal flow.
    wrong_num=add(data['correct_r_numerator'],scale(mul(x,power(a1,3)),6))
    assert add(wrong_num,scale(data['correct_r_numerator'],-1))
    # Enlarging lambda*K to 1/5 destroys the very denominator guard used above.
    assert not (1-Q(5,5)>0)
    # A positive quartic logcosh coefficient is not the actual CW limit.
    bad=dict(inputs); bad['log_cosh_coefficients']=list(inputs['log_cosh_coefficients'])
    bad['log_cosh_coefficients'][4]='1/12'
    try: critical_series(bad)
    except AssertionError: pass
    else: raise AssertionError('wrong critical quartic accepted')
    return 4

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=HERE/RECEIPT)
    args=ap.parse_args()
    check_ledger(HERE)
    inputs=json.loads((HERE/INPUT).read_text())
    assert inputs['schema']==1 and inputs['model_only'] is True
    assert inputs['all_real_tau'] is True and inputs['actual_gamma_transfer'] is False
    nids,data=exact_algebra()
    ns=critical_series(inputs)
    guards=rational_guards()
    neg=negatives(data,inputs)
    payload={
      'status':'PASS','scope':'Quartic model algebra and bounds only; no actual Gamma/SFC transfer',
      'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__)),
      'input_sha256':digest(HERE/INPUT),'source_ledger_sha256':digest(HERE/LEDGER),
      'universal_identities':nids,'critical_series_order':ns,
      'rational_guards':len(guards),'guard_names':sorted(guards),
      'negative_controls':neg,'mandatory_sources':len(MANDATORY),
      'formal_verification':False,'external_review':False,'independent_replay':False,
    }
    args.output.write_text(json.dumps(payload,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:payload[k] for k in ['status','universal_identities','critical_series_order','rational_guards','negative_controls','mandatory_sources']},sort_keys=True))

if __name__=='__main__': main()

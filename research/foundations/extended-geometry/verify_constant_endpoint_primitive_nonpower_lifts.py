#!/usr/bin/env python3
"""Literal original-variable controls for primitive nonpower Phi=1 lifts."""
if not __debug__:raise RuntimeError('Refusing -O: exact assertions are required')
from fractions import Fraction as Q
from pathlib import Path
from math import comb,factorial,gcd
from functools import reduce
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ZERO=(0,0,0,0)
ONE={ZERO:Q(1)}

def add(p,q,scale=1):
    out=dict(p)
    for m,v in q.items():out[m]=out.get(m,0)+scale*v
    return {m:v for m,v in out.items()if v}
def mul(p,q):
    out={}
    for m,v in p.items():
        for n,w in q.items():
            k=tuple(a+b for a,b in zip(m,n));out[k]=out.get(k,0)+v*w
    return {m:v for m,v in out.items()if v}
def power(p,n):
    out=ONE
    for _ in range(n):out=mul(out,p)
    return out
def divide(p,q):
    """Exact multivariate division in lex order, with full remainder."""
    work=dict(p);out={};rem={};lm=max(q);lc=q[lm]
    while work:
        m=max(work);c=work[m]
        if all(a>=b for a,b in zip(m,lm)):
            k=tuple(a-b for a,b in zip(m,lm));v=c/lc
            out[k]=out.get(k,0)+v
            work=add(work,{tuple(a+b for a,b in zip(k,n)):v*w for n,w in q.items()},-1)
        else:rem[m]=work.pop(m)
    return out,rem
def order(p,q):
    n=0
    while True:
        quotient,rem=divide(p,q)
        if rem:return n
        assert quotient and n<100
        n+=1;p=quotient

def peval(p,point):return sum((v*reduce(lambda a,b:a*b,(x**k for x,k in zip(point,m)),1)for m,v in p.items()),Q(0))
def contraction(p):return sum((v*factorial(i)*factorial(j)for(i,j,k,l),v in p.items()if i==k and j==l),Q(0))
def umul(p,q):
    out=[Q(0)]*(len(p)+len(q)-1)
    for i,v in enumerate(p):
        for j,w in enumerate(q):out[i+j]+=v*w
    return tuple(out)
def upower(p,n):
    out=(Q(1),)
    for _ in range(n):out=umul(out,p)
    return out

def raw(a,m,q,C,D):
    if a*m%2:return Q(0)
    k=m-q
    if k<0:return Q(0)
    out=Q(0)
    for j in range(m+1):
        coeff=umul(upower(C,m-j),upower(D,j))
        if k<len(coeff):out+=Q(comb(m,j),a*m+2*j+1)*coeff[k]
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output','--write-receipt',dest='output',type=Path)
    args=parser.parse_args()
    xiz={(1,0,1,0):Q(1)};xi2z2={(0,1,0,1):Q(1)}
    R=add(xiz,xi2z2);T=add(xiz,xi2z2,-1);Z={(1,0,0,1):Q(1)};W={(0,1,1,0):Q(2)}
    assert add(add(power(T,2),mul(Z,W),2),power(R,2),-1)=={}
    def homogeneous(p):
        out={};d=len(p)-1
        for j,v in enumerate(p):out=add(out,mul(power(R,d-j),power(Z,j)),v)
        return out
    inputs_path=HERE/'constant_endpoint_primitive_nonpower_lifts_inputs.json'
    inputs=json.loads(inputs_path.read_text());assert len(inputs['fixtures'])==3
    fixtures=[(row['a'],tuple(map(Q,row['P'])),tuple(map(Q,row['H'])),
               [(Q(alpha),int(order))for alpha,order in row['geometric_roots_and_H_orders']],
               row['balanced_degree'])for row in inputs['fixtures']]
    results=[];negative=False
    for a,P,H,roots,expected_degree in fixtures:
        r=len(P)-1;h=len(H)-1;d=(a+2)*(r+1)
        assert h==a*r+1 and P[0]==H[0]==1 and d==expected_degree
        Ph,Hh=homogeneous(P),homogeneous(H)
        quotient,remainder=divide(add(power(Ph,2),power(R,2*r),-1),Z)
        assert not remainder
        B=add(mul(W,power(R,2*r)),mul(power(T,2),quotient),Q(-1,2))
        B={m:2*v for m,v in B.items()}
        assert mul(Z,B)==add(power(R,2*r+2),mul(power(T,2),power(Ph,2)),-1)
        F=mul(mul(power(T,a),Hh),B)
        assert F and all(sum(m[:2])==sum(m[2:])==d for m in F)
        torder=order(F,T);rorder=order(F,R)
        assert torder==a and rorder==0
        horders=[]
        for alpha,expected in roots:
            alpha=Q(alpha);L=add(R,Z,-1/alpha)
            observed=order(F,L);assert observed==expected
            assert divide(mul(Z,B),L)[1]==divide(power(R,2*r+2),L)[1]
            assert order(B,L)==0
            horders.append(observed)
        assert sum(horders)==h and reduce(gcd,[torder]+horders)==1
        point=(Q(1),Q(1),Q(1),Q(-1));value=peval(F,point)
        expected=-H[-1]*P[-1]**2*Q(2)**(a+2)*Q(-1)**(h+2*r-1)
        assert peval(R,point)==0 and value==expected!=0
        # Known proper-power negative control: actual factor orders double.
        squared=power(F,2)
        squared_orders=[order(squared,T)]+[order(squared,add(R,Z,-1/Q(alpha)))for alpha,_ in roots]
        assert squared_orders==[2*torder]+[2*v for v in horders]
        assert reduce(gcd,squared_orders)==2
        D=tuple(-v for v in umul(H,upower(P,2)))
        contractions=[]
        for m in (1,2):
            Fm=power(F,m);pure=contraction(Fm);marked=contraction(mul(Z,Fm))
            assert pure==factorial(d*m+1)*raw(a,m,0,H,D)==0
            target=raw(a,m,1,H,D)
            assert marked==factorial(d*m+2)*target
            if m==2:assert marked!=0
            contractions.append({'m':m,'normalized_pure':str(pure/Q(factorial(d*m+1))),
                                 'normalized_marked':str(target),'expanded_power_terms':len(Fm)})
        try:assert torder==a+1
        except AssertionError:negative=True
        else:raise AssertionError('Corrupt factor order accepted')
        results.append({'a':a,'r':r,'balanced_degree':d,'expanded_terms':len(F),
                        'T_order':torder,'geometric_L_orders':horders,'R_order':rorder,
                        'R_zero_specialization':str(value),'proper_square_order_gcd':2,
                        'original_contractions':contractions})
    ledger_path=HERE/'constant_endpoint_primitive_nonpower_lifts_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    pins=[]
    for entry in ledger['retained_resources']:
        path=ROOT/entry['retained_path'];got=hashlib.sha256(path.read_bytes()).hexdigest()
        assert got==entry['sha256'];pins.append({'path':entry['retained_path'],'sha256':got})
    out={'passed':True,'scope':'Three literal original-variable fixtures and exact factor-order checks supplement the all-extension proof; no purity/nonpower conclusion is inferred universally from bounded tests.',
         'fixtures':results,'changed_factor_order_rejected':negative,'retained_source_pins':pins,
         'packet_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in
             [HERE/'CONSTANT_ENDPOINT_PRIMITIVE_NONPOWER_LIFTS.md',Path(__file__),inputs_path,ledger_path]}}
    output=args.output or HERE/'constant_endpoint_primitive_nonpower_lifts_verification.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'literal_fixtures':len(results),'original_powers':6,'output':str(output)}))
if __name__=='__main__':main()

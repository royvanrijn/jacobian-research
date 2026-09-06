#!/usr/bin/env python3
"""CAS-free Rabin checks for a three-prime rational irreducibility certificate."""
import argparse
from pathlib import Path
import retrospective as r
import surface_discriminant_irreducibility as source

OUTPUT=r.OUT/'rank_jump_surface_discriminant_modular_verification_v1.json'


def trim(a):
    a=list(a)
    while len(a)>1 and not a[-1]:a.pop()
    return a


def sub(a,b,p):
    return trim([((a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0))%p for i in range(max(len(a),len(b)))])


def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return trim(c)


def mod(a,b,p):
    a=trim(a);b=trim(b);inv=pow(b[-1],-1,p)
    while a!=[0] and len(a)>=len(b):
        shift=len(a)-len(b);c=a[-1]*inv%p
        for i,v in enumerate(b):a[shift+i]=(a[shift+i]-c*v)%p
        a=trim(a)
    return a


def gcd(a,b,p):
    while b!=[0]:a,b=b,mod(a,b,p)
    return [(x*pow(a[-1],-1,p))%p for x in a]


def power(a,n,f,p):
    out=[1]
    while n:
        if n&1:out=mod(mul(out,a,p),f,p)
        a=mod(mul(a,a,p),f,p);n//=2
    return out


def irreducible(f,p):
    n=len(f)-1;x=[0,1];xp=x;frobenius=[x]
    for _ in range(n):
        xp=power(xp,p,f,p);frobenius.append(xp)
    if mod(sub(xp,x,p),f,p)!=[0]:return False
    divisors=[q for q in range(2,n+1) if n%q==0 and all(q%d for d in range(2,q))]
    return all(gcd(f,sub(frobenius[n//q],x,p),p)==[1] for q in divisors)


def compute():
    data=r.read(source.OUTPUT);assert data['status']=='PASS_IRREDUCIBLE'
    for name,digest in data['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
    coeff=list(map(int,data['primitive_coefficients_ascending']))
    candidates=set(range(1,24));rows=[]
    for p in [167,181,191]:
        row=next(x for x in data['modular_certificate'] if x['prime']==p)
        factors=row['monic_factor_coefficients_ascending'];product=[1]
        assert len({tuple(f) for f in factors})==len(factors)
        for f in factors:assert irreducible(f,p);product=mul(product,f,p)
        assert product==[(c*pow(coeff[-1],-1,p))%p for c in coeff]
        degrees=[len(f)-1 for f in factors];sums={0}
        for d in degrees:sums |= {s+d for s in list(sums)}
        candidates &= sums-{0,24}
        rows.append({'prime':p,'factor_degrees':degrees,'irreducible_factors_verified':True,
            'product_verified':True,'remaining_proper_factor_degrees':sorted(candidates)})
    assert not candidates
    return {'schema':'rank-jump.surface-discriminant-modular-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),source.OUTPUT)},
        'status':'PASS_IRREDUCIBLE_OVER_Q','rows':rows,
        'method':'Exact finite-polynomial multiplication and Rabin irreducibility; Gauss lemma excludes every proper rational factor degree.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=compute()
    if a.mode=='check':assert r.read(OUTPUT)==data;print('PASS CAS-free degree-24 irreducibility')
    else:r.write_new(OUTPUT,data)

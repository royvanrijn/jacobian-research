#!/usr/bin/env python3
"""Independent finite fields and exhaustive square-set point counts, no CAS."""
import argparse
from itertools import product
from pathlib import Path
import retrospective as r
import auxiliary_elliptic_multiplicity as initial
import auxiliary_multiplicity_completion as completed

OUTPUT=r.OUT/'rank_jump_auxiliary_multiplicity_verification_v1.json'


class Field:
    def __init__(self,p,n):
        self.p=p;self.n=n;self.q=p**n
        self.digits=[tuple((a//p**i)%p for i in range(n)) for a in range(self.q)]
        # A degree-two or degree-three polynomial is irreducible iff it has no base-field root.
        self.modulus=None if n==1 else next(list(c)+[1] for c in product(range(p),repeat=n)
            if c[0] and all((sum(c[i]*a**i for i in range(n))+a**n)%p for a in range(p)))

    def add_constant(self,x,c):
        return x-x%self.p+(x%self.p+c)%self.p

    def mul(self,x,y):
        p,n=self.p,self.n
        if n==1:return x*y%p
        a=self.digits[x];b=self.digits[y];out=[0]*(2*n-1)
        for i in range(n):
            for j in range(n):out[i+j]+=a[i]*b[j]
        for k in range(2*n-2,n-1,-1):
            for j in range(n):out[k-n+j]-=out[k]*self.modulus[j]
        return sum((out[i]%p)*p**i for i in range(n))

    def evaluate(self,coeff,x):
        value=0
        for c in reversed(coeff):value=self.add_constant(self.mul(value,x),c%self.p)
        return value


def count(coeff,p,n):
    F=Field(p,n);squares={F.mul(x,x) for x in range(F.q)}
    zeros=0;nonzero=0
    for x in range(F.q):
        value=F.evaluate(coeff,x)
        if not value:zeros+=1
        elif value in squares:nonzero+=1
    infinity=2*int(coeff[-1]%p in squares)
    return {'extension_degree':n,'field_size':F.q,'affine_branch_points':zeros,
        'nonzero_square_values':nonzero,'points_at_infinity':infinity,'points':zeros+2*nonzero+infinity},F.modulus


def remainder(poly,divisor):
    out=list(map(r.F,poly))
    while len(out)>=len(divisor):
        a=out[-1]/divisor[-1];k=len(out)-len(divisor)
        for i,b in enumerate(divisor):out[k+i]-=a*b
        while out and out[-1]==0:out.pop()
    return out or [r.F(0)]


def compute():
    original=r.read(initial.OUTPUT);completion=r.read(completed.OUTPUT)
    for name,digest in completion['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
    rows=[]
    for record in [original['rows'][0],completion['row']]:
        assert record['status']=='BOTH_Q_HOM_SPACES_ZERO'
        for name,digest in record['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
        decisive=next(x for x in record['records'] if x['status']=='PASS' and all(y['factor_excluded'] for y in x['elliptic_checks']))
        p=decisive['prime'];coeff=list(map(int,record['product_polynomial_ascending']))
        counts=[];moduli=[]
        for n in [1,2,3]:
            value,modulus=count(coeff,p,n);assert value==decisive['counts'][n-1]
            counts.append(value);moduli.append(modulus)
        s1,s2,s3=[p**n+1-counts[n-1]['points'] for n in [1,2,3]]
        assert (s1*s1-s2)%2==0 and (s1**3-3*s1*s2+2*s3)%6==0
        e2=(s1*s1-s2)//2;e3=(s1**3-3*s1*s2+2*s3)//6
        W=[p**3,-p*p*s1,p*e2,-e3,e2,-s1,1]
        assert list(map(str,W))==decisive['Jacobian_Frobenius_coefficients_ascending']
        for row in decisive['elliptic_checks']:
            _,a,_,b,c=initial.MODELS[row['name']]
            N=1+sum((y*y-(x**3+a*x*x+b*x+c))%p==0 for x,y in product(range(p),repeat=2))
            assert N==row['point_count']
            P=[p,N-p-1,1]
            rem=remainder(W,P)
            assert list(map(str,rem))==row['remainder_coefficients_ascending'] and any(rem)
        rows.append({'alignment':record['alignment'],'prime':p,'status':'PASS',
            'independent_field_moduli_ascending':moduli,'independent_point_counts':[x['points'] for x in counts],
            'Frobenius_coefficients_ascending':W,'both_rational_elliptic_factors_excluded':True,
            'signed_genus_five_rational_elliptic_multiplicity':2})
    # Distinct traces at 11 exclude a Q-isogeny between the two elliptic controls.
    assert original['rows'][0]['records'][-1]['elliptic_checks'][0]['trace']!=original['rows'][0]['records'][-1]['elliptic_checks'][1]['trace']
    return {'schema':'rank-jump.auxiliary-multiplicity-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),initial.OUTPUT,completed.OUTPUT)},
        'rows':rows,'control_curves_not_Q_isogenous':True,
        'boundary':'Independent finite-field arithmetic. Multiplicity uses the retained V4 Jacobian decomposition and End_Q(E)=Z; it is not a Mordell-Weil rank or a geometric-Hom claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=compute()
    if a.mode=='check':assert r.read(OUTPUT)==data;print('PASS independent finite-field multiplicity verification')
    else:r.write_new(OUTPUT,data)

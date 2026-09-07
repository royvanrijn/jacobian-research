#!/usr/bin/env python3
"""Independent full fibre recount and exact Frobenius moment verification."""
import argparse
from fractions import Fraction as F
from math import prod
from pathlib import Path
import numpy as np
from sage.all import GF, PolynomialRing, EllipticCurve
import retrospective as r
from native_twist_frobenius import geometry

HERE=Path(__file__).resolve().parent
SOURCE=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
PARITY=r.OUT/'rank_jump_native_twist_moment_parity_v1.json'
OUTPUT=r.OUT/'rank_jump_native_twist_frobenius_verification_v2.json'


def verify():
    data=r.read(SOURCE);parity=r.read(PARITY)
    for old in (data,parity):
        for path,sha in old['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    assert geometry()==data['geometry']
    p=131;Rp=PolynomialRing(GF(p),'z');z=Rp.gen()
    k=GF(p*p,'alpha',modulus=z*z-2);alpha=k.gen()
    R=PolynomialRing(k,'t');A,B,q=[R(c) for c in data['geometry']['modular_coefficients']]
    def coords(x):
        c=list(x);return (int(c[0]) if c else 0,int(c[1]) if len(c)>1 else 0)
    # Build the character table by listing actual squares in Sage's finite field.
    characters=np.full(p*p,-1,dtype=np.int64);characters[0]=0
    for a in range(p):
        for b in range(p):
            if a or b:
                u,v=coords((k(a)+b*alpha)**2);characters[u+p*v]=1
    xr=np.tile(np.arange(p,dtype=np.int64),p);xi=np.repeat(np.arange(p,dtype=np.int64),p)
    # Expand the cubic directly; the C++ producer uses repeated field multiplication.
    cube_r=(xr*xr*xr+6*xr*xi*xi)%p;cube_i=(3*xr*xr*xi+2*xi*xi*xi)%p
    totals={n:[0,0] for n in(1,2)};checked=0;sage_checks=0
    for n,a,b,mult,old_trace,old_twist in data['fibre_trace_ledger']:
        if a==-1:aa,bb,qq=A[8],B[12],q[2]
        else:
            t=k(a)+b*alpha;aa,bb,qq=A(t),B(t),q(t)
        ar,ai=coords(aa);br,bi=coords(bb);qr,qi=coords(qq)
        if n==1:
            assert ai==bi==qi==b==0
            field=GF(p);c=[0 if x==0 else (1 if field(x).is_square() else -1) for x in range(p)]
            trace=-sum(c[(x**3+ar*x+br)%p] for x in range(p));sign=c[qr]
        else:
            vr=(cube_r+ar*xr+2*ai*xi+br)%p
            vi=(cube_i+ar*xi+ai*xr+bi)%p
            trace=-int(characters[vr+p*vi].sum());sign=int(characters[qr+p*qi])
        assert trace==old_trace and sign*trace==old_twist
        smooth=4*aa**3+27*bb**2!=0
        if smooth and (n==1 or checked%541==0 or a==-1):
            field=GF(p) if n==1 else k
            E=EllipticCurve(field,[field(ar) if n==1 else aa,field(br) if n==1 else bb])
            assert int(E.trace_of_frobenius())==trace;sage_checks+=1
        if not smooth:assert trace in(-1,1)
        totals[n][0]-=mult*trace;totals[n][1]-=mult*sign*trace;checked+=1
    assert checked==132+8647
    for j,row in enumerate(data['rows']):
        tr=[totals[n][j] for n in(1,2)];assert tr==row['Frobenius_traces']
        N=row['cohomology_dimension'];s=F(tr[0],p);ss=F(tr[1],p*p);M=(N+ss)/2
        c=F(row['quadratic_center'])
        bound=(M-2*c*s+N*c*c)/(1-c)**2
        assert str(bound)==row['exact_eigenvalue_one_bound'] and bound.numerator//bound.denominator==row['arithmetic_generic_rank_upper_bound']
    # Independently compute local split characters inside each residue field.
    R0=PolynomialRing(GF(p),'t');aa,bb,qq=[R0(c) for c in data['geometry']['modular_coefficients']]
    d=-16*(4*aa**3+27*bb**2);product=R0(1);Ws=[1,1]
    for row in parity['multiplicative_places']:
        f=R0(row['factor']);assert f.is_irreducible();product*=f
        if f.degree()==1:
            root=-f[0];field=GF(p)
        else:
            field=GF(p**f.degree(),'beta',modulus=f);root=field.gen()
        cb=PolynomialRing(field,'s')(bb)(root);cq=PolynomialRing(field,'s')(qq)(root)
        chars=[1 if x.is_square() else -1 for x in (864*cb,864*cb*cq**3)]
        assert chars==row['split_characters']
        for j,chi in enumerate(chars):Ws[j]*=-chi
    assert product==d.monic()
    assert Ws==[row['global_root_number'] for row in parity['rows']]
    for old,new in zip(data['rows'],parity['rows'],strict=True):
        m=old['arithmetic_generic_rank_upper_bound']
        while (-1)**m!=new['global_root_number']:m-=1
        assert m==new['parity_refined_upper_bound']
    return {'schema':'rank-jump.native-twist-frobenius-verification.v1','status':'PASS',
            'fibre_orbits_independently_recounted':checked,'Sage_elliptic_trace_crosschecks':sage_checks,
            'Frobenius_traces':{str(n):v for n,v in totals.items()},'refined_bounds':parity['rows'],
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (SOURCE,PARITY,Path(__file__),HERE/'native_twist_frobenius.py',HERE/'retrospective.py')},
            'boundary':'All counted fibre orbits are independently recounted. The cycle-class specialization, L-function purity/dimension, and local epsilon-factor formulas are mathematical dependencies stated in the companion note.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);mode=p.parse_args().mode
    result=verify()
    if mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS complete independent recount;',result['refined_bounds'],flush=True)

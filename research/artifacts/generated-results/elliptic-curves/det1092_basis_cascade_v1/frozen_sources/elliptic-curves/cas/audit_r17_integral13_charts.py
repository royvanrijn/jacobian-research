#!/usr/bin/env python3
"""Exact determinant13 parameter charts, without a new parameter or point search."""
import argparse,json
from fractions import Fraction as F
from math import comb
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ART/'compact_six_r17_atlas_v1.json';OUT=ART/'r17_integral13_parameter_charts_v1.json'

def reduced_basis(residue):
    u,v=((1,0),(0,13)) if residue=='infinity' else ((13,0),(int(residue),1))
    for _ in range(32):
        norm=lambda a:sum(x*x for x in a)
        if norm(v)<norm(u):u,v=v,u
        mu=round(F(sum(a*b for a,b in zip(u,v)),norm(u)))
        if not mu:break
        v=tuple(a-mu*b for a,b in zip(v,u))
    else:raise ArithmeticError('finite lattice reduction bound exceeded')
    if u[0]<0 or (u[0]==0 and u[1]<0):u=tuple(-x for x in u)
    if u[0]*v[1]-u[1]*v[0]<0:v=tuple(-x for x in v)
    return [u[0],v[0],u[1],v[1]]

def transform(coefficients,degree,matrix):
    a,b,c,d=map(F,matrix);out=[F(0)]*(degree+1)
    for i,k in enumerate(coefficients):
        for j in range(i+1):
            for l in range(degree-i+1):out[j+l]+=k*comb(i,j)*a**j*b**(i-j)*comb(degree-i,l)*c**l*d**(degree-i-l)
    return out

def size(a,b):return max(4*sum(abs(x) for x in a)**3,27*sum(abs(x) for x in b)**2)

def expected():
    rows=[]
    for f in cert.read(INPUT)['families']:
        a=list(map(F,f['A_coefficients_low_to_high']));b=list(map(F,f['B_coefficients_low_to_high']));old=size(a,b)
        for residue in [*range(13),'infinity']:
            m=reduced_basis(residue);x,y,z,w=m
            if x*w-y*z!=13:raise ArithmeticError('index13 lattice basis required')
            if residue=='infinity':
                if z%13 or w%13:raise ArithmeticError('infinity residue lattice differs')
            elif (x-residue*z)%13 or (y-residue*w)%13:raise ArithmeticError('affine residue lattice differs')
            aa=[v/13**4 for v in transform(a,8,m)];bb=[v/13**6 for v in transform(b,12,m)];integral=all(v.denominator==1 for v in aa+bb)
            r={'family':f['family'],'residue_mod13':residue,'matrix':m,'integral_after_curve_scale13':integral,'A_coefficients_low_to_high':list(map(str,aa)),'B_coefficients_low_to_high':list(map(str,bb))}
            if integral:
                aliases=[]
                for swap in (False,True):
                    for sn in (-1,1):
                        for sd in (-1,1):
                            q=[0,sn,sd,0] if swap else [sn,0,0,sd]
                            if transform(a,8,q)==aa and transform(b,12,q)==bb:aliases.append(q)
                r.update(weighted_coefficient_bound_ratio=str(size(aa,bb)/old),signed_permutation_self_presentations=aliases,remaining_common_affine_roots_mod13=[i for i in range(13) if sum(int(v)*i**j for j,v in enumerate(aa))%13==0 and sum(int(v)*i**j for j,v in enumerate(bb))%13==0],remaining_common_infinity_root_mod13=int(aa[-1])%13==0 and int(bb[-1])%13==0)
            rows.append(r)
    if len(rows)!=84:raise ArithmeticError('fixed six families and14 projective residue cells required')
    return {'schema':'elliptic-curves.r17-integral13-parameter-charts.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),INPUT,Path(cert.__file__).resolve()]},'rows':rows,'integral_charts':sum(r['integral_after_curve_scale13'] for r in rows),'scope':'An exact coefficient audit of all14 projective residue cells modulo13 in each of six compact R17 families. A Gauss-reduced integer matrix of determinant13 spans each residue lattice. Substitute the binary forms and divide A by13^4 and B by13^6, retaining every nonintegral outcome. Integral charts preserve the Q(t) family and generic rank by an invertible rational base change and explicit Weierstrass scaling. Primitive old parameter pairs additionally require M(u,v) not both divisible by13. Coefficient bounds and bounded signed-permutation comparisons do not establish optimality, new curve populations, new point visibility or rank gains. No parameter scan, trace extension or point search occurs.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('exact integral13 chart audit differs')
    else:
        if OUT.exists():raise FileExistsError('preserve integral13 chart proof')
        checkpoint(OUT,d)
    print('EXACT84 PARAMETER CELLS;',d['integral_charts'],'INTEGRAL13 CHARTS',flush=True)
    for r in d['rows']:
        if r['integral_after_curve_scale13']:print(r['family'],r['residue_mod13'],r['matrix'],r['weighted_coefficient_bound_ratio'],'self presentations',r['signed_permutation_self_presentations'],flush=True)

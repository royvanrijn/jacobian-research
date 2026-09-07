#!/usr/bin/env python3
"""Sage-free exact verification of the retrospective cover-point witnesses."""
import argparse
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import retrospective as r
SOURCE=r.OUT/'rank_jump_constructed_class_oracle_solubility_v1.json'
CLASSES=r.OUT/'rank_jump_constructed_class_compaction_v1.json'
REFERENCE=r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json'
ORACLE=r.OUT/'small_conductor_rank22_proof_v1.json'
OUTPUT=r.OUT/'rank_jump_constructed_cover_points_verification_v1.json'

def compute():
    source=r.read(SOURCE);classes=r.read(CLASSES);ref=r.read(REFERENCE);oracle=r.read(ORACLE)
    for path,sha in source['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    f=list(map(F,ref['cubic_ascending']));assert f[3]==1
    def mul(a,b):
        v=[sum((a[i]*b[j] for i in range(3) for j in range(3) if i+j==k),F(0)) for k in range(5)]
        for k in [4,3]:
            for j in range(3):v[k-3+j]-=v[k]*f[j]
        return v[:3]
    a1,a2,a3,a4,a6=map(F,oracle['integral_model']);assert (a1,a2,a3)==(1,0,0)
    rows=[]
    for witness,case in zip(source['cases'],classes['cases']):
        assert witness['column']==case['column'] and witness['status']=='RATIONAL_CERTIFIED_RETROSPECTIVELY'
        x,y=map(F,witness['point_on_original_model']);X,Y=map(F,witness['point_on_cubic_model'])
        assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
        assert X==4*x and Y==8*y+4*x and Y*Y==sum(c*X**j for j,c in enumerate(f))
        beta=list(map(F,case['beta_ascending']));xi=list(map(F,witness['cubic_field_square_root_ascending']))
        assert mul(beta,mul(xi,xi))==[X,F(-1),F(0)]
        u,v,w,s=map(int,witness['primitive_cover_point']);assert s and gcd(gcd(u,v),gcd(w,s))==1
        assert [F(u,s),F(v,s),F(w,s)]==xi
        Q=mul(beta,mul(list(map(F,[u,v,w])),list(map(F,[u,v,w]))))
        assert Q[2]==0 and Q[1]+s*s==0 and Q[0]==X*s*s
        rows.append({'column':case['column'],'rational_cover_point_verified':True,'exact_square_identity_verified':True})
    return {'schema':'rank-jump.constructed-cover-points-verification.v1','status':'PASS','cases':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),SOURCE,CLASSES,REFERENCE,ORACLE]},
        'boundary':'Independent rational arithmetic verifies the displayed point and cubic-field identities. Oracle-assisted retrospective solubility only; no prospective constructor consumes this output.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS two rational cover points, using exact fractions only')

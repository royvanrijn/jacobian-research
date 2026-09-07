#!/usr/bin/env python3
"""Independent norm determinants and private-good-prime obstruction proof."""
import argparse
from math import gcd,isqrt
from pathlib import Path
from sage.all import QQ,ZZ,matrix
import retrospective as r

INPUT=r.OUT/'rank_jump_generic_chord_ramification_v1.json'
OUTPUT=r.OUT/'rank_jump_generic_chord_ramification_verification_v1.json'


def verify():
    source=r.read(INPUT)
    for path,d in source['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==d
    panel=r.read(r.ROOT/'artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_panel_inputs_v1.json')
    rows=[]
    for row in source['rows']:
        assert row['status']=='PASS'
        A,B=int(row['A']),int(row['B']);S=row['S_finite'];disc=-4*A**3-27*B**2
        original=next(x for x in panel['cases'] if x['token']==row['token'])
        a1,a2,a3,a4,a6=map(QQ,original['model'])
        b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
        AA=b4/2-b2*b2/48;BB=b6/4-b2*b4/24+b2**3/864
        scale=AA.denominator().lcm(BB.denominator())
        assert AA*scale**4==A and BB*scale**6==B
        points=[((QQ(x)+b2/12)*scale**2,(QQ(y)+(a1*QQ(x)+a3)/2)*scale**3)
                for x,y in original['generic_sections']]
        assert all(y*y==x**3+A*x+B for x,y in points)
        expected={(i,i,1) for i in range(len(points))}
        expected|={(i,j,s) for i in range(len(points)) for j in range(i+1,len(points)) for s in (-1,1)}
        encountered=set();norms=[]
        for item in row['forms']:
            a,b=int(item['a']),int(item['b']);assert b>0 and gcd(a,b)==1
            determinant=matrix(ZZ,[[a,0,-b*B],[b,a,-b*A],[0,b,a]]).det()
            assert str(determinant)==item['norm'];norms.append(abs(int(determinant)))
            for i,j,sign in item['sources']:
                encountered.add((i,j,sign));x,y=points[i];u,v=points[j];v*=sign
                if i==j:
                    assert (a+b*x)*(3*x*x+A)==2*b*y*y
                else:
                    assert (a+b*x)*v==(a+b*u)*y
                    assert (v-y)!=0
        assert encountered==expected
        assert len(norms)==len(set((x['a'],x['b']) for x in row['forms']))
        atoms=list(map(int,row['atoms']));witnesses=[]
        # This proof does not use the producer's root partitions, parity
        # signatures, Gaussian elimination, or gcd-refinement algorithm.
        for i,item in enumerate(row['forms']):
            candidates=[]
            for j,c in enumerate(atoms):
                if isqrt(c)**2==c or gcd(c,disc)!=1 or any(c%p==0 or gcd(c,p)>1 for p in S):continue
                n=norms[i];e=0
                while n%c==0:n//=c;e+=1
                if e%2!=1 or gcd(n,c)!=1:continue
                if any(gcd(c,norms[k])!=1 for k in range(len(norms)) if k!=i):continue
                a,b=int(item['a']),int(item['b']);assert gcd(b,c)==1
                root=(-a*pow(b,-1,c))%c
                assert (root**3+A*root+B)%c==0 and gcd(3*root*root+A,c)==1
                candidates.append((c.bit_length(),j,e))
            assert candidates,'no independent private-atom certificate'
            bits,j,e=min(candidates)
            witnesses.append({'form_index':i,'atom_index':j,'norm_exponent':e,'atom_bits':bits})
        assert row['outside_S_parity_rank']==len(norms) and row['coefficient_kernel_dimension']==0
        rows.append({'token':row['token'],'status':'PASS','forms_verified':len(norms),
            'private_obstructions':witnesses,'outside_S_parity_rank':len(norms),
            'Selmer_dimension_beyond_generic_in_this_constructor':0})
    return {'schema':'rank-jump.generic-chord-ramification-verification.v1','status':'PASS',
            'method':'independent cubic multiplication determinants, source chord identities, and private coprime nonsquare atoms; no parity elimination or integer factorization',
            'rows':rows,'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),INPUT]}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','check']);args=p.parse_args()
    result=verify()
    if args.mode=='capture':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],[x['forms_verified'] for x in result['rows']])

#!/usr/bin/env python3
"""Pure-rational and polynomial-arithmetic separation replay."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import governing_third_class_separation as run
from verify_explicit_governing_octic import finite_polynomials

OUTPUT=r.OUT/'rank_jump_governing_third_class_separation_verification_v1.json'


def compute():
    out=r.read(run.OUTPUT);octics={x['token']:x for x in r.read(run.octics.OUTPUT)['rows']}
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    rows=[]
    for row in out['rows']:
        assert row['status']=='PASS';w=row['witness'];p=w['prime']
        assert p in r.primes(10009)
        mod=lambda z:F(z).numerator*pow(F(z).denominator,-1,p)%p
        ev=lambda cs,x:sum(c*pow(x,i,p) for i,c in enumerate(cs))%p
        x,y=map(F,row['third_generic_point']);f=list(map(F,row['cubic_ascending']));A=f[1]
        assert y*y==sum(c*x**i for i,c in enumerate(f))
        q=list(map(F,row['third_class_quartic_ascending']))
        assert q==[-3*x*x/16-A/4,-y,-3*x/2,0,1]
        hp=list(map(mod,octics[row['token']]['integral_octic_ascending']));qp=list(map(mod,q));fp=list(map(mod,f))
        assert len(hp)==9 and hp[-1]!=0 and len(set(w['octic_roots']))==8
        assert all(0<=z<p and ev(hp,z)==0 for z in w['octic_roots'])
        trim,sub,mul,power,gcd=finite_polynomials(p)
        assert len(gcd(qp,[i*qp[i] for i in range(1,5)]))==1
        assert all(ev(qp,z)!=0 for z in range(p))
        assert sub(power([0,1],p*p,qp),[0,1])==[]
        roots=w['cubic_roots'];assert len(set(roots))==3 and all(ev(fp,z)==0 for z in roots)
        sig=[int(pow((mod(x)-z)%p,(p-1)//2,p)==p-1) for z in roots]
        assert all((mod(x)-z)%p for z in roots) and sig==w['third_kummer_signature'] and sum(sig)==2
        rows.append({'token':row['token'],'prime':p,'octic_splits_completely':True,'third_quartic_splitting':[2,2],'Kummer_signature':sig})
    files=(Path(__file__),run.OUTPUT,run.octics.OUTPUT,Path(r.__file__),Path(__file__).with_name('verify_explicit_governing_octic.py'))
    return {'schema':'rank-jump.governing-third-class-separation-verification.v1','status':'PASS','rows':rows,
        'method':'Eight distinct roots certify full octic splitting; no quartic root and Frobenius^2=identity certify two quadratics, without a factorization routine.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS three independent finite-prime field-separation witnesses')

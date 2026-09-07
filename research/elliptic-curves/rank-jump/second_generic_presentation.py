#!/usr/bin/env python3
"""Finite equation-only j-map inversion on retained curves; no point search."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import generic_selmer_capacity as prior

PROTOCOL=Path(__file__).with_name('SECOND_GENERIC_PRESENTATION_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_second_generic_presentation_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_second_generic_presentation_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-second-generic-presentation-v1'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    data=r.read(prior.INPUT);cases=[]
    for c in data['cases']:
        model,_=r.short(c['model'],[])
        assert model[:3]==['0','0','0']
        cases.append({'token':c['token'],'a':model[3],'b':model[4],
            'control_family':c['family'],'control_parameter':c['parameter']})
    assert len(cases)==16 and len(data['families'])==7
    r.write_new(INPUT,{'schema':'rank-jump.second-generic-presentation-inputs.v1','cases':cases,
        'families':[{k:f[k] for k in ('family','A','B')} for f in data['families']],
        'bindings':bindings([Path(__file__),PROTOCOL,prior.INPUT,prior.VERIFICATION]),
        'whitelist':'Target short equations, source family equations and own-family rational parameter only. No generic or exceptional point, rank label, Selmer class or character.'})


def filename(token,family):return WORK/(token+'--'+family+'.json')


def pair(c,f):
    from sage.all import QQ,ZZ,GF,PolynomialRing,prime_range,pari
    R=PolynomialRing(QQ,'s');s=R.gen();A=R(f['A']);B=R(f['B']);a=QQ(c['a']);b=QQ(c['b'])
    assert a*b and 4*a**3+27*b**2 and A.degree()<=8 and B.degree()<=12
    F=b*b*A**3-a**3*B**2;assert F and F.degree()<=24
    enc=lambda q:list(map(str,q.list()))
    remainder=F;roots=[];control=None
    if f['family']==c['control_family']:
        t=QQ(c['control_parameter']);assert F(t)==0
        e=0
        while remainder(t)==0:remainder,rem=remainder.quo_rem(s-t);assert rem==0;e+=1
        roots.append((t,e));control={'parameter':str(t),'multiplicity':e}
    unit=QQ(remainder.leading_coefficient());remainder/=unit
    quotient=remainder;excluded=None;factors=[]
    if remainder.degree()>0:
        for p in prime_range(3,252):
            if any(x.denominator()%p==0 for x in remainder):continue
            Rp=PolynomialRing(GF(p),'s');q=Rp(remainder)
            assert q.degree()==remainder.degree()
            if q.gcd(pow(Rp.gen(),p,q)-Rp.gen()).degree()==0:
                excluded={'prime':int(p),'polynomial_mod_p':list(map(int,q.list()))};break
        if excluded is None:
            fac=pari(remainder).factor();prod=R(1)
            for i in range(len(fac[0])):
                q=R(fac[0][i]);e=int(fac[1][i]);q/=q.leading_coefficient();prod*=q**e
                factors.append({'coefficients_ascending':enc(q),'multiplicity':e,'degree':int(q.degree())})
                if q.degree()==1:roots.append((-q[0],e))
            assert prod==remainder
    preimages=[]
    if F.degree()<24:roots.append(('infinity',24-int(F.degree())))
    for t,e in roots:
        av,bv=(A[8],B[12]) if t=='infinity' else (A(t),B(t))
        rec={'parameter':str(t),'multiplicity':e,'source_a':str(av),'source_b':str(bv)}
        if 4*av**3+27*bv**2==0:
            rec.update({'status':'SINGULAR_SOURCE','rationally_isomorphic':False});preimages.append(rec);continue
        assert av*bv and b*b*av**3==a**3*bv**2
        u2=(b/bv)/(a/av);assert u2**2==a/av and u2**3==b/bv
        square=u2.is_square();rec.update({'status':'PASS','scaling_square':str(u2),
            'rationally_isomorphic':bool(square),'rational_scaling':str(u2.sqrt()) if square else None})
        preimages.append(rec)
    if control:assert any(x['parameter']==control['parameter'] and x['rationally_isomorphic'] for x in preimages)
    return {'token':c['token'],'family':f['family'],'status':'PASS','j_equation_ascending':enc(F),
        'finite_degree':int(F.degree()),'infinity_multiplicity':24-int(F.degree()),'control':control,
        'remaining_monic_polynomial':enc(quotient),'remaining_unit':str(unit),
        'root_free_modular_certificate':excluded,'rational_factorization':factors,'preimages':preimages,
        'bindings':bindings([Path(__file__),PROTOCOL,INPUT])}


def worker(token):
    from sage.all import pari
    pari.allocatemem(64000000,268435456,silent=True)
    data=r.read(INPUT);c=next(x for x in data['cases'] if x['token']==token)
    # Own-family first: control is always attempted before any cross-family case.
    families=sorted(data['families'],key=lambda f:(f['family']!=c['control_family'],f['family']))
    for f in families:
        result=pair(c,f);r.write_new(filename(token,f['family']),result)
        print(token,f['family'],result['status'],[(x['parameter'],x['rationally_isomorphic']) for x in result['preimages']],flush=True)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);data=r.read(INPUT);rows=[]
    for c in data['cases']:
        token=c['token'];log=WORK/(token+'.log')
        if not log.exists():
            with log.open('x') as h:
                try:
                    p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=h,stderr=h,timeout=15)
                    reason='Worker failure' if p.returncode else None
                except subprocess.TimeoutExpired:reason='15-second target worker cap'
            r.write_new(WORK/(token+'-terminal.json'),{'reason':reason})
        terminal=r.read(WORK/(token+'-terminal.json'))
        for f in data['families']:
            path=filename(token,f['family'])
            if path.exists():rows.append(r.read(path))
            else:
                assert terminal['reason'];rows.append({'token':token,'family':f['family'],'status':'UNKNOWN','reason':terminal['reason']})
        print(token,'completed',sum(x['token']==token and x['status']=='PASS' for x in rows),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.second-generic-presentation.v1','rows':rows,
        'bindings':bindings([Path(__file__),PROTOCOL,INPUT]),
        'boundary':'Equation-only rational j-preimages and rational-isomorphism tests. No generic-section union or new rank is inferred. A timeout is not absence of a presentation.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='worker':worker(args.token)
    else:globals()[args.mode]()

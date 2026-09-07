#!/usr/bin/env python3
"""Exact local scaling audit of the two saved higher-parameter retention pools."""
import argparse
from math import gcd
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'higher-displayed-reduction-scalings-v1';OUT=ART/'higher_displayed_reduction_scalings_v1.json';INPUTS=[LOCAL/n/'result.json' for n in ('higher32768-r17-extended-v1','higher32768-product-first-extended-v1')]
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',*INPUTS)}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve local scaling audit')
    for path in INPUTS:
        s=cert.read(path.parent/'replay.supervisor.json')
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('both retained score populations must replay first')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher-displayed-reduction-scalings.v1','sources':sources(),'prime_bound':65521,'maximum_distinct_models':10482,'gate':'Both score policies omit singular displayed reductions. Rationally isomorphic models can have different displayed good primes; audit whether p-power coefficient scaling hides good reduction in the actual retained population before attributing weak selection to rank incidence. The score policies and point cohorts are already frozen and remain unchanged.','scope':'For each saved integer short model and each prime5..65521 with p^4 dividing A and p^6 dividing B, divide repeatedly by p^4,p^6. Check the scaled discriminant modulo p. A nonsingular scaled model is an exact good-reduction witness for the same rational curve; a singular scaled model is reported only as still singular, without a reduction-type theorem. No new parameter, trace, point, factorization campaign or candidate reordering.','limits':{'wall_seconds':120,'rss_bytes':536870912,'workers':1},'boundaries':'Only the10482 retained address-models and score primes through65521. This does not audit the entire122368792-address population, primes2/3, minimal models globally, reduction types, rank scores as predictors, or actual ranks.'})
def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen local scaling inputs changed')
    models={}
    for path in INPUTS:
        for row in cert.read(path)['rows']:
            key=row['family'],row['parameter'];model=row['model']
            if key in models and models[key]!=model:raise ArithmeticError('same address has differing short model')
            models[key]=model
    if len(models)!=10482:raise ArithmeticError('fixed retained union differs')
    primes=_primes_up_to(65521);cache={};rows=[]
    for (family,parameter),model in sorted(models.items()):
        if model[:3]!=['0','0','0']:raise ArithmeticError('integer short models required')
        A,B=map(int,model[3:]);g=gcd(A,B)
        if g not in cache:
            rest=g;candidates=[]
            for q in primes:
                if q**4>rest:break
                v=0
                while rest%q==0:rest//=q;v+=1
                if q>=5 and v>=4:candidates.append(q)
            cache[g]=candidates
        for q in cache[g]:
            a,b=A,B;k=0
            while a%(q**4)==0 and b%(q**6)==0:a//=q**4;b//=q**6;k+=1
            if not k:continue
            if A!=a*q**(4*k) or B!=b*q**(6*k):raise ArithmeticError('exact rational curve scaling differs')
            delta=(-16*(4*pow(a,3,q)+27*pow(b,2,q)))%q
            rows.append({'family':family,'parameter':parameter,'prime':q,'scale_exponent':k,'scaled_A':str(a),'scaled_B':str(b),'scaled_discriminant_mod_p':delta,'status':'RECOVERED_GOOD_DISPLAY' if delta else 'SCALED_DISPLAY_STILL_SINGULAR'})
    return {'schema':'elliptic-curves.higher-displayed-reduction-scalings-result.v1','status':'PASS_EXACT_SCALING_AUDIT','sources':sources(),'protocol_hash':digest(p),'address_models':len(models),'scaling_cases':len(rows),'status_counts':dict(Counter(r['status'] for r in rows)),'prime_counts':dict(Counter(str(r['prime']) for r in rows)),'rows':rows,'claim_boundary':p['boundaries']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='replay':
            if cert.read(OUT)!=d:raise ArithmeticError('local scaling audit replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve local scaling audit result')
            checkpoint(OUT,d)
        print('LOCAL SCALING AUDIT',d['address_models'],d['scaling_cases'],d['status_counts'],d['prime_counts'],flush=True)

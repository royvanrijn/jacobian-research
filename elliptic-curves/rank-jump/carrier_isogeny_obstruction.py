#!/usr/bin/env python3
"""A scalar isogeny obstruction to one fixed everywhere-local genus-one carrier."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'CARRIER_ISOGENY_OBSTRUCTION_PROTOCOL.json'
SOURCE=r.OUT/'rank_jump_carrier_sha_class_v1.json'
RANK=r.OUT/'rank_jump_disjoint_soluble_carriers_v1.json'
OUTPUT=r.OUT/'rank_jump_carrier_isogeny_obstruction_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-carrier-isogeny-obstruction-v1'


def fingerprints(values):
    values=[r.F(str(x)) for x in values];primes=[];sigs=[0]*len(values);col=0
    for p in r.primes(503):
        if any(x.numerator%p==0 or x.denominator%p==0 for x in values):continue
        primes.append(p)
        for i,x in enumerate(values):sigs[i]|=int(pow(r.mod(x,p),(p-1)//2,p)==p-1)<<col
        col+=1
    return {'values':list(map(str,values)),'primes':primes,'fingerprints':sigs,'rank':r.rank(sigs)}


def worker():
    from sage.all import QQ,EllipticCurve,pari,PolynomialRing,FractionField
    from sage.version import version
    old=r.read(SOURCE);mapping=old['mapping'];E=EllipticCurve(QQ,list(map(QQ,mapping['Jacobian_model'])))
    torsion=[P for P in E.torsion_points() if P];assert len(torsion)==1
    e=torsion[0][0];assert torsion[0][1]==0
    a,b=3*e,3*e*e+E.a4();bp=a*a-4*b
    J=EllipticCurve(QQ,[0,a,0,b,0]);Jp=EllipticCurve(QQ,[0,-2*a,0,bp,0]);assert J.is_isomorphic(E)
    ring=PolynomialRing(QQ,'x');x=ring.gen();field=FractionField(ring)
    xp=x+a+b/x;mult=1-b/x**2
    assert (x**3+a*x*x+b*x)*mult**2==xp**3-2*a*xp*xp+bp*xp
    known=[E([QQ(x) for x in P]) for P in old['transported_rational_points_including_torsion']]
    alpha=[b if P==torsion[0] else P[0]-e for P in known]
    beta=sum(QQ(c)*e**i for i,c in enumerate(mapping['beta']))
    original=fingerprints(alpha+[beta]);base_rank=r.rank(original['fingerprints'][:-1])
    Em=Jp.minimal_model();pari.allocatemem(67108864,silent=True)
    ans=pari.ellrank(pari.ellinit(Em.a_invariants()),0)
    iso=Em.isomorphism_to(Jp)
    ps=[iso(Em([QQ(x) for x in P])) for P in ans[3]]+[Jp(0,0)]
    vals=[bp if P[0]==0 else P[0] for P in ps]
    other=fingerprints(vals)
    upper=r.read(RANK)['rows'][0]['descent']['rank_upper_bound'];assert upper==2
    full=base_rank==2 and other['rank']>=2
    obstructed=full and original['rank']>base_rank
    result={'schema':'rank-jump.carrier-isogeny-obstruction.v1',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,SOURCE,RANK,Path(__file__),HERE/'retrospective.py')},
            'software':{'sage':version,'pari':str(pari('version()'))},
            'Jacobian_root':str(e),'a':str(a),'b':str(b),'b_prime':str(bp),
            'original_scalar_classes_with_beta_last':original,'original_known_scalar_rank':base_rank,
            'isogenous_model':list(map(str,Jp.a_invariants())),
            'isogenous_minimal_model':list(map(str,Em.a_invariants())),
            'isogenous_points_including_2torsion':[list(map(str,P[:2])) for P in ps],
            'isogenous_scalar_classes':other,'raw_isogenous_ellrank':str(ans),
            'original_rank_upper_bound':upper,'original_scalar_image_complete':full,
            'original_scalar_image_dimension':2 if full else 'UNKNOWN',
            'carrier_nonzero_Sha_class_proved':obstructed,
            'carrier_global_solubility':'NO' if obstructed else 'UNKNOWN',
            'reason':'Two scalar classes on each isogeny side exhaust rank+2=4. Beta outside the original scalar image cannot come from a rational point on the carrier.'}
    r.write_new(WORK/'result.json',result)
    print('scalar dimensions',base_rank,other['rank'],'with beta',original['rank'],'Sha obstruction',obstructed,flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log'
    if not log.exists():
        with log.open('x') as f:
            try:
                p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],cwd=r.ROOT,stdout=f,stderr=f,timeout=30)
                status='COMPLETE' if p.returncode==0 else 'FAILED'
            except subprocess.TimeoutExpired:status='TIMEOUT'
        r.write_new(WORK/'execution.json',{'status':status})
    result=r.read(WORK/'result.json') if (WORK/'result.json').exists() else {'status':'UNKNOWN','execution':r.read(WORK/'execution.json'),'log':log.read_text()}
    r.write_new(OUTPUT,result);print(log.read_text());print(r.read(WORK/'execution.json'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);a=p.parse_args()
    worker() if a.mode=='worker' else run()

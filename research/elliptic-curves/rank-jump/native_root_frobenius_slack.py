#!/usr/bin/env python3
"""A bounded exact Frobenius test of shared native-twist capacity slack."""
import argparse
from pathlib import Path
from math import isqrt
import subprocess
import sys
import retrospective as r
import branch_divisibility_capacity as branch
import verify_branch_divisibility_capacity as verified
import root_curve_frobenius_capacity as count

PROTOCOL=Path(__file__).with_name('NATIVE_ROOT_FROBENIUS_SLACK_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_native_root_frobenius_slack_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-root-frobenius-slack-v1'


def compute():
    from sage.all import QQ,PolynomialRing
    d=r.read(branch.INPUT);v=r.read(verified.OUTPUT)
    assert v['status']=='PASS' and v['generic_mod_two_dimension']==17 and v['global_pool_upper_bound']==19
    for obj in (d,v):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    R=PolynomialRing(QQ,'t');A=R(d['A']);B=R(d['B']);D=-4*A**3-27*B**2
    f={'A':d['A'],'B':d['B'],'discriminant':list(map(str,D.list())),
       'squarefree_factors':[{'coefficients_ascending':list(map(str,D.list())),'multiplicity':1}]}
    rows=[]
    for p in r.primes(1009):
        ok=count.eligible(f,p)
        if ok is None:continue
        a,b,delta,nodes=ok;assert not nodes
        cubes=[x*x*x%p for x in range(p)];fibres=[]
        for t in range(p):
            aa,bb=count.value(a,t,p),count.value(b,t,p)
            fibres.append(sum((cubes[x]+aa*x+bb)%p==0 for x in range(p)))
        inf=[x for x in range(p) if (cubes[x]+a[8]*x+b[12])%p==0]
        total=sum(fibres)+len(inf);assert abs(total-p-1)<=isqrt(400*p)
        rows.append({'p':p,'finite_fibre_counts':fibres,'infinity_roots':inf,'point_count':total,
            'frobenius_trace':p+1-total,'odd_point_count':bool(total%2)})
        if len(rows)==12:break
    odd=[x['p'] for x in rows if x['odd_point_count']];cap=18 if odd else 19
    return {'schema':'rank-jump.native-root-frobenius-slack.v1','status':'PASS' if len(rows)==12 else 'UNKNOWN',
        'rows':rows,'odd_count_primes':odd,'global_pool_upper_bound':cap,
        'single_support_twist_upper_bound':cap-16,'multiple_support_twist_upper_bound':cap-17,
        'capacity_improved':bool(odd),'bindings':branch.bindings([Path(__file__),PROTOCOL,branch.INPUT,verified.OUTPUT,Path(count.__file__),Path(r.__file__)]),
        'boundary':'Even counts leave the existing bound unchanged. They do not prove rationality of the full two-torsion or determine the Galois image. No fibre-rank or point-search statistic.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        reason=None
        with (WORK/'worker.log').open('x') as log:
            try:
                proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=30)
                if proc.returncode:reason='worker failure'
            except subprocess.TimeoutExpired:reason='bounded timeout'
        if reason:r.write_new(OUTPUT,{'status':'UNKNOWN','reason':reason});return
    result=r.read(path);r.write_new(OUTPUT,result)
    print([(x['p'],x['point_count']) for x in result['rows']], 'capacity',result['global_pool_upper_bound'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'worker.json',compute())
    elif a.mode=='check':assert compute()==r.read(OUTPUT);print('PASS native root Frobenius replay')
    else:capture()

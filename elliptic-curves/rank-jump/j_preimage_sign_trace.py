#!/usr/bin/env python3
"""Modular certificates for irreducibility and the isomorphism sign obstruction."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import second_generic_presentation as source
import verify_second_generic_presentation as verified

PROTOCOL=Path(__file__).with_name('J_PREIMAGE_SIGN_TRACE_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_j_preimage_sign_trace_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_j_preimage_sign_trace_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-j-preimage-sign-trace-v1'


def export():
    from sage.all import QQ,PolynomialRing
    R=PolynomialRing(QQ,'s');s=R.gen();data=r.read(source.INPUT);out=r.read(source.OUTPUT)
    assert r.read(verified.OUTPUT)['status']=='PASS';rows=[]
    for x in out['rows']:
        q=R(x['j_equation_ascending'])
        for p in x['preimages']:
            assert p['parameter']!='infinity';q,rem=q.quo_rem((s-QQ(p['parameter']))**p['multiplicity']);assert rem==0
        q/=q.leading_coefficient();assert q.degree() in (23,24)
        rows.append({'token':x['token'],'family':x['family'],'polynomial_ascending':list(map(str,q.list()))})
    r.write_new(INPUT,{'schema':'rank-jump.j-preimage-sign-trace-inputs.v1','families':data['families'],
        'cases':[{k:c[k] for k in ('token','a','b')} for c in data['cases']],'rows':rows,
        'bindings':source.bindings([Path(__file__),PROTOCOL,source.INPUT,source.OUTPUT,verified.OUTPUT])})


def filename(row):return WORK/(row['token']+'--'+row['family']+'.json')


def cell(row,c,f):
    from sage.all import QQ,GF,PolynomialRing,prime_range
    q=list(map(QQ,row['polynomial_ascending']));a,b=QQ(c['a']),QQ(c['b']);A=list(map(QQ,f['A']));B=list(map(QQ,f['B']))
    n=len(q)-1;possible=set(range(1,n));certs=[];sign=None
    for p in prime_range(3,252):
        if any(v.denominator()%p==0 for v in q+A+B+[a,b]):continue
        R=PolynomialRing(GF(p),'s');x=R.gen();qp=R(q);ap,bp=GF(p)(a),GF(p)(b)
        if not ap*bp or qp.degree()!=n or qp.gcd(qp.derivative())!=1:continue
        fac=qp.factor();parts=[(h,int(e)) for h,e in fac];assert all(e==1 for h,e in parts)
        degrees=[int(h.degree()) for h,e in parts];sums={0}
        for d in degrees:sums|={s+d for s in list(sums)}
        possible&=sums-{0,n}
        cert={'prime':int(p),'factors_ascending':[list(map(int,h.list())) for h,e in parts],
            'factor_degrees':degrees,'remaining_possible_factor_degrees':sorted(possible)}
        certs.append(cert)
        if sign is None:
            num=bp*R(A);den=ap*R(B)
            for i,(h,e) in enumerate(parts):
                if not num%h or not den%h:continue
                residue=(num*den.inverse_mod(h))%h
                leg=pow(residue,(int(p)**int(h.degree())-1)//2,h)
                assert leg in (R(1),R(-1))
                if leg==R(-1):
                    sign={'prime':int(p),'factor_index':i,'residue_factor_ascending':list(map(int,h.list())),
                        'scaling_class_residue_ascending':list(map(int,residue.list())),
                        'finite_field_degree':int(h.degree()),'quadratic_character':-1};break
        if not possible and sign is not None:break
    return {'token':row['token'],'family':row['family'],'degree':n,'status':'PASS' if not possible and sign is not None else 'UNKNOWN',
        'irreducibility_certificates':certs,'remaining_possible_factor_degrees':sorted(possible),
        'nonsquare_isomorphism_class_witness':sign,
        'rational_trace_span_dimension':0 if not possible and sign is not None else None,
        'bindings':source.bindings([Path(__file__),PROTOCOL,INPUT])}


def worker():
    from sage.all import pari
    pari.allocatemem(64000000,268435456,silent=True)
    data=r.read(INPUT);cases={x['token']:x for x in data['cases']};families={x['family']:x for x in data['families']}
    for row in data['rows']:
        result=cell(row,cases[row['token']],families[row['family']]);r.write_new(filename(row),result)
        print(row['token'],row['family'],result['status'],len(result['irreducibility_certificates']),flush=True)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);data=r.read(INPUT)
    with (WORK/'worker.log').open('x') as h:
        try:
            p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=h,stderr=h,timeout=30)
            reason='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:reason='30-second total worker cap'
    rows=[]
    for row in data['rows']:
        path=filename(row)
        if path.exists():rows.append(r.read(path))
        else:
            assert reason;rows.append({'token':row['token'],'family':row['family'],'status':'UNKNOWN','reason':reason})
    r.write_new(OUTPUT,{'schema':'rank-jump.j-preimage-sign-trace.v1','rows':rows,
        'bindings':source.bindings([Path(__file__),PROTOCOL,INPUT]),
        'boundary':'Point-free modular evidence for a trace-cancellation obstruction on retained preimage fields. Does not exclude other rational points or other block constructors.'})
    print('completed',sum(x['status']=='PASS' for x in rows),'of',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);args=p.parse_args();globals()[args.mode]()

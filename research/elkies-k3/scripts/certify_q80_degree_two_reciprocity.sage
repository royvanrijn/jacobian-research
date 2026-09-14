#!/usr/bin/env sage
"""Freeze/replay the complete degree-two Q80 reciprocity reduction gate."""
from sage.all import GF, PolynomialRing, QQ
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import argparse
import json
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
CHECKER='elkies-k3/scripts/verify_q80_degree_two_reciprocity.py'
OLD_CHECKER='elkies-k3/scripts/verify_q80_single_branch_reciprocity.py'
OLD_INPUT='artifacts/generated-results/elkies-k3-q80-single-branch-reciprocity-v1/input.json'
def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def write(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def freeze():
    files=[SOURCE,CHECKER,OLD_CHECKER,OLD_INPUT,str(Path(__file__).relative_to(ROOT))]
    packet={'schema':1,'prime':131,'field_modulus_low_to_high':[1,0,1],
            'section_indices':list(range(17)),'cpu_seconds':60,'memory_bytes':4*1024**3,
            'source':SOURCE,'bindings':{s:digest(ROOT/s) for s in files},
            'scope':'Complete quadratic residue norm characters; collision boundary retained; no coefficient lift or positive MW17 endpoint.'}
    write('input.json',packet)
    print(json.dumps({'status':'FROZEN','input_sha256':digest(OUT/'input.json')}),flush=True)
def run():
    packet=read(OUT/'input.json')
    for s,h in packet['bindings'].items():assert digest(ROOT/s)==h
    assert packet['prime']==131 and packet['section_indices']==list(range(17))
    resource.setrlimit(resource.RLIMIT_CPU,(packet['cpu_seconds'],packet['cpu_seconds']+5))
    resource.setrlimit(resource.RLIMIT_AS,(packet['memory_bytes'],packet['memory_bytes']))
    start=time.process_time();source=read(ROOT/SOURCE);p=131
    F=GF(p);R=PolynomialRing(F,'t');K=R.fraction_field()
    F2=GF(p*p,'a',modulus=PolynomialRing(F,'a').gen()**2+1)
    Q=PolynomialRing(F2,'x');x=Q.gen()
    def poly(values):return R([F(QQ(v)) for v in values])
    A=poly(source['weierstrass_model']['A_coefficients_low_to_high'])
    B=poly(source['weierstrass_model']['B_coefficients_low_to_high'])
    sections=[]
    for row in source['sections']['records']:
        coords=[]
        for key in ['X','Y']:
            v=row[key]
            coords.append(K(poly(v['numerator_coefficients_low_to_high']))/poly(v['denominator_coefficients_low_to_high']))
        assert coords[1]**2==coords[0]**3+A*coords[0]+B
        sections.append(coords)
    def encode(z):
        c=z.polynomial().list()
        return int(c[0])+(p*int(c[1]) if len(c)>1 else 0) if c else 0
    def at(f,b,weight):
        n,d=f.numerator(),f.denominator()
        if b is not None:return n(b)/d(b) if d(b) else None
        gap=d.degree()+weight-n.degree()
        if gap<0:return None
        return F2(0) if gap else F2(n.leading_coefficient()/d.leading_coefficient())
    def chars(b,roots,a,bb,degree):
        points=[]
        for X,Y in sections:
            u,v=at(X,b,4),at(Y,b,6)
            assert (u is None)==(v is None)
            assert u is None or v*v==u**3+a*u+bb
            points.append((u,v))
        codes=[]
        for e in roots:
            code=0
            for i,(u,v) in enumerate(points):
                value=F2(1) if u is None else (3*e*e+a if u==e else u-e)
                assert value, 'An actual inherited section cannot meet an I1 node.'
                if degree==1:value=F(value)
                if not value.is_square():code|=1<<i
            codes.append(code)
        return codes
    def row(b,degree):
        a,bb=(F2(A[8]),F2(B[12])) if b is None else (F2(A(b)),F2(B(b)))
        roots=sorted([v for v,m in (x**3+a*x+bb).roots() for _ in range(int(m))],key=encode)
        r={'t':None if b is None else encode(b),'smooth':bool(4*a**3+27*bb**2),
           'roots':[encode(v) for v in roots], 'norm_codes':chars(b,roots,a,bb,2)}
        if degree==1:
            rational_roots=[v for v in roots if v**p==v]
            r['rational_roots']=[encode(v) for v in rational_roots]
            r['rational_codes']=chars(b,rational_roots,a,bb,1)
        return r
    rational=[row(F2(i),1) for i in range(p)]+[row(None,1)]
    write('rational-fibres.json',{'rows':rational})
    quadratic=[row(F2(i)+j*F2.gen(),2) for i in range(p) for j in range(1,66)]
    write('quadratic-fibres.json',{'rows':quadratic})
    groups=defaultdict(list); split=[]
    for r in rational:
        assert r['smooth']
        if len(r['rational_roots'])==3:
            split.append(r['t'])
            for i in range(3):
                for j in range(3):
                    if i!=j:groups[(r['rational_codes'][i],r['rational_codes'][j])].append([r['t'],r['roots'][i],r['roots'][j]])
    assert len(groups)==96 and all(len(v)==1 for v in groups.values())
    assert all(r['norm_codes'].count(0)<2 for r in quadratic if len(r['roots'])==3)
    assert all(r['norm_codes'].count(0)<2 for r in rational if len(r['roots'])==3 and len(r['rational_roots'])==1)
    hist=Counter((len(r['roots']),r['smooth']) for r in quadratic)
    result={'status':'PASS','input_sha256':digest(OUT/'input.json'),
        'records':{s:digest(OUT/s) for s in ['rational-fibres.json','quadratic-fibres.json']},
        'rational_base_values':132,'quadratic_base_orbits':8515,'rational_full_splitting_values':split,
        'rational_ordered_pair_codes':96,'rational_ordered_pair_collisions':0,
        'quadratic_histogram':[[n,s,v] for (n,s),v in sorted(hist.items())],
        'bad_quadratic_fibres':[r for r in quadratic if not r['smooth']],
        'quadratic_two_zero_code_triples':0,'rational_nonsplit_two_zero_code_triples':0,
        'genus0_branch_reduction':'double finite point in the sixteen full-splitting residues',
        'collision_strata_excluded':False,'coefficient_lift_complete':False,'positive_mw17_target_complete':False}
    write('result.json',result)
    write('execution.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),
        'result_sha256':digest(OUT/'result.json'),'cpu_seconds':time.process_time()-start,
        'cpu_limit_seconds':60,'memory_limit_bytes':4*1024**3})
    print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['freeze','run'])
    args=parser.parse_args();freeze() if args.action=='freeze' else run()

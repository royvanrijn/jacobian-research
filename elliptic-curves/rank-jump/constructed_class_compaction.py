#!/usr/bin/env python3
"""Expand certified strict classes after exact half-ideal square reduction."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import constructed_class_half_ideal as half

PROTOCOL=Path(__file__).with_name('CONSTRUCTED_CLASS_COMPACTION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_constructed_class_compaction_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-constructed-class-compaction-v1'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    ref=r.read(half.REFERENCE);c=r.read(half.CONSTRUCTION);h=r.read(half.OUTPUT);policy=r.read(PROTOCOL)
    pari.allocatemem(64000000,policy['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    factors={('generic',i):pari.Mod(pari(R(a['beta_ascending'])),pari(f)) for i,a in enumerate(ref['generic_classes'])}
    factors.update({('projected_atom',a['index']):ZZ(a['norm'])*pari.Mod(pari(R(a['alpha_ascending'])),pari(f)) for a in c['atoms']})
    cases=[]
    for record in h['columns'][6:]:
        bases={};exponents={}
        def include(value,e):
            key=str(pari.nfalgtobasis(nf,value));bases[key]=pari(key);exponents[key]=exponents.get(key,0)+e
        for label in record['factor_labels']:include(factors[(label['kind'],label['index'])],1)
        for step in record['half_ideal_reduction_steps']:include(pari.nfbasistoalg(nf,pari(step['principal_multiplier_GP'])),-2)
        keys=[k for k,e in exponents.items() if e]
        path=WORK/('column_%02d_factorization.json'%record['column'])
        r.write_new(path,{'column':record['column'],'factor_basis_GP':keys,'exponents':[exponents[k] for k in keys]})
        print('EXPANDING',record['column'],'FACTORS',len(keys),flush=True)
        vector=pari.nffactorback(nf,[bases[k] for k in keys],[exponents[k] for k in keys])
        beta=pari.nfbasistoalg(nf,vector);J=pari(record['final_reduced_half_ideal_hnf'])
        assert pari.idealhnf(nf,beta)==pari.idealpow(nf,J,2)
        N=QQ(pari.nfeltnorm(nf,beta));assert N==QQ(pari.idealnorm(nf,J))**2
        coeff=[QQ(pari.lift(beta).polcoef(i)) for i in range(3)]
        size=max(max(abs(a.numerator()).nbits(),a.denominator().nbits()) for a in coeff)
        result={'column':record['column'],'beta_ascending':list(map(str,coeff)),
            'max_coefficient_bits':int(size),'norm':str(N),'half_ideal_hnf':str(J),
            'square_equivalence_factorization':str(path.relative_to(r.ROOT)),
            'square_equivalence_sha256':r.digest(path.read_bytes())}
        if size<=policy['bounds']['maximum_cover_coefficient_bits']:
            V=PolynomialRing(QQ,['u','v','w','s']);u,v,w,s=V.gens();T=PolynomialRing(V,'theta');theta=T.gen()
            ff=T(f.list());b=T(coeff);q=(b*(u+v*theta+w*theta**2)**2)%ff
            result['cover_quadrics']=[str(q[2]),str(q[1]+s*s)]
            result['x_numerator']=str(q[0]);result['status']='PASS'
        else:result['status']='REPRESENTATIVE_ONLY_COVER_OUTPUT_BOUND'
        r.write_new(WORK/('column_%02d.json'%record['column']),result);cases.append(result)
        print('COMPACT',record['column'],'BITS',size,'STATUS',result['status'],flush=True)
    return {'schema':'rank-jump.constructed-class-compaction.v1','status':'PASS','cases':cases,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,half.REFERENCE,half.CONSTRUCTION,half.OUTPUT]},
        'boundary':'Explicit square-equivalent representatives and, where bounded, two-cover equations. Rational solubility remains UNKNOWN.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()

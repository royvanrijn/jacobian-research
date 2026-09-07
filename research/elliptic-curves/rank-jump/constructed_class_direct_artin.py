#!/usr/bin/env python3
"""Exact local Artin evaluation of fixed half ideals, without coprimality heuristics."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import constructed_class_half_ideal as half

PROTOCOL=Path(__file__).with_name('CONSTRUCTED_CLASS_DIRECT_ARTIN_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_constructed_class_direct_artin_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-constructed-class-direct-artin-v1'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    h=r.read(half.OUTPUT);c=r.read(half.CONSTRUCTION);ref=r.read(half.REFERENCE)
    for path,sha in h['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    pari.allocatemem(64000000,r.read(PROTOCOL)['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    factors={('generic',i):pari.Mod(pari(R(a['beta_ascending'])),pari(f)) for i,a in enumerate(ref['generic_classes'])}
    factors.update({('projected_atom',a['index']):ZZ(a['norm'])*pari.Mod(pari(R(a['alpha_ascending'])),pari(f)) for a in c['atoms']})
    cache={};columns=[]
    for column in h['columns']:
        I=pari(column['final_reduced_half_ideal_hnf']);fac=pari.idealfactor(nf,I)
        assert pari.idealhnf(nf,pari.idealfactorback(nf,fac))==I
        parts=[];bits=[0]*len(h['character_labels'])
        for j in range(fac.nrows()):
            P,e=fac[j,0],int(fac[j,1]);key=str(pari.idealhnf(nf,P))
            if key not in cache:
                pi=pari.nfbasistoalg(nf,pari.idealappr(nf,P));assert int(pari.idealval(nf,pi,P))==1
                symbols={key:int(pari.nfhilbert(nf,a,pi,P)==-1) for key,a in factors.items()}
                row=[sum(symbols[(a['kind'],a['index'])] for a in char)%2 for char in h['character_labels']]
                cache[key]={'p':int(P[0]),'hnf':key,'e':int(P[2]),'f':int(P[3]),'uniformizer_GP':str(pi),'artin_bits':row}
            parts.append({'prime_hnf':key,'valuation':e})
            if e%2:bits=[x^y for x,y in zip(bits,cache[key]['artin_bits'])]
        if column['status']=='PASS':assert bits==column['artin_bits']
        record={'column':column['column'],'role':column['role'],'reduced_half_ideal_hnf':str(I),
            'ideal_norm':str(pari.idealnorm(nf,I)),'prime_factors':parts,'artin_bits':bits}
        columns.append(record);r.write_new(WORK/('column_%02d.json'%column['column']),record)
        print('COLUMN',column['column'],'ARTIN',bits,flush=True)
    packed=[r.pack(x['artin_bits']) for x in columns];grank=r.rank(packed[:6]);rank=r.rank(packed)
    return {'schema':'rank-jump.constructed-class-direct-artin.v1','status':'PASS','columns':columns,
        'prime_artin_certificates':list(cache.values()),'generic_half_ideal_artin_rank':grank,
        'combined_half_ideal_artin_rank':rank,'detected_relative_half_ideal_rank':rank-grank,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,half.OUTPUT,half.CONSTRUCTION,half.REFERENCE]},
        'boundary':'An Artin rank increase proves extra ordinary half-ideal information beyond the generic strict images. Vanishing does not prove principality, unit representation, CT vanishing or rational solubility.'}

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

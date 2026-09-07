#!/usr/bin/env python3
"""An exact affine obstruction against the entire retained norm dictionary."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import early_relation_pool as pool
import relation_root_class as root
import complete_relation_root_class as completed

PROTOCOL=Path(__file__).with_name('RELATION_ROOT_AFFINE_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_relation_root_affine_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-relation-root-affine-v1'


def dictionary(d):
    blocks={}
    for i,c in enumerate(d['columns']):blocks.setdefault(c['p'],[]).append(i)
    rows=[];S=set(d['S_finite'])
    for rel in d['relations']:
        vals=dict(rel['ideal_factorization']);v=0
        for p in {d['columns'][i]['p'] for i in vals}-S:
            n=sum(d['columns'][i]['f']*vals.get(i,0) for i in blocks[p])
            for i in blocks[p]:v|=((vals.get(i,0)+d['columns'][i]['e']*n)%2)<<i
        rows.append(v)
    return rows,blocks


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    pari.allocatemem(64000000,268435456,silent=True)
    d=r.read(pool.INPUT);old=r.read(pool.OUTPUT);inp=r.read(root.INPUT);out=r.read(completed.OUTPUT)
    for obj in (old,inp,out):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    rows,blocks=dictionary(d);piv=r.basis(rows)
    assert len(rows)==old['element_count']==4134 and len(piv)==old['outside_S_parity_rank']==4133
    R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);nf=pari.nfinit([pari(f),d['S_finite']])
    th=pari.Mod('z',pari(f));beta=pari(R(out['projection_ascending']))(th)
    gamma=[pari(R(g))(th) for g in inp['generic_classes_ascending']]
    S=set(d['S_finite']);valuations=[None]*len(d['columns']);ideals={};target=0
    for p,ids in sorted(blocks.items()):
        if p in S:continue
        dec=list(pari.idealprimedec(nf,p));assert len(dec)==len(ids)
        for i,P in zip(ids,dec):
            c=d['columns'][i];assert str(pari.idealhnf(nf,P))==c['hnf']
            assert int(P[2])==c['e']==1 and int(P[3])==c['f']
            v=int(pari.idealval(nf,beta,P));valuations[i]=v;ideals[i]=P;target|=(v%2)<<i
    residual=r.reduce(target,piv);dual=0;witness=[]
    if residual:
        j=residual.bit_length()-1;assert j not in piv;dual=1<<j
        for k in sorted(piv):
            if (dual&piv[k]).bit_count()%2:dual|=1<<k
        assert all((dual&row).bit_count()%2==0 for row in rows) and (dual&target).bit_count()%2==1
        for i,c in enumerate(d['columns']):
            if not dual>>i&1:continue
            assert c['p'] not in S and valuations[i] is not None
            gv=[int(pari.idealval(nf,g,ideals[i])) for g in gamma]
            assert all(v%2==0 for v in gv)
            witness.append({'column':i,**c,'projection_valuation':valuations[i],'generic_valuations':gv})
    return {'schema':'rank-jump.relation-root-affine.v1','status':'PASS',
        'affine_status':'INCONSISTENT' if residual else 'RESTRICTED_SYSTEM_CONSISTENT',
        'dictionary_rows':len(rows),'dictionary_parity_rank':len(piv),
        'augmented_parity_rank':len(piv)+int(bool(residual)),
        'retained_prime_blocks':len(blocks),'retained_prime_ideal_columns':len(d['columns']),
        'projection_valuations_by_column':valuations,'odd_projection_columns':target.bit_count(),
        'dual_witness':witness,'dual_dictionary_values':[int((dual&row).bit_count()%2) for row in rows],
        'dual_projection_value':int((dual&target).bit_count()%2),
        'new_unramified_classes_from_root_dictionary_coset':0 if residual else 'UNKNOWN',
        'boundary':'Inconsistency excludes every retained norm-dictionary and generic correction of the fixed relation root. This is not a full class-group bound, not a new fibre search, and not a solubility obstruction inside Selmer.',
        'bindings':root.bindings([Path(__file__),PROTOCOL,pool.INPUT,pool.OUTPUT,root.INPUT,completed.OUTPUT,
            Path(r.__file__)])}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        reason=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=60)
                if p.returncode:reason='Worker failure'
            except subprocess.TimeoutExpired:reason='Bounded worker timeout'
        if reason:r.write_new(OUTPUT,{'status':'UNKNOWN','reason':reason});return
    result=r.read(path);r.write_new(OUTPUT,result)
    print(result['affine_status'],'dual support',len(result['dual_witness']),'primes',sorted({x['p'] for x in result['dual_witness']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'worker.json',compute())
    elif a.mode=='check':assert compute()==r.read(OUTPUT);print('PASS affine replay')
    else:capture()

#!/usr/bin/env python3
"""Exact ramification capacity of a retained equation-only principal dictionary."""
import argparse
from pathlib import Path
import re
import subprocess
import sys
from math import prod
import retrospective as r
import bounded_gain_reference as ref

PROTOCOL=Path(__file__).with_name('RETAINED_NORM_BATCH_CAPACITY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_retained_norm_batch_capacity_inputs_v1.json'
PROVENANCE=r.OUT/'rank_jump_retained_norm_batch_capacity_provenance_v1.json'
OUTPUT=r.OUT/'rank_jump_retained_norm_batch_capacity_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-retained-norm-batch-capacity-v1'
REL=r.OUT/'small_conductor_norm_batch_relations_v1.json'
FORM=r.OUT/'small_conductor_reduced_norm_form_v1.json'
BATCH=r.OUT/'small_conductor_norm_batch_v1.json'


def export():
    rel=r.read(REL);form=r.read(FORM);batch=r.read(BATCH);base=r.read(ref.OUTPUT)['stages']
    assert rel['status']==form['status']==batch['status']=='PASS'
    assert batch['protocol']['box']==512 and batch['protocol']['smooth_bound']==400000
    cols=rel['columns'];ps={cols[i]['p'] for row in rel['relations'] for i,e in row['ideal_factorization']}
    ids=[i for i,c in enumerate(cols) if c['p'] in ps];index={old:new for new,old in enumerate(ids)}
    rows=[{'m':x['m'],'n':x['n'],'alpha_ascending':x['beta_power_basis'],
           'ideal_factorization':[[index[i],e] for i,e in x['ideal_factorization']]} for x in rel['relations']]
    scale=r.F(base['factor']['scale'])
    r.write_new(INPUT,{'schema':'rank-jump.retained-norm-batch-capacity-inputs.v1',
        'cubic_ascending':list(reversed(form['original_monic_cubic_descending'])),
        'field_discriminant':form['field_discriminant'],
        'masked_cubic_ascending':base['factor']['integral_cubic_ascending'],
        'affine_map_old_root_from_masked':['-1/3',str(4/scale**2)],
        'S_finite':base['local']['S_finite'],'norm_generator':form['integral_norm_generator'],
        'sl2_matrix':form['sl2_matrix'],'columns':[cols[i] for i in ids],'relations':rows})
    files=(REL,FORM,BATCH,ref.OUTPUT,INPUT,PROTOCOL,
        r.ROOT/'elliptic-curves/cas/extend_small_conductor_norm_batch.py',
        r.ROOT/'elliptic-curves/cas/prepare_small_conductor_norm_form.py',
        r.ROOT/'elliptic-curves/cas/audit_small_conductor_norm_batch.sage')
    r.write_new(PROVENANCE,{'schema':'rank-jump.retained-norm-batch-capacity-provenance.v1',
        'selection':'Retrospective fixed-box relation batch; later targeted/protected/point-derived waves excluded.',
        'whitelist':'Equation and norm-map data, all principal relation rows, complete prime ideal blocks at their norm primes. No class bound, outcome label, exceptional point, anchor or character.',
        'retained_relation_count':len(rows),'source_fixed_box':batch['protocol']['box'],
        'original_primitive_pair_count':batch['primitive_pairs'],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT,Path(r.__file__))}


def setup():
    from sage.all import QQ,ZZ,pari,PolynomialRing
    data=r.read(INPUT);R=PolynomialRing(QQ,'x');f=R(data['cubic_ascending']);F=R(data['masked_cubic_ascending'])
    b,a=map(QQ,data['affine_map_old_root_from_masked']);assert f(a*R.gen()+b)==a**3*F
    remaining=abs(ZZ(f.discriminant()))
    for p in data['S_finite']:
        assert ZZ(p).is_prime(proof=True)
        while remaining%p==0:remaining//=p
    assert remaining==1
    nf=pari.nfinit([pari(f),data['S_finite']]);assert str(nf.disc())==data['field_discriminant']
    th=pari.Mod('x',pari(f));cols=data['columns'];ideals=[];blocks={};decompositions={}
    for i,c in enumerate(cols):
        assert re.fullmatch(r'[\[\]0-9;, \-]+',c['hnf'])
        H=pari(c['hnf']);ideals.append(H);blocks.setdefault(c['p'],[]).append(i)
    for p,ids in blocks.items():
        assert ZZ(p).is_prime(proof=True);dec=list(pari.idealprimedec(nf,p));decompositions[p]=dec
        assert len(dec)==len(ids)
        for i,P in zip(ids,dec):
            assert pari.idealhnf(nf,P)==ideals[i] and int(P[2])==cols[i]['e'] and int(P[3])==cols[i]['f']
        assert sum(cols[i]['e']*cols[i]['f'] for i in ids)==3
    return data,f,nf,th,ideals,blocks,decompositions


def compute():
    from sage.all import QQ,ZZ,pari
    pari.allocatemem(64000000,268435456,silent=True)
    data,f,nf,th,ideals,blocks,dec=setup();cols=data['columns'];S=set(data['S_finite']);rows=[];masks=[]
    a=QQ(data['norm_generator']['fixed_a']);w=pari(f.parent()(data['norm_generator']['w_power_basis']))(th);M=data['sl2_matrix']
    for rel in data['relations']:
        m,n=rel['m'],rel['n'];assert abs(m)<=512 and 1<=n<=512 and ZZ(m).gcd(n)==1
        alpha=pari(f.parent()(rel['alpha_ascending']))(th)
        assert alpha==pari(a)*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
        I=pari.idealhnf(nf,1);vals=dict(rel['ideal_factorization']);norm=ZZ(pari.nfeltnorm(nf,alpha));assert norm
        assert len(vals)==len(rel['ideal_factorization']) and all(e>0 for e in vals.values())
        for i,e in vals.items():I=pari.idealmul(nf,I,pari.idealpow(nf,ideals[i],e))
        assert pari.idealhnf(nf,I)==pari.idealhnf(nf,alpha)
        assert prod(cols[i]['p']**(cols[i]['f']*e) for i,e in vals.items())==abs(norm)
        projection=pari(norm)*alpha;assert pari.nfeltnorm(nf,projection)==norm**4
        mask=0;sparse=[]
        for p in sorted({cols[i]['p'] for i in vals}):
            if p in S:continue
            v=sum(cols[i]['f']*vals.get(i,0) for i in blocks[p])
            for i in blocks[p]:
                e=vals.get(i,0)+cols[i]['e']*v
                if e%2:mask|=1<<i;sparse.append(i)
        masks.append(mask);rows.append({'m':m,'n':n,'norm_alpha':str(norm),'outside_parity_columns':sparse})
    rank=r.rank(masks)
    return {'schema':'rank-jump.retained-norm-batch-capacity.v1','status':'PASS','bindings':bindings(),
        'rows':rows,'relation_count':len(rows),'complete_prime_ideal_blocks':len(blocks),'column_count':len(cols),
        'outside_S_parity_rank':rank,'coefficient_kernel_dimension':len(rows)-rank,
        'additional_strict_class_dimension_in_projected_dictionary':0 if rank==len(rows) else 'UNKNOWN',
        'boundary':'Full parity rank excludes all nonempty norm-projection products from strict Selmer, including after generic correction. This is a dictionary obstruction, not a class-group upper bound or a Sha obstruction.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        error=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=60)
                if p.returncode:error='Worker failure'
            except subprocess.TimeoutExpired:error='Bounded worker timeout'
        if error:r.write_new(path,{'status':'UNKNOWN','reason':error,'bindings':bindings()})
    result=r.read(path);assert result['bindings']==bindings();r.write_new(OUTPUT,result)
    print({k:v for k,v in result.items() if k not in ('rows','bindings')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/'worker.json',compute())
    else:globals()[args.mode]()

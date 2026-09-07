#!/usr/bin/env python3
"""Bounded replay of the entire retained pre-protected norm dictionary."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from math import prod,gcd
import retrospective as r
import retained_norm_batch_capacity as prior
import retained_norm_inherited_hit as hit

PROTOCOL=Path(__file__).with_name('EARLY_RELATION_POOL_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_early_relation_pool_inputs_v1.json'
PROVENANCE=r.OUT/'rank_jump_early_relation_pool_provenance_v1.json'
OUTPUT=r.OUT/'rank_jump_early_relation_pool_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-early-relation-pool-v1'


def paths():
    names=['small_conductor_norm_batch_relations_v1.json','small_conductor_special_prime_relations_v1.json','small_conductor_small_base_relations_v1.json']
    names += [f'small_conductor_class_target_wave_{i:03d}_v1.json' for i in range(1,5)]
    names += [f'small_conductor_class_target_strip_wave_{i:03d}_v1.json' for i in range(1,6)]
    return [r.OUT/n for n in names]


def export():
    sources={};seen={};phases=[];orig=r.read(prior.REL);cols=orig['columns']
    def read(p):
        raw=p.read_bytes();sources[str(p.relative_to(r.ROOT))]=r.digest(raw);return json.loads(raw)
    for p in paths():
        d=read(p);assert d['status']=='PASS';occ=0
        for name,sha in d['sources'].items():
            q=r.ROOT/name;assert r.digest(q.read_bytes())==sha;sources[name]=sha
        sets=[d['relations']] if 'relations' in d else []
        for chunk in d.get('chunks',[]):
            q=r.ROOT/chunk['path'];assert r.digest(q.read_bytes())==chunk['sha256'];sets.append(read(q)['relations'])
        for rows in sets:
            for v in rows:
                m,n=v['m'],v['n'];assert gcd(m,n)==1 and n
                sign=-1 if n<0 else 1;key=(sign*m,sign*n)
                row={'m':key[0],'n':key[1],'alpha_ascending':[str(sign*r.F(c)) for c in v['beta_power_basis']],
                     'ideal_factorization':v['ideal_factorization']}
                if key in seen:assert row==seen[key]
                else:seen[key]=row
                occ+=1
        phases.append({'source':p.name,'retained_occurrences':occ,'cumulative_unique_elements':len(seen),
                       'reported_stop_dimension':d.get('matrix',{}).get('quotient_dimension')})
    rows=list(seen.values());support={cols[i]['p'] for v in rows for i,e in v['ideal_factorization']}
    ids=[i for i,c in enumerate(cols) if c['p'] in support];mapping={c:i for i,c in enumerate(ids)}
    for row in rows:row['ideal_factorization']=[[mapping[i],e] for i,e in row['ideal_factorization']]
    d=read(prior.INPUT);d.update(schema='rank-jump.early-relation-pool-inputs.v1',columns=[cols[i] for i in ids],relations=rows,phases=phases)
    # Compact immutable JSON avoids repeating indentation for thousands of rows.
    with INPUT.open('x') as f:json.dump(d,f,separators=(',',':'),sort_keys=True);f.write('\n')
    sources[str(INPUT.relative_to(r.ROOT))]=r.digest(INPUT.read_bytes())
    r.write_new(PROVENANCE,{'schema':'rank-jump.early-relation-pool-provenance.v1',
        'bindings':sources,'protocol_sha256':r.digest(PROTOCOL.read_bytes()),
        'selection_audit':{'early_matrix_sources':'principal norm relations only; code paths baseline/prior_state/selection inspected',
            'equation_sources':'maximal-order norm form, discriminant hints, Hessian and real roots',
            'outcome_dependence':'target dimension16 and chosen research effort are retrospective; no phase reached target16',
            'excluded':'all protected, capped and residual waves, including point-derived character anchors and principal rows',
            'not_replayed':'original exhaustive sieves, all target-selection decisions and adaptive misses'},
        'phases':phases})


def setup():
    old=prior.INPUT
    try:prior.INPUT=INPUT;return prior.setup()
    finally:prior.INPUT=old


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),PROTOCOL,INPUT,Path(prior.__file__),Path(r.__file__),hit.OUTPUT)}


def elimination(rows):
    piv={};kernel=[]
    for i,row in enumerate(rows):
        mask=row;comb=1<<i
        while mask:
            p=mask.bit_length()-1
            if p not in piv:piv[p]=(mask,comb);break
            v,c=piv[p];mask^=v;comb^=c
        if not mask:kernel.append(comb)
    return len(piv),kernel


def compute():
    from sage.all import QQ,ZZ,pari
    pari.allocatemem(64000000,268435456,silent=True)
    d,f,nf,th,ideals,blocks,dec=setup();cols=d['columns'];S=set(d['S_finite']);masks=[];norms=[]
    a=QQ(d['norm_generator']['fixed_a']);w=pari(f.parent()(d['norm_generator']['w_power_basis']))(th);M=d['sl2_matrix']
    for row in d['relations']:
        m,n=row['m'],row['n'];alpha=pari(f.parent()(row['alpha_ascending']))(th)
        assert alpha==pari(a)*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
        I=pari.idealhnf(nf,1);vals=dict(row['ideal_factorization']);N=ZZ(pari.nfeltnorm(nf,alpha));assert N
        for i,e in vals.items():
            assert e>0;I=pari.idealmul(nf,I,pari.idealpow(nf,ideals[i],e))
        assert pari.idealhnf(nf,I)==pari.idealhnf(nf,alpha)
        assert prod(cols[i]['p']**(cols[i]['f']*e) for i,e in vals.items())==abs(N)
        mask=0
        for p in {cols[i]['p'] for i in vals}-S:
            bit=sum(cols[i]['f']*vals.get(i,0) for i in blocks[p])%2
            for i in blocks[p]:mask^=((vals.get(i,0)+bit*cols[i]['e'])%2)<<i
        masks.append(mask);norms.append(str(N))
    rank,kernel=elimination(masks);phases=[]
    for phase in d['phases']:
        n=phase['cumulative_unique_elements'];rr,kk=elimination(masks[:n])
        phases.append({**phase,'outside_S_parity_rank':rr,'kernel_dimension':len(kk)})
    assert len(kernel)==1 and kernel[0].bit_count()==1
    index=kernel[0].bit_length()-1;row=d['relations'][index];known=r.read(hit.OUTPUT)
    assert [row['m'],row['n']]==known['address'] and row['alpha_ascending']==known['alpha_ascending']
    assert norms[index]==known['norm_alpha'] and known['additional_squareclass_dimension']==0
    # The independent checker below replays the stored square-root identity too.
    return {'schema':'rank-jump.early-relation-pool.v1','status':'PASS','bindings':bindings(),
        'element_count':len(masks),'complete_prime_blocks':len(blocks),'ideal_column_count':len(cols),
        'outside_S_parity_rank':rank,'kernel_supports':[[index]],'kernel_addresses':[[row['m'],row['n']]],
        'norms_alpha':norms,'phases':phases,'additional_strict_class_capacity':0,
        'inherited_generic_mask':known['generic_product_mask'],
        'boundary':'All coefficient-kernel generators are inherited; the whole retained norm-projection dictionary adds no strict class modulo G. No full class-group or rational-solubility upper bound.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        error=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=60)
                if p.returncode:error='Worker failure'
            except subprocess.TimeoutExpired:error='Bounded worker timeout'
        if error:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':error})
    value=r.read(path);assert value['bindings']==bindings();r.write_new(OUTPUT,value)
    print({k:v for k,v in value.items() if k not in ('bindings','norms_alpha','phases')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/'worker.json',compute())
    else:globals()[args.mode]()

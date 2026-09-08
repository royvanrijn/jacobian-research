#!/usr/bin/env python3
"""Disjoint-prime diagnostics on a fixed, previously prepared population.

Freeze natural prime prefixes and two disjoint validation blocks before any
new trace computation. Shared residue-equation tables amortize point counts.
No parameter is added and no validation statistic is a rank prediction.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from math import log, sqrt
from pathlib import Path
import statistics
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint,digest
from mod2_reduction_independence import _is_prime

ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/bounded-prime-selectors-v1'
POP=ART/'fibre_height_population_v1.json.gz'
SCALE=10**12


def read(path):
    data=path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix=='.gz' else data)


def short_mod(model,p):
    a1,a2,a3,a4,a6=(int(a)%p for a in model)
    inv2,inv3=pow(2,-1,p),pow(3,-1,p)
    c2=(a2+a1*a1*inv2**2)%p
    c4=(a4+a1*a3*inv2)%p
    c6=(a6+a3*a3*inv2**2)%p
    return (c4-c2*c2*inv3)%p,(c6-c2*c4*inv3+2*c2**3*pow(27,-1,p))%p


def trace(a,b,p,chi):
    if (4*a**3+27*b*b)%p==0:return None
    return -sum(chi[(x*x*x+a*x+b)%p] for x in range(p))


def contributions(t,p):
    if t is None:return (0,0)
    n=p+1-t
    return round((2-t)/n*log(p)*SCALE),round(log(n/p)*SCALE)


def ranks(values):
    ordered=sorted(set(values)); counts=Counter(values); start=0; out={}
    for value in ordered:
        out[value]=start+(counts[value]-1)/2;start+=counts[value]
    return [out[v] for v in values]


def correlation(x,y):
    if len(x)<2:return None
    x,y=ranks(x),ranks(y);mx,my=statistics.mean(x),statistics.mean(y)
    xx=sum((v-mx)**2 for v in x);yy=sum((v-my)**2 for v in y)
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/sqrt(xx*yy) if xx*yy else None


def freeze():
    path=LOCAL/'protocol.json'
    if path.exists():raise FileExistsError(path)
    pop=read(POP);primes=[p for p in range(5,3000) if _is_prime(p)][:384]
    assert len(primes)==384
    panel={r['id']:dict(family=f,model=r['minimal_model']) for f,rows in pop['families'].items() for r in rows['rows']}
    split={}
    for family in pop['families']:
        ids=sorted((i for i,r in panel.items() if r['family']==family),key=lambda i:digest(['prime-selector-split-v1',i]))
        split.update({i:'development' if k<len(ids)//2 else 'holdout' for k,i in enumerate(ids)})
    p=dict(schema='elliptic-curves.bounded-prime-selector-protocol.v1',population_sha256=sha256(POP.read_bytes()).hexdigest(),
        panel=panel,split=split,prefixes=[25,64,128,256],training_primes=primes[:256],
        validation_blocks=[primes[256:320],primes[320:384]],scores=['nagao','mestre_log'],
        score_formulas=['sum_good ((2-a_p)/(p+1-a_p))*log(p)','sum_good log((p+1-a_p)/p)'],
        bad_prime_rule='Omit singular reductions of the retained global minimal model, record all counts.',
        quantization='Round each prime contribution to nearest 1e-12, then sum integers.',
        top_per_full_family=8,top_per_split=4,top_per_measured_family=4,
        recovery_panel='The fixed 31 curves of the completed fibre-height pilot, with equal masked exposure. '
        'Separately report its original full-subgroup search outcomes. Never impute outcomes for unsearched curves.',
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        claim_boundary='Fixed finite comparison, no policy promotion or true-rank enrichment claim. '
        'Validation primes are disjoint from selection primes; curve holdouts share the same two families.')
    p['protocol_hash']=digest(p);checkpoint(path,p)
    print('FROZEN',len(panel),len(primes),p['protocol_hash'],flush=True)


def protocol():
    p=read(LOCAL/'protocol.json');h=p.pop('protocol_hash');assert digest(p)==h;p['protocol_hash']=h
    assert p['population_sha256']==sha256(POP.read_bytes()).hexdigest()
    assert p['source_sha256']==sha256(Path(__file__).read_bytes()).hexdigest()
    assert not set(p['training_primes'])&set(sum(p['validation_blocks'],[]))
    return p


def compute():
    p=protocol();start=time.monotonic(); tables={}
    for prime in p['training_primes']+sum(p['validation_blocks'],[]):
        path=LOCAL/'traces'/f'{prime}.json'
        if path.exists():table=read(path)
        else:
            chi=[0]*prime
            for x in range(1,prime):chi[x]=-1
            for x in range(1,prime):chi[x*x%prime]=1
            entries={}; lookup={}
            for identifier,row in p['panel'].items():
                a,b=short_mod(row['model'],prime);key=f'{a},{b}'
                if key not in entries:entries[key]=trace(a,b,prime,chi)
                lookup[identifier]=key
            table=dict(protocol_hash=p['protocol_hash'],prime=prime,traces=entries,candidate_lookup=lookup)
            checkpoint(path,table)
        if table['protocol_hash']!=p['protocol_hash'] or table['prime']!=prime:raise ArithmeticError('trace checkpoint drift')
        tables[str(prime)]=table
        if len(tables)%32==0:print('PRIMES',len(tables),'wall',round(time.monotonic()-start,2),flush=True)
    output=dict(protocol=p,tables=tables,unique_point_counts=sum(len(t['traces']) for t in tables.values()),
        candidate_prime_lookups=len(p['panel'])*len(tables),elapsed_seconds=time.monotonic()-start)
    raw=(json.dumps(output,sort_keys=True,separators=(',',':'))+'\n').encode()
    (ART/'bounded_prime_selector_traces_v1.json.gz').write_bytes(gzip.compress(raw,mtime=0))


def analyze():
    data=read(ART/'bounded_prime_selector_traces_v1.json.gz');p=protocol();assert data['protocol']==p
    old=read(ART/'fibre_height_search_v1.json.gz')['results']
    masked=read(ART/'ordinary_masked_controls_v1.json.gz')['results']
    pop=read(POP);rows={r['id']:r for f in pop['families'].values() for r in f['rows']}
    values={}
    def score(identifier,primes,kind):
        tracevals=[(q,data['tables'][str(q)]['traces'][data['tables'][str(q)]['candidate_lookup'][identifier]]) for q in primes]
        return sum(contributions(t,q)[kind] for q,t in tracevals),sum(t is not None for q,t in tracevals)
    for identifier in p['panel']:
        values[identifier]=dict(family=p['panel'][identifier]['family'],split=p['split'][identifier],scores={},validation={})
        for kind,name in enumerate(p['scores']):
            values[identifier]['scores'][name]={str(k):score(identifier,p['training_primes'][:k],kind) for k in p['prefixes']}
            values[identifier]['validation'][name]=[score(identifier,b,kind) for b in p['validation_blocks']]
    comparisons=[]
    for family in pop['families']:
        for subset in ('all','development','holdout','measured'):
            ids=[i for i,v in values.items() if v['family']==family and (subset=='all' or subset=='measured' and i in masked or v['split']==subset)]
            top=p['top_per_full_family'] if subset=='all' else p['top_per_split']
            for name in p['scores']:
                val=[sum(t[0] for t in values[i]['validation'][name]) for i in ids]
                vr=dict(zip(ids,ranks(val)))
                for prefix in p['prefixes']:
                    selected=sorted(ids,key=lambda i:(-values[i]['scores'][name][str(prefix)][0],i))[:top]
                    r=dict(family=family,subset=subset,n=len(ids),score=name,prefix=prefix,selected=selected,
                        spearman_validation=correlation([values[i]['scores'][name][str(prefix)][0] for i in ids],val),
                        selected_mean_validation_percentile=statistics.mean(vr[i]/(len(ids)-1) for i in selected),
                        selected_validation_mean_units=statistics.mean(sum(t[0] for t in values[i]['validation'][name]) for i in selected),
                        population_validation_mean_units=statistics.mean(val),
                        selected_mean_coefficient_bits=statistics.mean(rows[i]['minimal_maximum_coefficient_bits'] for i in selected),
                        validation_block_correlations=[correlation([values[i]['scores'][name][str(prefix)][0] for i in ids],
                            [values[i]['validation'][name][b][0] for i in ids]) for b in (0,1)])
                    if subset=='measured':
                        r.update(masked_direct_signed_recoveries=sum(masked[i]['direct_signed_recovery'] for i in selected),
                            masked_finite_points=sum(masked[i]['finite_point_count'] for i in selected),
                            masked_completed_charts=sum(c['status']=='bounded_search_complete' for i in selected for c in masked[i]['charts']),
                            masked_allocated_charts=12*len(selected),
                            original_certified_new_directions=sum(old[i]['certified_new_directions'] for i in selected),
                            original_completed_charts=sum(c['status']=='bounded_search_complete' for i in selected for c in old[i]['charts']))
                    comparisons.append(r)
    out=dict(schema='elliptic-curves.bounded-prime-selector-comparison.v1',protocol=p,values=values,comparisons=comparisons,
        input_hashes={str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in
            (ART/'bounded_prime_selector_traces_v1.json.gz',ART/'ordinary_masked_controls_v1.json.gz',ART/'fibre_height_search_v1.json.gz')},
        unique_point_counts=data['unique_point_counts'],candidate_prime_lookups=data['candidate_prime_lookups'],
        status='COMPLETE_DIAGNOSTIC_NO_PROMOTION',claim_boundary=p['claim_boundary'])
    checkpoint(ART/'bounded_prime_selector_comparison_v1.json',out)
    for r in comparisons:
        if r['subset']=='holdout':print(r['family'],r['score'],r['prefix'],round(r['spearman_validation'],3),round(r['selected_mean_validation_percentile'],3),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=('freeze','compute','analyze'))
    a=parser.parse_args();{'freeze':freeze,'compute':compute,'analyze':analyze}[a.mode]()

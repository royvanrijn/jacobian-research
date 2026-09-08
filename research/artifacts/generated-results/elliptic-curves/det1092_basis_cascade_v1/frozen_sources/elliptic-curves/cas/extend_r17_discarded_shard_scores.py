#!/usr/bin/env python3
"""Extend primes on768 already saved sign-shard finalists discarded by the global short-prime cut."""
import argparse
from pathlib import Path
from math import log
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import select_compact_r17_wide as original
import benchmark_r17_extended_prime_traces as engine
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];D=engine.LOCAL/'r17-discarded-shards-extended-v1';PARENT=engine.PARENT

def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(original.__file__).resolve(),ROOT/'elliptic-curves/tests/test_r17_prime_extension.py')}}

def sums(traces):
    low=high=low_good=high_good=0
    for p,t in traces:
        if t is None:continue
        value=round((2-t)*log(p)/(p+1-t)*10**12)
        if p<=32749:low+=value;low_good+=1
        else:high+=value;high_good+=1
    return {'extension_selection_units':low,'validation_units':high,'extension_good':low_good,'validation_good':high_good}

def choose(rows):
    # No validation values or measured ranks enter this function.
    remaining=list(rows)
    ordered=sorted(remaining,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
    return {'extended_top_two':[r['retained_index'] for r in ordered[:2]]}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve extension protocol')
    bench=cert.read(engine.D/'result.json')
    if bench['status']!='PASS' or len(bench['rows'])!=6:raise ArithmeticError('trace calibration incomplete')
    prior=cert.read(ROOT/'artifacts/generated-results/elliptic-curves/next24_portable_replay_v1.json')
    if prior['status']!='PASS':raise ArithmeticError('completed prior experiment gate failed')
    rows=[];pops={};shards={};parent=cert.read(PARENT/'protocol.json')
    for f in cert.read(spec.ATLAS)['families']:
        family=f['family'];path=PARENT/family/'population.json';population=cert.read(path);pops[family]=cert.hashed(path);discarded=[]
        seen={r['parameter'] for r in population['retained_candidates']}
        for sign in (-1,1):
            sp=PARENT/family/f'scan-{sign}.json';lp=PARENT/family/f'scan-{sign}.log';sd=cert.read(sp);parsed,summary=original.parse(lp.read_text(),sign)
            if parsed!=sd['rows'] or summary!=sd['summary'] or sd['protocol_hash']!=digest(parent) or sd['supervision']['outcome']!='completed' or sd['supervision']['returncode']!=0 or sd['supervision']['log_sha256']!=cert.hashed(lp):raise ArithmeticError('saved sign shard differs')
            shards[str(sp.relative_to(ROOT))]=cert.hashed(sp);shards[str(lp.relative_to(ROOT))]=cert.hashed(lp)
            discarded.extend(r for r in sd['rows'] if r['parameter'] not in seen)
        discarded.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
        if len(discarded)!=128 or len({r['parameter'] for r in discarded})!=128:raise ArithmeticError('fixed discarded128 roster differs')
        for index,r in enumerate(discarded):
            model,pts=spec.specialize(f,r['parameter'])
            if any(q.denominator!=1 for q in model):raise ArithmeticError('integral specialization required')
            rows.append({**r,'family':family,'retained_index':index,'model':list(map(str,model))})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-discarded-shards-extended.v1','sources':sources(),'rows':rows,'population_hashes':pops,'shard_hashes':shards,'prior_completed_replay_sha256':cert.hashed(ROOT/'artifacts/generated-results/elliptic-curves/next24_portable_replay_v1.json'),'benchmark_sha256':cert.hashed(engine.D/'result.json'),'gp_sha256':cert.hashed(engine.GP),'prime_roster':engine.PRIMES,'selection_extension':[4099,32749],'disjoint_validation':[32771,65521],'per_curve_seconds':2,'outer_seconds':600,'rss_bytes':1073741824,'maximum_workers':1,'score':'Previously declared quantized562-prime score plus the same selection extension through32749. Keep32771..65521 exclusively for validation. Same good-prime and denominator/numerator ties.','finite_gate':'The previous extended selector produced three catalogue-unmatched rank27 curves, but only re-ranked128 finalists per family. Existing complete scans retained128 per sign before merging. Reuse the128 additional saved addresses per family without another122million-address scan. No claim that this recovers all discarded candidates.','selection':'Exactly768 addresses absent from the previous768 extended pool. Freeze the best two per family by the extended selection score only. No catalogue, public point, record equation or validation score enters selection.','prospective_point_scope':'No point search in this script. At most12 candidates may enter a separately frozen, equally exposed point protocol after exact trace/selection replay and the current height-control outcome. No automatic continuation.','boundaries':'Saved sign-shard finalists remain a truncated population. A better finite score is not a rank theorem. This extension tests omitted saved candidates, not a new broad parameter campaign.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP):raise ArithmeticError('extension sources changed')
    for f,h in p['population_hashes'].items():
        if h!=cert.hashed(PARENT/f/'population.json'):raise ArithmeticError('population changed')
    for name,h in p['shard_hashes'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('saved shard changed')
    return p

def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve finite extension attempt')
    data={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,data)
    for r in p['rows']:
        folder=D/r['family']/f"candidate-{r['retained_index']:03}";cmd=engine.program(r['model']);cap=capture([str(engine.GP),'-q','-s','256000000'],input_text=cmd,limits=Limits(p['per_curve_seconds'],536870912),log_path=folder/'gp.log',separate_stderr=True,check=False)
        raw={'program':cmd,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision};checkpoint(folder/'raw.json',raw)
        if cap.supervision['outcome']!='completed' or cap.supervision['returncode']!=0 or cap.stderr:raise ArithmeticError('trace failure/censor, retain partial batch')
        traces,ms=engine.parse(cap.stdout,r['model']);s=sums(traces);row={**r,**s,'combined_selection_units':r['score_units']+s['extension_selection_units'],'combined_good':r['good_primes']+s['extension_good'],'trace_cpu_ms':ms,'raw_sha256':cert.hashed(folder/'raw.json')};data['rows'].append(row);checkpoint(out,data)
        if len(data['rows'])%32==0:print('EXTENDED DISCARDED R17',len(data['rows']),'of768',flush=True)
    data['selection']={f:choose([r for r in data['rows'] if r['family']==f]) for f in p['population_hashes']};data['status']='COMPLETE_FROZEN_TRACE_EXTENSION';checkpoint(out,data)

def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['protocol_hash']!=digest(p) or d['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or len(d['rows'])!=len(p['rows']):raise ArithmeticError('incomplete trace roster')
    for original,r in zip(p['rows'],d['rows']):
        if any(r[k]!=v for k,v in original.items()):raise ArithmeticError('input roster differs')
        path=D/r['family']/f"candidate-{r['retained_index']:03}"/'raw.json';raw=cert.read(path)
        if cert.hashed(path)!=r['raw_sha256'] or raw['program']!=engine.program(r['model']) or raw['stderr'] or raw['supervision']['returncode']!=0 or raw['supervision']['outcome']!='completed':raise ArithmeticError('trace source failed')
        traces,ms=engine.parse(raw['stdout'],r['model']);s=sums(traces)
        if any(r[k]!=v for k,v in s.items()) or r['combined_selection_units']!=r['score_units']+s['extension_selection_units'] or r['combined_good']!=r['good_primes']+s['extension_good'] or r['trace_cpu_ms']!=ms:raise ArithmeticError('score differs')
    selection={f:choose([r for r in d['rows'] if r['family']==f]) for f in p['population_hashes']}
    if selection!=d['selection']:raise ArithmeticError('final selection differs')
    print('REPLAYED768 DISCARDED TRACE ROSTERS AND FIXED SCORE SELECTION',selection,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();globals()[a.stage]()

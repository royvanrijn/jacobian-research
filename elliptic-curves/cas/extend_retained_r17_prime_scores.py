#!/usr/bin/env python3
"""Frozen disjoint-prime extension of the768 retained compact R17 addresses."""
import argparse
from pathlib import Path
from math import log
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import benchmark_r17_extended_prime_traces as engine
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];D=engine.LOCAL/'r17-retained-extended-primes-v1';PARENT=engine.PARENT

def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/tests/test_r17_prime_extension.py')}}

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
    remaining=[r for r in rows if r['retained_index']>=4]
    ordered=sorted(remaining,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
    return {'extended_top_two':[r['retained_index'] for r in ordered[:2]],'original_next_two':[4,5]}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve extension protocol')
    bench=cert.read(engine.D/'result.json')
    if bench['status']!='PASS' or len(bench['rows'])!=6:raise ArithmeticError('trace calibration incomplete')
    rows=[];pops={}
    for f in cert.read(spec.ATLAS)['families']:
        path=PARENT/f['family']/'population.json';population=cert.read(path);pops[f['family']]=cert.hashed(path)
        if len(population['retained_candidates'])!=128:raise ArithmeticError('retained population differs')
        for index,r in enumerate(population['retained_candidates']):
            model,pts=spec.specialize(f,r['parameter'])
            if any(q.denominator!=1 for q in model):raise ArithmeticError('integral specialization required')
            rows.append({**r,'family':f['family'],'retained_index':index,'model':list(map(str,model))})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-retained-extended-primes.v1','sources':sources(),'rows':rows,'population_hashes':pops,'benchmark_sha256':cert.hashed(engine.D/'result.json'),'gp_sha256':cert.hashed(engine.GP),'prime_roster':engine.PRIMES,'selection_extension':[4099,32749],'disjoint_validation':[32771,65521],'per_curve_seconds':2,'outer_seconds':600,'rss_bytes':1073741824,'maximum_workers':1,'score':'Retained quantized all562-prime score plus sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)) over nonsingular displayed reductions at4099..32749. Ties: combined good-prime count descending, denominator ascending, numerator ascending. Validation uses only32771..65521.','finite_gate':'Allsix predetermined benchmark curves completed5978 extra prime traces each with48 direct character-sum checks. Earlier short-prime selection lost high-rank controls; test extension on the fixed768 retained addresses, not a new parameter sweep.','selection':'Exclude indices0..3, which received the prior fixed batch. Per family select two by combined score and separately retain original next indices4,5 as a paired baseline. Shared addresses receive one point attempt, no refill. Selection never consumes validation scores, catalogue data, public points or measured ranks.','prospective_point_scope':'Freeze a separate point protocol after all768 trace records and exact score replay pass. At most24 distinct family-address pairs, all exact maximum generic parity classes43or49 at H100000 with PARI maps. No point campaign launched by this script.','boundaries':'The768 pool is already truncated by the earlier4093-prime score. No extension recovers candidates outside that pool. Numerical scores schedule finite work, not rank claims. Previously measured24 may be analyzed only after all scores/selection are frozen; no policy promotion from their lower bounds.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP):raise ArithmeticError('extension sources changed')
    for f,h in p['population_hashes'].items():
        if h!=cert.hashed(PARENT/f/'population.json'):raise ArithmeticError('population changed')
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
        if len(data['rows'])%32==0:print('EXTENDED RETAINED R17',len(data['rows']),'of768',flush=True)
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
    print('REPLAYED768 EXTENDED TRACE ROSTERS AND FIXED SCORE SELECTION',selection,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();globals()[a.stage]()

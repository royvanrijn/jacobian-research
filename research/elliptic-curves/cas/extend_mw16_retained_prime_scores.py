#!/usr/bin/env python3
"""Extended-prime selection on all1280 already saved MW16 sign-shard finalists."""
import argparse
from pathlib import Path
from math import log,gcd
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
import benchmark_mw16_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'mw16-retained-extended-primes-v1';PARENT=LOCAL/'prospective-mw16-h4096-v1';TABLES=LOCAL/'prospective-mw16-h1024-v1'
def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(scoring.__file__).resolve(),spec.ATLAS)}}
def short_rows(family):
    rows=[];bindings={}
    for sign in (-1,1):
        path=PARENT/family/f'scan-{sign}.json';raw=cert.read(path)
        if len(raw['rows'])!=128 or raw['summary']!=[4096,4096,128,10200039,10200039,128]:raise ArithmeticError('fixed saved MW16 shard differs')
        if any(gcd(abs(r['numerator']),r['denominator'])!=1 or r['numerator']*sign<=0 or not 1<=r['denominator']<=4096 or abs(r['numerator'])>4096 for r in raw['rows']):raise ArithmeticError('primitive saved address differs')
        rows+=raw['rows'];bindings[str(path.relative_to(ROOT))]=cert.hashed(path)
    rows.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
    if len({r['parameter'] for r in rows})!=256:raise ArithmeticError('saved family addresses not unique')
    return rows,bindings
def verify_short(family,rows):
    f=next(r for r in cert.read(spec.ATLAS)['families'] if r['fibration_id']==family);model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]};scores=[0]*256;goods=[0]*256;hashes={}
    for p in (p for p in range(5,4094) if _is_prime(p)):
        path=TABLES/family/'trace-tables'/f'{p}.json';d=cert.read(path)
        if d['input']!={'family':family,'model_hash':digest(model),'prime':p}:raise ArithmeticError('canonical MW16 table differs')
        hashes[str(path.relative_to(ROOT))]=cert.hashed(path)
        for i,r in enumerate(rows):
            t=r['numerator']*pow(r['denominator'],-1,p)%p if r['denominator']%p else p
            if d['good'][t]:
                ap=d['traces'][t];scores[i]+=round((2-ap)/(p+1-ap)*log(p)*10**12);goods[i]+=1
    if any((r['score_units'],r['good_primes'])!=(s,g) for r,s,g in zip(rows,scores,goods)):raise ArithmeticError('all562 short scores differ')
    return hashes
def model_at(f,t):
    t=cert.F(t);d=t.denominator;model=[cert.F(0)]*3+[spec.polynomial(f['A_coefficients_low_to_high'],t)*d**8,spec.polynomial(f['B_coefficients_low_to_high'],t)*d**12]
    if any(q.denominator!=1 for q in model) or not 4*model[3]**3+27*model[4]**2:raise ArithmeticError('integral nonsingular specialization required')
    return list(map(str,model))
def choose(rows):
    ordered=sorted(rows,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
    return [r['retained_index'] for r in ordered[:4]]
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve MW16 prime protocol')
    b=cert.read(engine.D/'result.json');replay=cert.read(engine.D/'replay.json');coordinate=cert.read(LOCAL/'mw16-top25-pari-followup-v1/ledger.json')
    if b['status']!='PASS' or replay['status']!='PASS' or coordinate['status']!='PASS':raise ArithmeticError('trace/control gates incomplete')
    families=cert.read(spec.ATLAS)['families'];rows=[];shards={};tables={}
    for f in families:
        family=f['fibration_id'];pool,h=short_rows(family);shards.update(h);tables.update(verify_short(family,pool))
        for i,r in enumerate(pool):rows.append({**r,'family':family,'retained_index':i,'model':model_at(f,r['parameter'])})
    if len(rows)!=1280:raise ArithmeticError('fixed1280 population differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-retained-prime-extension.v1','sources':sources(),'rows':rows,'shard_hashes':shards,'canonical_table_hashes':tables,'benchmark_sha256':cert.hashed(engine.D/'result.json'),'benchmark_replay_sha256':cert.hashed(engine.D/'replay.json'),'coordinate_gate_sha256':cert.hashed(LOCAL/'mw16-top25-pari-followup-v1/ledger.json'),'prime_roster':engine.PRIMES,'gp_sha256':cert.hashed(engine.GP),'seconds_per_curve':5,'rss_bytes_per_curve':536870912,'outer_seconds':600,'maximum_workers':1,'checkpoint_block_size':16,'gate':'All five compactMW16 families retain256 saved signed H4096 finalists. Their prior selectors used only562 primes and the old incomplete metric-coordinate point boxes. The five-address extended-trace benchmark passes40 exact character-sum checks; both existing new25 curves now complete and replay all43 PARI-coordinate boxes at unchanged100000 height. Extend all1280 saved addresses without a new broad scan or prior-rank/catalogue filtering.','selection':'Four per family by original quantized short score plus extension through32749, then combined selection-band good count, denominator, signed numerator. Disjoint validation32771..65521 never enters ordering. No known-record targets, public points, prior ranks or catalogue labels enter the trace worker or selection.','future_point_scope':'At most20 equally exposed generic16-only point attempts, using the same43 recorded generic parity labels per family at100000/10seconds, after complete score replay and a separate frozen point protocol. These43 computed-class choices are not claimed to exhaust every generic or specialized parity maximum.','checkpoints':'Immutable raw per-address calls and aggregate16-address checkpoints. Failed/censored calls remain; no automatic retry. Explicit resume may reuse only source-bound successful records.','boundaries':'Finite short-prime-truncated saved population; no global score optimum, rank prediction theorem, upper rank bound, exact rank or universal novelty.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP):raise ArithmeticError('MW16 trace inputs changed')
    for name,h in p['shard_hashes'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('saved shard changed')
    return p
def checked_row(row,p,create):
    folder=D/row['family']/f"candidate-{row['retained_index']:03}";path=folder/'raw.json';program=engine.program(row['model'])
    if not path.exists():
        if not create:raise ArithmeticError('missing raw trace')
        cap=capture([str(engine.GP),'-q','-s','256000000'],input_text=program,limits=Limits(p['seconds_per_curve'],p['rss_bytes_per_curve']),log_path=folder/'gp.log',separate_stderr=True,check=False);checkpoint(path,{'program':program,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision})
    raw=cert.read(path)
    if raw['program']!=program or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('failed/censored immutable trace; no automatic retry')
    traces,ms=engine.parse(raw['stdout'],row['model']);scores=scoring.sums(traces)
    return {**row,**scores,'combined_selection_units':row['score_units']+scores['extension_selection_units'],'combined_good':row['good_primes']+scores['extension_good'],'cpu_ms':ms,'raw_sha256':cert.hashed(path)}
def run(resume=False):
    p=protocol();out=D/'result.json'
    if out.exists():
        if not resume:raise FileExistsError('explicit resume required')
        data=cert.read(out)
        if data['status']!='RUNNING' or data['protocol_hash']!=digest(p):raise ArithmeticError('invalid resume checkpoint')
        for r,s in zip(p['rows'],data['rows']):
            if checked_row(r,p,False)!=s:raise ArithmeticError('committed row differs')
    else:
        if resume:raise FileNotFoundError('no checkpoint')
        data={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,data)
    for start in range(len(data['rows']),1280,16):
        data['rows'] += [checked_row(r,p,True) for r in p['rows'][start:start+16]];checkpoint(out,data)
        if len(data['rows'])%128==0:print('MW16 EXTENDED',len(data['rows']),'of1280',flush=True)
    data['selection']={f:choose([r for r in data['rows'] if r['family']==f]) for f in sorted({r['family'] for r in p['rows']})};data['status']='COMPLETE_FROZEN_TRACE_EXTENSION';checkpoint(out,data)
def replay():
    p=protocol();d=cert.read(D/'result.json');families={f['fibration_id']:f for f in cert.read(spec.ATLAS)['families']};expected=[];tables={}
    for family in families:
        rows,_=short_rows(family);tables.update(verify_short(family,rows));expected += [{**r,'family':family,'retained_index':i,'model':model_at(families[family],r['parameter'])} for i,r in enumerate(rows)]
    if expected!=p['rows'] or tables!=p['canonical_table_hashes'] or d['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or d['protocol_hash']!=digest(p) or len(d['rows'])!=1280:raise ArithmeticError('full1280 trace roster differs')
    for row,r in zip(p['rows'],d['rows']):
        if checked_row(row,p,False)!=r:raise ArithmeticError('trace/score replay differs')
    selection={f:choose([r for r in d['rows'] if r['family']==f]) for f in sorted(families)}
    if selection!=d['selection']:raise ArithmeticError('fixed20 selector differs')
    print('REPLAYED1280 SHORT SCORES, EXTENDED TRACE ROSTERS AND FIXED20 SELECTION',selection,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','resume','replay']);a=p.parse_args();run(True) if a.stage=='resume' else globals()[a.stage]()

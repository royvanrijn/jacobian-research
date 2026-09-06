#!/usr/bin/env python3
"""Eight frozen signed slices beyond the completely scanned11952 square."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import scan_mw16_outer_bands as annulus
import benchmark_periodic_nagao_scanner as r17
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
ROOT=annulus.ROOT;CAS=annulus.CAS;ART=annulus.ART
D=annulus.LOCAL/'11952-new-annulus-v1';INNER=131072;OUTER=524288;MODULUS=16384;KEEP=4096
SALT='11952-new-annulus-beyond-complete131072-v1'
PRIOR=annulus.LOCAL/'full11952-h131072-short-v1/protocol.json'
RESULT=ART/'full11952_late64_experiment_v1.json'
PORTABLE=ART/'full11952_late64_point_portable_replay_v1.json'

def sources():
    paths=[Path(__file__).resolve(),Path(annulus.__file__),Path(r17.__file__),annulus.BINARY,annulus.REFERENCE,annulus.CPP,r17.spec.ATLAS,PRIOR,RESULT,PORTABLE,annulus.D/'benchmark.json',annulus.D/'replay.json',Path(cert.__file__),CAS/'verify_periodic_nagao_scanner.py',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def roster():
    prior=cert.read(PRIOR);rows=[]
    for si,sign in enumerate((-1,1)):
        table=next(r['table'] for r in prior['rows'] if r['sign']==sign)
        for quarter in range(4):
            h=sha256(f'{SALT}|{sign}|{quarter}'.encode()).hexdigest()
            shard=quarter*4096+2*(int(h,16)%2048)+(si+quarter)%2
            count=annulus.counts.population(OUTER,OUTER,shard,MODULUS)-annulus.counts.population(INNER,INNER,shard,MODULUS)
            rows.append({'id':('neg' if sign<0 else 'pos')+f'-q{quarter}','family':'11952','sign':sign,'quarter':quarter,'inner':INNER,'outer':OUTER,'shard':shard,'shards':MODULUS,'sha256_choice':h,'primitive_population':count,'table':table})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve new annulus')
    prior=cert.read(PRIOR)
    if prior['numerator_bound']!=INNER or prior['denominator_bound']!=INNER or prior['total_primitive_population']!=20888422894:raise ArithmeticError('previous complete square differs')
    if cert.read(RESULT)['status']!='PASS' or cert.read(PORTABLE)['status']!='PASS' or cert.read(PORTABLE)['logical_stages']!=194:raise ArithmeticError('completed late11952 point proofs required')
    if cert.read(annulus.D/'benchmark.json')['status']!='PASS' or cert.read(annulus.D/'replay.json')['status']!='PASS':raise ArithmeticError('existing annular binary proof gates required')
    rows=roster()
    if len(rows)!=8 or len({(r['sign'],r['shard']) for r in rows})!=8 or sum(r['shard']%2 for r in rows)!=4:raise ArithmeticError('eight distinct balanced slices required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.11952-new-annulus.v1','sources':sources(),'rows':rows,'keep_per_slice':KEEP,'maximum_workers':2,'seconds_per_slice':60,'first_slice_cost_gate_seconds':25,'outer_seconds':600,'total_primitive_population':sum(r['primitive_population'] for r in rows),'gate':'The completed two11952 point cohorts each supply a new independent27-point curve from generic17-only charts. Their parameters lie inside the now completely short-scanned131072 square. A new annulus tests candidate incidence beyond that coverage; no presumption of higher rank, density or point visibility at greater coefficient size is made. Reuse the proven annular scanner only after new signed11952 frame regressions.','population':'131072<max(abs(n),d)<=524288, nonzero signed primitive n/d. Four denominator slices per sign modulo16384, one SHA256 choice in each quarter of the residue interval, with balanced parity. Every address lies outside the former complete square. Neither known targets nor previously successful parameter neighbourhoods define the slices. Parameter disjointness alone does not prove equation novelty.','retention':'Keep4096 per signed slice under unchanged562-prime quantized S1, good count, denominator and absolute numerator. All32768 survivors pass canonical projective-table score replay. No compact-square rescan, denominator refill or retry.','validation':'The same annular binary already passes40 MW16 frames and their full retained-score replay. Before this population, check complete reference and top7 annular ordering on all eight actual signed11952 large-modulus frames plus both complete small squares. Recompute every returned score from the11952 canonical projective trace tables and independently count primitive addresses. First fullslice has a25-second runtime gate and is reused once.','future_scope':'Only after full scanner and score replay may a separate protocol use the existing exact4099..32749 cache on all32768 survivors, freeze4096 distinct equations for scalar scoring through65521, and choose64 finalists with wholly disjoint65537..131071 validation. A later separate generic17-only point protocol may expose49 charts per fibre at125000/10seconds, all maps before points. No point search is launched here.','boundaries':'One bounded stratified annulus experiment, not complete coverage, a global optimum, a rank prediction, point absence or universal novelty. Low detection at higher height remains distinct from low incidence of rank jumps.'})
    print('FROZEN8 NEW11952 SLICES',sum(r['primitive_population'] for r in rows),'ADDRESSES32768 RETENTION',flush=True)

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen annulus inputs changed')
    for row in p['rows']:
        if cert.hashed(ROOT/row['table']['path'])!=row['table']['sha256']:raise ArithmeticError('signed table changed')
    return p

def exact_scores(rows):return r17.exact_scores(rows,{'family':'11952'})

def benchmark(check=False):
    p=protocol();records=[]
    for row in p['rows']:
        cases=[(31,2*MODULUS,MODULUS,row['shard'],MODULUS)]
        if row['quarter']==0:cases.append((31,29,7,0,1))
        for ci,(N,M,inner,shard,modulus) in enumerate(cases):
            folder=D/'benchmark'/row['id']/str(ci)
            full,_=annulus.invocation(row,N,M,0,1000,shard,modulus,folder/'reference.json',not check,True)
            expected={(row['sign']*n,d) for n in range(1,N+1) for d in range(1+shard,M+1,modulus) if annulus.gcd(n,d)==1}
            if {(r['numerator'],r['denominator']) for r in full}!=expected:raise ArithmeticError('complete small reference population differs')
            filtered=[r for r in full if max(abs(r['numerator']),r['denominator'])>inner]
            for keep in (1000,7):
                rows,_=annulus.invocation(row,N,M,inner,keep,shard,modulus,folder/f'annulus-{keep}.json',not check)
                if rows!=filtered[:keep]:raise ArithmeticError('new11952 annulus/top7 differs')
            tables=exact_scores(full)
            records.append({'id':row['id'],'case':[N,M,inner,shard,modulus],'population':len(full),'canonical_tables_hash':digest(tables),'raw_hashes':{q.name:cert.hashed(q) for q in sorted(folder.glob('*.json'))}})
    result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'records':records}
    if check:
        if cert.read(D/'benchmark.json')!=result:raise ArithmeticError('new11952 benchmark replay differs')
    else:
        if (D/'benchmark.json').exists():raise FileExistsError('preserve new11952 benchmark')
        checkpoint(D/'benchmark.json',result)
    print('NEW11952 ALL10 FRAME REGRESSIONS PASS',flush=True)

def checked(row,p,create):
    path=D/row['id']/'raw.json'
    rows,raw=annulus.invocation(row,OUTER,OUTER,INNER,KEEP,row['shard'],MODULUS,path,create)
    if len(rows)!=KEEP:raise ArithmeticError('full4096 retention required')
    return {**row,'rows':rows,'raw_sha256':cert.hashed(path),'canonical_hash':digest(exact_scores(rows)),'wall_seconds':raw['supervision']['wall_seconds']}

def run():
    p=protocol();benchmark(True)
    if (D/'result.json').exists():raise FileExistsError('preserve new11952 run')
    data={'status':'RUNNING','protocol_hash':digest(p),'shards':[]};checkpoint(D/'result.json',data);got={}
    try:
        got[0]=checked(p['rows'][0],p,True);data['shards']=[got[0]];checkpoint(D/'result.json',data)
        if got[0]['wall_seconds']>p['first_slice_cost_gate_seconds']:raise ArithmeticError('first new11952 slice runtime gate failed')
        with ThreadPoolExecutor(max_workers=2) as pool:
            pending={pool.submit(checked,r,p,True):i for i,r in enumerate(p['rows']) if i}
            for future in as_completed(pending):
                i=pending[future];got[i]=future.result();data['shards']=[got[j] for j in sorted(got)];checkpoint(D/'result.json',data);print('NEW11952 SLICE',got[i]['id'],got[i]['wall_seconds'],flush=True)
        data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(D/'result.json',data)
    except Exception as exc:
        data.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'result.json',data);raise

def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_hash']!=digest(p) or len(d['shards'])!=8:raise ArithmeticError('full new11952 slice cohort required')
    for r,s in zip(p['rows'],d['shards']):
        if checked(r,p,False)!=s:raise ArithmeticError('new11952 score/population replay differs')
    out=D/'replay.json'
    result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'retained_scores':32768,'primitive_population':p['total_primitive_population'],'disjoint_previous_complete_square':True}
    if out.exists():
        if cert.read(out)!=result:raise ArithmeticError('immutable replay differs')
    else:checkpoint(out,result)
    print('REPLAYED32768 NEW11952 RETAINED SCORES',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()

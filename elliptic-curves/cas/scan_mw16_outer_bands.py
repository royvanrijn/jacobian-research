#!/usr/bin/env python3
"""Target-free MW16 disjoint parameter bands, with bounded stratified coverage."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
from math import gcd, log
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
import select_prospective_mw16_wide as old
import verify_periodic_nagao_scanner as counts
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits
ROOT=Path(__file__).resolve().parents[2]; CAS=ROOT/'elliptic-curves/cas'
LOCAL=ROOT/'artifacts/local/elliptic-curves'; ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=LOCAL/'mw16-outer-bands-v1'; TABLES=LOCAL/'prospective-mw16-h1024-v1'
CPP=CAS/'newfamily/scan_rational_nagao_annulus.cpp'; BINARY=D/'scanner'
REFERENCE=LOCAL/'periodic-nagao-scanner-benchmark-v1/scanner-v2'
BANDS=[(4096,16384,16),(16384,65536,256),(65536,262144,4096)]
SALT='target-free-mw16-disjoint-outer-bands-v1'
PRIMES=[p for p in range(5,4094) if _is_prime(p)]

def sources():
    paths=[Path(__file__).resolve(),CPP,spec.ATLAS,Path(spec.__file__),Path(old.__file__),Path(cert.__file__),Path(counts.__file__),CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py',CAS/'mod2_reduction_independence.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def roster():
    parent=cert.read(TABLES/'protocol.json'); rows=[]
    for band,(inner,outer,modulus) in enumerate(BANDS):
        for i,family in enumerate(parent['families']):
            for j,sign in enumerate((-1,1)):
                h=sha256(f'{SALT}|{band}|{family}|{sign}'.encode()).hexdigest()
                shard=2*(int(h,16)%(modulus//2))+(i+j)%2
                table=TABLES/family/f'tables-{sign}.txt'
                count=counts.population(outer,outer,shard,modulus)-counts.population(inner,inner,shard,modulus)
                rows.append({'id':f'b{band+1}-{family}-'+('neg' if sign<0 else 'pos'),'band':band+1,'family':family,'sign':sign,'inner':inner,'outer':outer,'shard':shard,'shards':modulus,'sha256_choice':h,'primitive_population':count,'table':{'path':str(table.relative_to(ROOT)),'sha256':cert.hashed(table)}})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve protocol')
    if BINARY.exists():raise FileExistsError('preserve existing binary')
    D.mkdir(parents=True,exist_ok=True)
    command=['g++','-O3','-std=c++17',str(CPP),'-o',str(BINARY)]
    c=capture(command,limits=Limits(60,1073741824),log_path=D/'build.log',separate_stderr=True,check=False)
    checkpoint(D/'build.json',{'command':command,'supervision':c.supervision,'stdout':c.stdout,'stderr':c.stderr})
    if c.supervision['outcome']!='completed' or c.supervision['returncode']!=0:raise ArithmeticError('build failed')
    parent=cert.read(old.DIRECTORY/'protocol.json')
    if parent['height']!=4096 or len(parent['families'])!=5:raise ArithmeticError('compact population differs')
    gates=[LOCAL/'mw16-top25-pari-followup-v1/ledger.json',LOCAL/'mw16-extended-prime-benchmark-v1/replay.json']
    if any(cert.read(p)['status']!='PASS' for p in gates):raise ArithmeticError('existing target-free point/trace gates required')
    rows=roster()
    if len(rows)!=30 or sum(r['shard']%2 for r in rows)!=15:raise ArithmeticError('balanced roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-outer-bands.v1','sources':sources(),'binary_sha256':cert.hashed(BINARY),'reference_sha256':cert.hashed(REFERENCE),'build_sha256':cert.hashed(D/'build.json'),'prior_parameter_protocol_sha256':cert.hashed(old.DIRECTORY/'protocol.json'),'gate_hashes':{str(p.relative_to(ROOT)):cert.hashed(p) for p in gates},'rows':rows,'keep_per_slice':512,'maximum_workers':2,'seconds_per_slice':60,'first_slice_per_band_cost_gate_seconds':20,'outer_seconds':900,'total_primitive_population':sum(r['primitive_population'] for r in rows),'gate':'Five established compact MW16 models carry sixteen exact generic sections and their target-free prospective fibres have yielded independent exceptional points with the current PARI coordinate backend. Broad coverage currently ends at4096. This fixed experiment measures incidence in new, disjoint parameter-height bands with equal budgets across all five families; it does not assume higher height improves rank or visibility.','population':'4096 < H <=16384,16384 < H <=65536,65536 < H <=262144, H=max(abs(n),d), primitive nonzero n/d. One denominator slice per band/family/sign modulo16,256,4096 respectively, SHA256 chosen with balanced parity. The inner square is excluded before heap admission. Not complete outer-square coverage or uniform sampling.','selection':'Retain512 per signed slice by the unchanged562-prime quantized score, good count, denominator and absolute numerator; preserve all15360 survivors for subsequent selection. The duplicated p5 held field is unused.','target_free_boundary':'No catalogue equations, known-record parameters, points, ranks, j-invariants or jump labels are read during selection or prospective point execution. Post-terminal novelty comparison only.','validation':'Before broad scanning, compare complete reference populations and filtered annulus ordering on ten signed MW16 frames, both ordinary and actual large-modulus frames. Independently recalculate every returned score using canonical projective trace tables. Independently count primitive addresses by inclusion-exclusion. Replay retained raw transcripts.','future_scope':'After scanner replay, a separately frozen extension may score15360 survivors through65521 and use disjoint primes65537..131071 for validation on at most60 finalists, four per family per band. Freeze prospective equation deduplication and maps before a fixed generic16-only point batch using43 charts at125000/10seconds. No public-record filtering; repeats of own measured equations may be reported after completion.','boundaries':'Candidate incidence and conditional point visibility only. No rank inferred from score, exact-rank or absence theorem, universal novelty, or complete parameter coverage.'})
    print('FROZEN30 OUTER MW16 SLICES',sum(r['primitive_population'] for r in rows),flush=True)

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['binary_sha256']!=cert.hashed(BINARY) or p['reference_sha256']!=cert.hashed(REFERENCE) or p['build_sha256']!=cert.hashed(D/'build.json'):raise ArithmeticError('source/build changed')
    if p['rows']!=roster():raise ArithmeticError('population changed')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['gate_hashes'].items()):raise ArithmeticError('gate changed')
    return p

def exact_scores(rows,family):
    f=next(r for r in cert.read(spec.ATLAS)['families'] if r['fibration_id']==family)
    model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    scores=[0]*len(rows);goods=[0]*len(rows); hashes={}
    for p in PRIMES:
        path=TABLES/family/'trace-tables'/f'{p}.json';d=cert.read(path)
        if d['input']!={'family':family,'model_hash':digest(model),'prime':p}:raise ArithmeticError('canonical model binding differs')
        hashes[str(path.relative_to(ROOT))]=cert.hashed(path)
        for i,r in enumerate(rows):
            t=r['numerator']*pow(r['denominator'],-1,p)%p if r['denominator']%p else p
            if d['good'][t]:
                ap=d['traces'][t];scores[i]+=round((2-ap)/(p+1-ap)*log(p)*10**12);goods[i]+=1
    if any((r['score_units'],r['good_primes'])!=(s,g) for r,s,g in zip(rows,scores,goods)):raise ArithmeticError('canonical retained score differs')
    return hashes

def invocation(row,N,M,inner,keep,shard,modulus,path,create,reference=False):
    command=[str(REFERENCE if reference else BINARY),str(ROOT/row['table']['path']),str(N),str(M),str(keep),str(shard),str(modulus)]
    if not reference:command.append(str(inner))
    if create:
        if path.exists():raise FileExistsError('preserve raw invocation')
        c=capture(command,limits=Limits(60,1073741824),log_path=path.with_suffix('.log'),separate_stderr=True,check=False)
        checkpoint(path,{'command':command,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=command or raw['supervision']['command']!=command or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('immutable scanner failed/censored')
    rows,summary=old.parse(raw['stdout'],row['sign'])
    count=counts.population(N,M,shard,modulus)-(0 if reference else counts.population(min(N,inner),min(M,inner),shard,modulus))
    if summary!=[N,M,keep,count,count,min(keep,count)] or raw['stdout'].splitlines().count(f'R {shard} {modulus}')!=1:raise ArithmeticError('population frame differs')
    if not reference and raw['stdout'].splitlines().count(f'I {inner}')!=1:raise ArithmeticError('inner bound differs')
    if len({r['parameter'] for r in rows})!=len(rows) or any(gcd(abs(r['numerator']),r['denominator'])!=1 or r['numerator']*row['sign']<=0 or not 1<=abs(r['numerator'])<=N or not 1<=r['denominator']<=M or (r['denominator']-1)%modulus!=shard or (not reference and max(abs(r['numerator']),r['denominator'])<=inner) for r in rows):raise ArithmeticError('address outside frozen population')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('order differs')
    return rows,raw

def benchmark(check=False):
    p=protocol();records=[]
    for row in p['rows']:
        # Each actual slice frame, with at least two denominators, plus one square per signed family.
        cases=[(31,2*row['shards'],row['shards'],row['shard'],row['shards'])]
        if row['band']==1:cases.append((31,29,7,0,1))
        for ci,(N,M,inner,shard,modulus) in enumerate(cases):
            folder=D/'benchmark'/row['id']/str(ci)
            full,raw=invocation(row,N,M,0,1000,shard,modulus,folder/'reference.json',not check,True)
            expected={(row['sign']*n,d) for n in range(1,N+1) for d in range(1+shard,M+1,modulus) if gcd(n,d)==1}
            if {(r['numerator'],r['denominator']) for r in full}!=expected:raise ArithmeticError('reference full small population differs')
            filtered=[r for r in full if max(abs(r['numerator']),r['denominator'])>inner]
            allrows=[]
            for keep in (1000,7):
                got,_=invocation(row,N,M,inner,keep,shard,modulus,folder/f'annulus-{keep}.json',not check)
                if got!=filtered[:keep]:raise ArithmeticError('annular complete population/topK differs')
                allrows+=got
            hashes=exact_scores(full+allrows,row['family'])
            records.append({'id':row['id'],'case':[N,M,inner,shard,modulus],'checked_population':len(full),'canonical_hash':digest(hashes),'raw_hashes':{q.name:cert.hashed(q) for q in sorted(folder.glob('*.json'))}})
    result={'status':'PASS','protocol_hash':digest(p),'cases':len(records),'records':records}
    if check:
        if cert.read(D/'benchmark.json')!=result:raise ArithmeticError('benchmark replay differs')
    else:checkpoint(D/'benchmark.json',result)
    print('MW16 ANNULUS BENCHMARK PASS',len(records),'FRAMES',flush=True)

def checked(row,p,create):
    path=D/row['id']/'raw.json'
    rows,raw=invocation(row,row['outer'],row['outer'],row['inner'],p['keep_per_slice'],row['shard'],row['shards'],path,create)
    return {**row,'rows':rows,'raw_sha256':cert.hashed(path),'canonical_hash':digest(exact_scores(rows,row['family'])),'wall_seconds':raw['supervision']['wall_seconds']}

def run():
    p=protocol(); benchmark(True)
    if (D/'result.json').exists():raise FileExistsError('preserve fixed campaign')
    data={'status':'RUNNING','protocol_hash':digest(p),'shards':[]};checkpoint(D/'result.json',data);got={}
    for band in (1,2,3):
        i=next(i for i,r in enumerate(p['rows']) if r['band']==band)
        got[i]=checked(p['rows'][i],p,True);data['shards']=[got[j] for j in sorted(got)];checkpoint(D/'result.json',data)
        if got[i]['wall_seconds']>p['first_slice_per_band_cost_gate_seconds']:raise ArithmeticError('band cost gate failed; remaining slices not launched')
        print('MW16 OUTER BAND COST PASS',band,got[i]['wall_seconds'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(checked,r,p,True):i for i,r in enumerate(p['rows']) if i not in got}
        for f in as_completed(pending):
            i=pending[f];got[i]=f.result();data['shards']=[got[j] for j in sorted(got)];checkpoint(D/'result.json',data);print('MW16 OUTER SLICE',got[i]['id'],flush=True)
    data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(D/'result.json',data)

def replay():
    p=protocol();data=cert.read(D/'result.json')
    if data['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or data['protocol_hash']!=digest(p):raise ArithmeticError('incomplete population')
    for r,s in zip(p['rows'],data['shards']):
        if checked(r,p,False)!=s:raise ArithmeticError('slice replay differs')
    if len(data['shards'])!=30:raise ArithmeticError('missing slice')
    checkpoint(D/'replay.json',{'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'retained_scores':15360,'primitive_population':p['total_primitive_population']})
    print('REPLAYED15360 MW16 OUTER SCORES',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=parser.parse_args()
    benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()

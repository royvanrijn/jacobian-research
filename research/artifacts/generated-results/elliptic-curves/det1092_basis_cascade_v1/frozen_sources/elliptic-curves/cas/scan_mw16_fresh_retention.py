#!/usr/bin/env python3
"""Fresh disjoint outer MW16 slices with eightfold short-score retention."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import scan_mw16_outer_bands as core
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
ROOT=core.ROOT;CAS=core.CAS;ART=core.ART;D=core.LOCAL/'mw16-fresh-outer-retention-v1';SALT='mw16-fresh-annulus-retention4096-v1'
PRIOR=core.D/'protocol.json';RESULT=ART/'outer60_mw16_experiment_v1.json';PORTABLE=ART/'outer60_mw16_point_portable_replay_v2.json'
def sources():
    paths=[Path(__file__).resolve(),Path(core.__file__),core.BINARY,core.REFERENCE,PRIOR,RESULT,PORTABLE]
    return {**core.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}
def roster():
    old=cert.read(PRIOR);rows=[]
    for r in old['rows']:
        if r['band']!=1:continue
        choices=[a for a in range(r['shard']%2,16,2) if a!=r['shard']];h=sha256(f"{SALT}|{r['family']}|{r['sign']}".encode()).hexdigest();shard=choices[int(h,16)%len(choices)]
        rows.append({'id':r['family']+('-neg' if r['sign']<0 else '-pos'),'family':r['family'],'band':1,'sign':r['sign'],'inner':4096,'outer':16384,'shard':shard,'shards':16,'excluded_previous_shard':r['shard'],'sha256_choice':h,'table':r['table'],'primitive_population':core.counts.population(16384,16384,shard,16)-core.counts.population(4096,4096,shard,16)})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fresh retention protocol')
    if cert.read(RESULT)['status']!='PASS' or cert.read(PORTABLE)['status']!='PASS' or cert.read(PORTABLE)['logical_stages']!=182:raise ArithmeticError('complete outer point/proof cohort required')
    rows=roster()
    if len(rows)!=10 or sum(r['shard']%2 for r in rows)!=5 or any(r['shard']==r['excluded_previous_shard'] for r in rows):raise ArithmeticError('balanced fresh slice roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-fresh-retention.v1','sources':sources(),'rows':rows,'keep_per_slice':4096,'maximum_workers':2,'seconds_per_slice':60,'first_slice_cost_gate_seconds':20,'outer_seconds':600,'total_primitive_population':sum(r['primitive_population'] for r in rows),'gate':'The completed outer60 trial supplies two certified23-point fibres in4096<H<=16384. One is504th of512 in its signed short-score retention and1016th of1024 in the merged family pool; a smaller short prefix would miss it. This finite observation motivates deeper retention on fresh parameter addresses, not a rank-density or calibrated selector claim. The higher two bands remain bounded low-yield observations, not rank exclusions.','population':'Exactly ten new signed denominator residues modulo16 in4096<max(abs(n),d)<=16384, all five families equally represented. Each SHA256 choice is among seven same-parity residues excluding that frame of the previous outer-band experiment. Thus every scanned primitive address is disjoint from all its previous three bands and from the compact4096 square. Parameter novelty alone is not equation novelty.','retention':'Retain4096 per signed slice by unchanged562-prime quantized S1, good count, denominator and absolute numerator. Preserve all40960 candidates for separately frozen longer-prime scoring. No rank, public target, catalogue model or exceptional point enters shard selection or the scanner.','validation':'Reuse the pinned annular binary after its completed40-frame exact tests. Additionally compare complete reference and annular outputs, including top7 order, in each of the ten new modulo16 frames. Recompute all40960 retained scores from canonical projective trace tables and independently count primitive addresses. One full-slice runtime gate is reused once; no retries or replacement slices.','future_scope':'After complete scanner replay, a separately frozen40960-curve extension may select twelve per family,60 total, by combined S1 through65521, with disjoint65537..131071 validation. A later fixed generic16-only cohort may use43 charts per fibre at125000/10seconds, all maps before points. Known-record data stay outside selection and point execution. No original point search runs under this scanner protocol.','boundaries':'Finite disjoint parameter search and exact retained scores. No complete annulus coverage, global score optimum, rank prediction, point absence, upper rank bound or universal novelty.'})
    print('FROZEN10 FRESH MW16 SLICES',sum(r['primitive_population'] for r in rows),'ADDRESSES40960 RETENTION',flush=True)
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('fresh frozen inputs differ')
    for row in p['rows']:
        if cert.hashed(ROOT/row['table']['path'])!=row['table']['sha256']:raise ArithmeticError('signed table changed')
    return p

def benchmark(check=False):
    p=protocol();records=[]
    for row in p['rows']:
        folder=D/'benchmark'/row['id'];full,_=core.invocation(row,31,32,0,1000,row['shard'],16,folder/'reference.json',not check,True)
        expected={(row['sign']*n,d) for n in range(1,32) for d in range(1+row['shard'],33,16) if core.gcd(n,d)==1}
        if {(r['numerator'],r['denominator']) for r in full}!=expected:raise ArithmeticError('complete small reference population differs')
        filtered=[r for r in full if max(abs(r['numerator']),r['denominator'])>16]
        for keep in (1000,7):
            rows,_=core.invocation(row,31,32,16,keep,row['shard'],16,folder/f'annulus-{keep}.json',not check)
            if rows!=filtered[:keep]:raise ArithmeticError('fresh annulus/top7 differs')
        tables=core.exact_scores(full,row['family']);records.append({'id':row['id'],'population':len(full),'table_hash':digest(tables),'raw_hashes':{q.name:cert.hashed(q) for q in sorted(folder.glob('*.json'))}})
    result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'records':records}
    if check:
        if cert.read(D/'benchmark.json')!=result:raise ArithmeticError('fresh benchmark replay differs')
    else:
        if (D/'benchmark.json').exists():raise FileExistsError('preserve fresh benchmark')
        checkpoint(D/'benchmark.json',result)
    print('ALL10 FRESH FRAME TESTS PASS',flush=True)
def checked(row,p,create):
    path=D/row['id']/'raw.json';rows,raw=core.invocation(row,16384,16384,4096,4096,row['shard'],16,path,create)
    if len(rows)!=4096:raise ArithmeticError('full4096 retention required')
    hashes=core.exact_scores(rows,row['family'])
    return {**row,'rows':rows,'raw_sha256':cert.hashed(path),'canonical_hash':digest(hashes),'wall_seconds':raw['supervision']['wall_seconds']}
def run():
    p=protocol();benchmark(True)
    if (D/'result.json').exists():raise FileExistsError('preserve fresh parameter run')
    data={'status':'RUNNING','protocol_hash':digest(p),'shards':[]};checkpoint(D/'result.json',data);got={}
    got[0]=checked(p['rows'][0],p,True);data['shards']=[got[0]];checkpoint(D/'result.json',data)
    if got[0]['wall_seconds']>p['first_slice_cost_gate_seconds']:raise ArithmeticError('first fresh slice runtime gate failed')
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(checked,r,p,True):i for i,r in enumerate(p['rows']) if i}
        for future in as_completed(pending):
            i=pending[future];got[i]=future.result();data['shards']=[got[j] for j in sorted(got)];checkpoint(D/'result.json',data);print('FRESH MW16 SLICE',got[i]['id'],got[i]['wall_seconds'],flush=True)
    data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(D/'result.json',data)
def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_hash']!=digest(p) or len(d['shards'])!=10:raise ArithmeticError('full fresh slice cohort required')
    for r,s in zip(p['rows'],d['shards']):
        if checked(r,p,False)!=s:raise ArithmeticError('fresh score/population replay differs')
    checkpoint(D/'replay.json',{'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'retained_scores':40960,'primitive_population':p['total_primitive_population'],'disjoint_previous_band_slices':True})
    print('REPLAYED40960 FRESH RETAINED MW16 SCORES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()

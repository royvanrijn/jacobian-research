#!/usr/bin/env python3
"""A fixed twelve-slice H32768 prospective population using exact cached scores."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import benchmark_periodic_nagao_scanner as engine
import verify_periodic_nagao_scanner as strict
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=engine.ROOT;CAS=engine.CAS;ART=engine.ART;LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'higher32768-r17-stratified-v1';BINARY=engine.D/'scanner-v2';SALT='elliptic-stratified-height-32768-v1'
def sources():
    return {**engine.sources(),str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve()),str(Path(strict.__file__).resolve().relative_to(ROOT)):cert.hashed(Path(strict.__file__).resolve()),str(BINARY.relative_to(ROOT)):cert.hashed(BINARY)}
def roster():
    parent=cert.read(engine.old.DIRECTORY/'protocol.json');rows=[]
    for i,family in enumerate(parent['families']):
        for j,sign in enumerate((-1,1)):
            h=sha256(f'{SALT}|{family}|{sign}'.encode()).hexdigest();shard=2*(int(h,16)%32)+(i+j)%2
            rows.append({'id':family+('-negative' if sign<0 else '-positive'),'family':family,'sign':sign,'sha256_choice':h,'shard':shard,'shards':64,'table':parent['trace_tables'][family][str(sign)],'primitive_population':strict.population(32768,32768,shard,64)})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve stratified population')
    gates=[ART/'public_compact_parameter_heights_v1.json',ART/'periodic_nagao_scanner_benchmark_v1.json',ART/'periodic_nagao_scanner_strict_v1.json']
    if any(cert.read(p)['status']!='PASS' for p in gates):raise ArithmeticError('exact population/engine gates required')
    rows=roster()
    if len(rows)!=12 or sum(r['shard']%2 for r in rows)!=6:raise ArithmeticError('fixed balanced twelve-slice roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher32768-r17-stratified.v1','sources':sources(),'gate_hashes':{str(p.relative_to(ROOT)):cert.hashed(p) for p in gates},'rows':rows,'numerator_bound':32768,'denominator_bound':32768,'keep_per_slice':512,'prime_bound':4093,'seconds_per_slice':60,'rss_bytes_per_worker':1073741824,'maximum_workers':2,'outer_seconds':900,'total_primitive_population':sum(r['primitive_population'] for r in rows),'selection_salt':SALT,'gate':'Every mapped public rank-at-least28 calibration fibre lies outside the earlier16384 box. That selected historical sample supplies an incidence diagnostic, not a probability law. The periodic scanner passes exact output, population-count and all-twelve-frame small exhaustive checks, with a measured4.7-fold speedup on one higher-height slice. Test a new population at32768 with roughly the same total number of primitive addresses as the earlier full4096 experiment.','population':'For each of the six atlas families and each sign, keep one denominator residue modulo64. SHA256 of the literal salt, family and sign selects one of32 residues of a preassigned parity; atlas index plus sign index alternates parity, giving six odd and six even denominator slices. This is a deterministic stratified population, not uniform sampling or a complete32768 box. Public parameters, equations, points and ranks never enter the slice selector or scanner. No exclusion or refill of previously tested addresses.','score':'Retain512 by all562 cached prime scores per signed slice. Preserve within-slice tie order by denominator then absolute numerator. Merge both signs only to assign1024 per-family indices; do not discard any of the6144 returned addresses.','future_scope':'After complete short-score replay, a separately frozen extension may evaluate all6144 saved addresses through65521 with selection through32749 and disjoint validation afterwards, and at most four generic17-only point attempts per family. No point search is authorized by this scanner protocol itself.','checkpoints':'Twelve immutable raw scanner calls, aggregate per-slice checkpoints, no automatic retries after failure/censoring. At most two selector workers alongside the already fixed one-curve MW16 point follow-up.','boundaries':'Exact returned scores and primitive population counts; full ranking trusts the pinned worker. No rank prediction theorem, global score optimum, density, point absence or rank upper bound.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen stratified sources or roster changed')
    for r in p['rows']:
        if cert.hashed(ROOT/r['table']['path'])!=r['table']['sha256']:raise ArithmeticError('signed table changed')
    return p

def checked(row,p,create):
    folder=D/row['id'];path=folder/'raw.json';cmd=[str(BINARY),str(ROOT/row['table']['path']),'32768','32768','512',str(row['shard']),'64']
    if create:
        if path.exists():raise FileExistsError('preserve signed slice execution')
        c=capture(cmd,limits=Limits(p['seconds_per_slice'],p['rss_bytes_per_worker']),log_path=folder/'scanner.log',separate_stderr=True,check=False);checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('immutable scanner failed/censored or differs')
    rows,summary=engine.old.parse(raw['stdout'],row['sign']);count=row['primitive_population']
    if summary!=[32768,32768,512,count,count,512] or raw['stdout'].splitlines().count(f"R {row['shard']} 64")!=1:raise ArithmeticError('population framing differs')
    if len(rows)!=512 or len({r['parameter'] for r in rows})!=512 or any(engine.gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=32768 or not 1<=r['denominator']<=32768 or r['numerator']*row['sign']<=0 or (r['denominator']-1)%64!=row['shard'] for r in rows):raise ArithmeticError('primitive signed slice roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('ordered slice differs')
    hashes=engine.exact_scores(rows,{'family':row['family']})
    return {**row,'status':'PASS','rows':rows,'summary':summary,'raw_sha256':cert.hashed(path),'canonical_table_hashes':hashes,'wall_seconds':raw['supervision']['wall_seconds']}
def merge(shards):
    rows=[]
    for family in sorted({r['family'] for r in shards}):
        pool=[r for s in shards if s['family']==family for r in s['rows']];pool.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
        if len(pool)!=1024 or len({r['parameter'] for r in pool})!=1024:raise ArithmeticError('per-family fixed1024 pool differs')
        rows += [{**r,'family':family,'retained_index':i} for i,r in enumerate(pool)]
    return rows

def run():
    p=protocol();path=D/'result.json'
    if path.exists():raise FileExistsError('preserve stratified batch')
    data={'status':'RUNNING','protocol_hash':digest(p),'shards':[]};checkpoint(path,data);got={}
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(checked,r,p,True):i for i,r in enumerate(p['rows'])}
        for f in as_completed(pending):
            i=pending[f];got[i]=f.result();data['shards']=[got[j] for j in sorted(got)];checkpoint(path,data);print('HIGHER R17 SLICE',got[i]['id'],got[i]['primitive_population'],got[i]['wall_seconds'],flush=True)
    data['rows']=merge(data['shards']);data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(path,data)
def replay():
    p=protocol();d=cert.read(D/'result.json');shards=[checked(r,p,False) for r in p['rows']]
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_hash']!=digest(p) or d['shards']!=shards or d['rows']!=merge(shards):raise ArithmeticError('exact stratified population replay differs')
    print('REPLAYED12 SLICES,6144 EXACT562-PRIME SCORES AND',p['total_primitive_population'],'PRIMITIVE ADDRESSES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();globals()[a.stage]()

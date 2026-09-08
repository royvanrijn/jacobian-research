#!/usr/bin/env python3
"""Twelve fixed H131072 denominator slices at the prior primitive-work scale."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import benchmark_periodic_nagao_scanner as engine
import verify_periodic_nagao_scanner as strict
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=engine.ROOT;CAS=engine.CAS;ART=engine.ART;LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'outer131072-r17-stratified-v1';BINARY=engine.D/'scanner-v2';REFERENCE=engine.D/'scanner-v1';SALT='elliptic-stratified-height-131072-v1'
def sources():
    return {**engine.sources(),str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve()),str(Path(strict.__file__).resolve().relative_to(ROOT)):cert.hashed(Path(strict.__file__).resolve()),str(BINARY.relative_to(ROOT)):cert.hashed(BINARY),str(REFERENCE.relative_to(ROOT)):cert.hashed(REFERENCE)}
def roster():
    parent=cert.read(engine.old.DIRECTORY/'protocol.json');rows=[]
    for i,family in enumerate(parent['families']):
        for j,sign in enumerate((-1,1)):
            h=sha256(f'{SALT}|{family}|{sign}'.encode()).hexdigest();shard=2*(int(h,16)%512)+(i+j)%2
            rows.append({'id':family+('-negative' if sign<0 else '-positive'),'family':family,'sign':sign,'sha256_choice':h,'shard':shard,'shards':1024,'table':parent['trace_tables'][family][str(sign)],'primitive_population':strict.population(131072,131072,shard,1024)})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve outer131072 population')
    gates=[ART/'public_compact_parameter_heights_v1.json',ART/'periodic_nagao_scanner_benchmark_v1.json',ART/'periodic_nagao_scanner_strict_v1.json',ART/'compact192_experiment_v1.json',ART/'compact192_specialized_followup_v1.json']
    if any(cert.read(path)['status']!='PASS' for path in gates):raise ArithmeticError('exact population/engine gates required')
    control=LOCAL/'native11952-height125-control-v1/125000/verification.json';c=cert.read(control);gates.append(control)
    if c['status']!='PASS' or c['rank_lower_bound']!=28 or c['completed_boxes']!=49:raise ArithmeticError('complete generic17-only high-parameter visibility control required')
    rows=roster()
    if len(rows)!=12 or sum(r['shard']%2 for r in rows)!=6:raise ArithmeticError('fixed balanced twelve-slice roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.outer131072-r17-stratified.v1','sources':sources(),'gate_hashes':{str(path.relative_to(ROOT)):cert.hashed(path) for path in gates},'rows':rows,'numerator_bound':131072,'denominator_bound':131072,'prior_inner_bound':32768,'keep_per_slice':512,'prime_bound':4093,'seconds_per_slice':120,'first_slice_cost_gate_seconds':15,'small_case':{'numerator_bound':31,'denominator_bound':2048,'keep':1000,'seconds_per_call':20},'rss_bytes_per_worker':1073741824,'maximum_workers':2,'outer_seconds':1200,'total_primitive_population':sum(r['primitive_population'] for r in rows),'selection_salt':SALT,'gate':'The completed compact192 cohort adds84 curves but no27 result, and its five own26-point follow-ups certify no stronger lower bound. The existing generic17-only native rank29 control certifies28 at point height125000; its independently recovered compact parameter height89074 lies outside the earlier32768 square and within131072. This motivates a wider parameter envelope with an already successful finite visibility control, not a rank-density law. Modulus1024 denominator slices keep roughly the former primitive-address work while changing the sampled population.','population':'One denominator residue modulo1024 per family and sign, SHA256-selected from512 residues of a preassigned parity using only the literal salt, family and sign. Six odd and six even slices. No public parameter, equation, point, rank or target neighbourhood enters the slice selector or scanner. This is deterministic stratification, not full131072 coverage or uniform sampling. Inner32768 addresses remain in the raw score roster for provenance; prospective outer-cohort selection must exclude them.','score':'Retain512 by the same562 cached prime scores per signed slice and unchanged denominator/absolute-numerator ties; merge both signs without dropping any of6144 saved addresses. No score alteration or validation-based tie breaking.','cost_and_replay':'Before any full slice, compare pinned old and periodic scanner binaries on the entire31-by2048 population in each of the12 selected modulo1024 slices, checking every returned score and exact primitive counts. Then execute the first full slice once as part of the cohort; require completion within15 seconds before the remaining11 slices. At most two full selector workers. Preserve all raw calls, failures and censoring; no retry or replacement shard.','future_scope':'Only after all twelve slices and6144 returned scores replay, a separate bounded protocol may extend selection through32749 and disjoint validation through65521. A later outer-only point roster may contain at most eight candidates per family,48 total, at the existing control-gated point budget; all model deduplication and map files must precede points. No point search or trace-extension campaign is launched by this scanner protocol.','boundaries':'Exact retained scores and independently counted primitive populations. Full ranking trusts the pinned worker plus its exact regression gates. No complete parameter coverage, inferred rank, density estimate, exclusion of higher-rank curves, rank upper bound, record or universal novelty.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen stratified sources or roster changed')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['gate_hashes'].items()):raise ArithmeticError('frozen population gates changed')
    for r in p['rows']:
        if cert.hashed(ROOT/r['table']['path'])!=r['table']['sha256']:raise ArithmeticError('signed table changed')
    return p

def checked(row,p,create):
    folder=D/row['id'];path=folder/'raw.json';cmd=[str(BINARY),str(ROOT/row['table']['path']),'131072','131072','512',str(row['shard']),'1024']
    if create:
        if path.exists():raise FileExistsError('preserve signed slice execution')
        c=capture(cmd,limits=Limits(p['seconds_per_slice'],p['rss_bytes_per_worker']),log_path=folder/'scanner.log',separate_stderr=True,check=False);checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('immutable scanner failed/censored or differs')
    rows,summary=engine.old.parse(raw['stdout'],row['sign']);count=row['primitive_population']
    if summary!=[131072,131072,512,count,count,512] or raw['stdout'].splitlines().count(f"R {row['shard']} 1024")!=1:raise ArithmeticError('population framing differs')
    if len(rows)!=512 or len({r['parameter'] for r in rows})!=512 or any(engine.gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=131072 or not 1<=r['denominator']<=131072 or r['numerator']*row['sign']<=0 or (r['denominator']-1)%1024!=row['shard'] for r in rows):raise ArithmeticError('primitive signed slice roster differs')
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

def benchmark(check=False):
    p=protocol();path=D/'small-benchmark.json';records=[]
    if not check and path.exists():raise FileExistsError('preserve stride1024 benchmark')
    for row in p['rows']:
        variants=[];case=p['small_case'];count=strict.population(31,2048,row['shard'],1024)
        for label,binary in [('reference',REFERENCE),('periodic',BINARY)]:
            rawpath=D/'small-benchmark'/row['id']/(label+'.json')
            cmd=[str(binary),str(ROOT/row['table']['path']),'31','2048','1000',str(row['shard']),'1024']
            if not check:
                if rawpath.exists():raise FileExistsError('preserve individual stride benchmark')
                result=capture(cmd,limits=Limits(case['seconds_per_call'],p['rss_bytes_per_worker']),log_path=rawpath.with_suffix('.log'),separate_stderr=True,check=False)
                checkpoint(rawpath,{'command':cmd,'stdout':result.stdout,'stderr':result.stderr,'supervision':result.supervision})
            raw=cert.read(rawpath);points,summary=engine.old.parse(raw['stdout'],row['sign'])
            if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or raw['stderr'] or summary!=[31,2048,1000,count,count,count]:raise ArithmeticError('complete small stride benchmark failed/censored')
            expected={(row['sign']*n,d) for d in range(1+row['shard'],2049,1024) for n in range(1,32) if engine.gcd(n,d)==1}
            if {(r['numerator'],r['denominator']) for r in points}!=expected or len(points)!=count:raise ArithmeticError('full small primitive roster differs')
            hashes=engine.exact_scores(points,{'family':row['family']});variants.append(points)
            records.append({'id':row['id'],'binary':label,'raw_sha256':cert.hashed(rawpath),'population':count,'canonical_table_hashes':hashes})
        if variants[0]!=variants[1]:raise ArithmeticError('old/periodic complete stride outputs differ')
    result={'status':'PASS','protocol_hash':digest(p),'records':records,'cases':24}
    if check:
        if cert.read(path)!=result:raise ArithmeticError('stride benchmark replay differs')
    else:checkpoint(path,result)
    print('ALL12 STRIDE1024 FRAMES,24 EXACT SMALL BINARY CHECKS PASS',flush=True)

def run():
    p=protocol();path=D/'result.json';benchmark(True)
    if path.exists():raise FileExistsError('preserve outer131072 batch')
    data={'status':'RUNNING','protocol_hash':digest(p),'shards':[]};checkpoint(path,data);got={}
    got[0]=checked(p['rows'][0],p,True);data['shards']=[got[0]];checkpoint(path,data)
    if got[0]['wall_seconds']>p['first_slice_cost_gate_seconds']:raise ArithmeticError('first slice cost gate failed; no remaining slices')
    print('OUTER131072 FIRST SLICE COST PASS',got[0]['wall_seconds'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(checked,r,p,True):i for i,r in enumerate(p['rows']) if i!=0}
        for future in as_completed(pending):
            i=pending[future];got[i]=future.result();data['shards']=[got[j] for j in sorted(got)];checkpoint(path,data);print('OUTER131072 SLICE',got[i]['id'],got[i]['primitive_population'],got[i]['wall_seconds'],flush=True)
    data['rows']=merge(data['shards']);data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(path,data)

def replay():
    p=protocol();benchmark(True);d=cert.read(D/'result.json');shards=[checked(r,p,False) for r in p['rows']]
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_hash']!=digest(p) or d['shards']!=shards or d['rows']!=merge(shards):raise ArithmeticError('exact stratified population replay differs')
    print('REPLAYED12 SLICES,6144 EXACT562-PRIME SCORES AND',p['total_primitive_population'],'PRIMITIVE ADDRESSES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()

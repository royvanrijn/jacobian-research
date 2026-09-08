#!/usr/bin/env python3
"""Complete short-score population on11952, retaining512 per signed denominator slice."""
import argparse,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import benchmark_periodic_nagao_scanner as engine
import verify_periodic_nagao_scanner as strict
import benchmark_retained_short_cache as short
import benchmark_retained_extended_cache as extended
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=engine.ROOT;CAS=engine.CAS;ART=engine.ART;D=ROOT/'artifacts/local/elliptic-curves/full11952-h131072-short-v1'
BINARY=engine.D/'scanner-v2';GATE=ART/'million_retained_cache_benchmark_v1.json';N=131072;MODULUS=1024;KEEP=512

def sources():
    paths=[Path(__file__).resolve(),Path(strict.__file__).resolve(),GATE,ROOT/'artifacts/local/elliptic-curves/million-retained-cache-benchmark-v1/check.supervisor.json',short.OUT,extended.OUT,short.BINARY,extended.BINARY,ART/'outer48_r17_results_v1.json']
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},str(BINARY.relative_to(ROOT)):cert.hashed(BINARY)}

def slice_counts():
    spf=list(range(N+1))
    for p in range(2,N+1):
        if spf[p]==p:
            for k in range(p*p,N+1,p):
                if spf[k]==k:spf[k]=p
    counts=[0]*MODULUS
    for d in range(1,N+1):
        z=d;primes=[]
        while z>1:
            p=spf[z];primes.append(p)
            while z%p==0:z//=p
        divisors=[(1,1)]
        for p in primes:divisors += [(v*p,-s) for v,s in list(divisors)]
        counts[(d-1)%MODULUS]+=sum(s*(N//v) for v,s in divisors)
    mu=[1]*(N+1);mu[0]=0
    for p in _primes_up_to(N):
        for k in range(p,N+1,p):mu[k]=-mu[k]
        for k in range(p*p,N+1,p*p):mu[k]=0
    whole=sum(mu[k]*(N//k)**2 for k in range(1,N+1))
    if sum(counts)!=whole or 2*whole!=20888422894:raise ArithmeticError('independent full primitive count differs')
    return counts

def roster():
    old=cert.read(engine.old.DIRECTORY/'protocol.json');counts=slice_counts();rows=[]
    for sign in [-1,1]:
        for shard in range(MODULUS):
            rows.append({'id':('negative' if sign<0 else 'positive')+f'-{shard:04}','family':'11952','sign':sign,'shard':shard,'shards':MODULUS,'table':old['trace_tables']['11952'][str(sign)],'primitive_population':counts[shard]})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve full11952 population')
    gate=cert.read(GATE);check=cert.read(ROOT/'artifacts/local/elliptic-curves/million-retained-cache-benchmark-v1/check.supervisor.json')
    if gate['status']!='PASS' or not gate['cost_gate_passed'] or check['outcome']!='completed' or check['returncode']!=0:raise ArithmeticError('full-size cached-score cost and exact agreement gates required')
    rows=roster();counts=[r['primitive_population'] for r in rows[:MODULUS]]
    for shard in [0,1,511,512,1022,1023]:
        if counts[shard]!=strict.population(N,N,shard,MODULUS):raise ArithmeticError('independent slice boundary regression differs')
    bands={}
    for label,module in [('short',short),('extended',extended)]:
        proof=cert.read(module.OUT)
        if proof['status']!='PASS' or cert.hashed(module.TABLE)!=proof['cache_sha256']:raise ArithmeticError('exact cached-score source changed')
        bands[label]={'binary':str(module.BINARY.relative_to(ROOT)),'binary_sha256':cert.hashed(module.BINARY),'cache':str(module.TABLE.relative_to(ROOT)),'cache_sha256':proof['cache_sha256']}
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.full11952-h131072-short.v1','sources':sources(),'recorded_root':str(ROOT),'rows':rows,'numerator_bound':N,'denominator_bound':N,'keep_per_slice':KEEP,'short_prime_bound':4093,'total_primitive_population':20888422894,'retained_rows':1048576,'maximum_workers':8,'checkpoint_block_size':16,'seconds_per_slice':120,'rss_bytes_per_worker':1073741824,'first_eight_block_gate_seconds':60,'scan_wall_seconds':7200,'cached_score_bands':bands,'gate':'The known native11952 rank29 control has a successful generic17-only point detector at125000 but shifts from short-score position435 to extended-score position2 against the same967 saved outer candidates. The earlier scan omitted its denominator residue. The completed selection-prime cache and exact full-size repeated-fixture benchmark now make it affordable to keep512 candidates from every signed denominator slice until the stronger score is available.11952 is chosen for this exact positive-control and retention diagnosis, not generic rank alone.','population':'All nonzero signed primitive n/d with abs(n),d<=131072 in11952, partitioned by sign and(d-1)mod1024 into2048 disjoint slices. The denominator-partition inclusion-exclusion count and independent Mobius square count agree on20888422894 addresses. No public parameter, equation, point or rank enters the enumeration, per-slice heap or prospective score ordering. Zero and infinity are outside this declared population.','retention':'Use the unchanged pinned periodic short scanner and signed562-prime tables, retaining exactly512 per slice under its original score and ties. Keep all1048576 survivors for cached short-score replay and extended S1 through32749; no global short-prefix truncation. This is complete short-population coverage, not complete extended-score or point coverage.','cost_gate':'Execute the first eight declared negative slices once as part of the population with eight workers; require their whole block to finish within60 seconds before any remaining slice. Later blocks contain at most16 calls. Checkpoint every completed slice and stop after the current block on any failure/censoring, with no retry, refill or changed slice.','future_scope':'Only after all2048 scanner invocations and all1048576 retained short scores replay may a separate fixed protocol apply the existing cached extended score, exclude rational isomorphs of all789 previously measured equations and earlier selected rows, and freeze at most64 new11952 point candidates. Every map must precede points; generic17 seeds,49 generic maximum classes, height125000, ten seconds per chart and provisional28 target remain the point limits. Validation and catalogue labels never enter selection. No point search is launched by this scanner protocol.','boundaries':'The full ranking within each slice trusts the pinned scanner and its exact regressions. Cached short-score replay is required before calling the retained values verified. Heuristic scores are not rank bounds, sampling/density laws, exact ranks, point absence or universal novelty.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen complete population changed')
    for row in [p['rows'][0],p['rows'][MODULUS]]:
        if cert.hashed(ROOT/row['table']['path'])!=row['table']['sha256']:raise ArithmeticError('signed table changed')
    return p

def checked(row,p,create=False):
    folder=D/row['id'];path=folder/'raw.json';cmd=[str(BINARY),str(ROOT/row['table']['path']),str(N),str(N),str(KEEP),str(row['shard']),str(MODULUS)]
    if create:
        if path.exists():raise FileExistsError('preserve one signed-slice scan')
        c=capture(cmd,limits=Limits(p['seconds_per_slice'],p['rss_bytes_per_worker']),log_path=folder/'scanner.log',separate_stderr=True,check=False)
        checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    expected=[str(Path(p['recorded_root'])/BINARY.relative_to(ROOT)),str(Path(p['recorded_root'])/row['table']['path']),str(N),str(N),str(KEEP),str(row['shard']),str(MODULUS)]
    if raw['command']!=expected or raw['supervision']['command']!=expected or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('slice failed/censored or command differs')
    rows,summary=engine.old.parse(raw['stdout'],row['sign']);count=row['primitive_population']
    if summary!=[N,N,KEEP,count,count,KEEP] or raw['stdout'].splitlines().count(f"R {row['shard']} {MODULUS}")!=1:raise ArithmeticError('complete slice framing differs')
    if len(rows)!=KEEP or len({r['parameter'] for r in rows})!=KEEP or any(engine.gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=N or not 1<=r['denominator']<=N or r['numerator']*row['sign']<=0 or (r['denominator']-1)%MODULUS!=row['shard'] for r in rows):raise ArithmeticError('primitive signed slice differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('per-slice score order differs')
    return {'id':row['id'],'status':'COMPLETE_PENDING_CACHED_SCORE_REPLAY','raw_sha256':cert.hashed(path),'primitive_population':count,'retained_rows':KEEP,'wall_seconds':raw['supervision']['wall_seconds']},rows

def run():
    p=protocol();path=D/'ledger.json'
    if path.exists():raise FileExistsError('preserve full11952 scan')
    ledger={'status':'RUNNING','protocol_hash':digest(p),'rows':[{'id':r['id'],'status':'PENDING'} for r in p['rows']]};checkpoint(path,ledger)
    def one(i):
        try:return checked(p['rows'][i],p,True)[0]
        except Exception as error:return {'id':p['rows'][i]['id'],'status':'FAILED_OR_CENSORED','reason':str(error)}
    try:
        with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
            blocks=[list(range(8))]+[list(range(i,min(i+p['checkpoint_block_size'],len(p['rows'])))) for i in range(8,len(p['rows']),p['checkpoint_block_size'])]
            for number,block in enumerate(blocks):
                start=time.monotonic();jobs={pool.submit(one,i):i for i in block}
                for future in as_completed(jobs):
                    i=jobs[future];ledger['rows'][i]=future.result();checkpoint(path,ledger)
                elapsed=time.monotonic()-start
                if any(ledger['rows'][i]['status']!='COMPLETE_PENDING_CACHED_SCORE_REPLAY' for i in block):raise ArithmeticError('failed/censored full-population block; no later block')
                if number==0:
                    ledger['first_eight_block_wall_seconds']=elapsed;checkpoint(path,ledger)
                    if elapsed>p['first_eight_block_gate_seconds']:raise ArithmeticError('first-eight concurrency cost gate failed')
                    print('FULL11952 FIRST EIGHT COST GATE PASS',elapsed,flush=True)
                if number%16==0:print('FULL11952 SLICES',max(block)+1,'OF',len(p['rows']),flush=True)
        ledger['status']='COMPLETE_FIXED_SCAN_PENDING_SCORE_REPLAY';checkpoint(path,ledger)
    except Exception as error:
        ledger['status']='FAILED_OR_CENSORED';ledger['reason']=str(error);checkpoint(path,ledger);raise
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','run']);v=a.parse_args();globals()[v.stage]()

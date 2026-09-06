#!/usr/bin/env python3
"""Fixed extended-prime S1 scoring of2048 retained skew-rectangle addresses."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import scan_skew_r17_boxes as parent
from math import log
import benchmark_r17_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=ROOT/'artifacts/local/elliptic-curves/skew-r17-extended-v1';OLD=ROOT/'artifacts/local/elliptic-curves/higher32768-r17-extended-v1'
def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(scoring.__file__).resolve(),spec.ATLAS)}}
def model_at(f,t):
    t=cert.F(t);q=t.denominator;model=[cert.F(0)]*3+[spec.polynomial(f['A_coefficients_low_to_high'],t)*q**8,spec.polynomial(f['B_coefficients_low_to_high'],t)*q**12]
    if any(x.denominator!=1 for x in model) or 4*model[3]**3+27*model[4]**2==0:raise ArithmeticError('integral nonsingular specialization required')
    return list(map(str,model))
CACHE_RESULTS=[OLD/'result.json',ROOT/'artifacts/local/elliptic-curves/higher32768-product-first-extended-v1/result.json']
def cache_bindings():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in CACHE_RESULTS}
def attach_cache(rows):
    cached={}
    for source in CACHE_RESULTS:
        old=cert.read(source)
        if old['status']!='COMPLETE_FROZEN_TRACE_EXTENSION':raise ArithmeticError('complete old trace caches required')
        for r in old['rows']:
            path=ROOT/r['raw_path'] if 'raw_path' in r else source.parent/r['family']/f"candidate-{r['retained_index']:04}"/'raw.json'
            value={'path':str(path.relative_to(ROOT)),'sha256':r['raw_sha256']}
            key=r['family'],r['parameter']
            if key in cached and cached[key][0]!=r['model']:raise ArithmeticError('overlap cache models differ')
            cached.setdefault(key,(r['model'],value))
    result=[]
    for row in rows:
        hit=cached.get((row['family'],row['parameter']));value=None
        if hit:
            if hit[0]!=row['model'] or cert.hashed(ROOT/hit[1]['path'])!=hit[1]['sha256']:raise ArithmeticError('cached model/trace differs')
            value=hit[1]
        result.append({**row,'cached_trace':value})
    return result

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher-height trace extension')
    p=parent.protocol();d=cert.read(parent.D/'result.json');replay=cert.read(parent.D/'replay.supervisor.json')
    if d['status']!='COMPLETE_FROZEN_SKEW_POPULATION' or replay['outcome']!='completed' or replay['returncode']!=0 or len(d['rows'])!=2048:raise ArithmeticError('complete2048-address short-score replay required')
    atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};rows=[{**r,'model':model_at(atlas[r['family']],r['parameter'])} for r in d['rows']];rows=attach_cache(rows);samples=[next(i for i,r in enumerate(rows) if r['family']==family and r['cached_trace'] is None) for family in sorted({r['family'] for r in rows})]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.skew-r17-prime-extension.v1','cache_result_hashes':cache_bindings(),'cached_rows':sum(r['cached_trace'] is not None for r in rows),'new_rows':sum(r['cached_trace'] is None for r in rows),'sources':sources(),'parent_protocol_sha256':cert.hashed(parent.D/'protocol.json'),'parent_result_sha256':cert.hashed(parent.D/'result.json'),'parent_replay_sha256':cert.hashed(parent.D/'replay.supervisor.json'),'rows':rows,'prime_roster':engine.PRIMES,'gp_sha256':cert.hashed(engine.GP),'benchmark_indices':samples,'direct_validation_indices':[0,1,2,len(engine.PRIMES)//4,len(engine.PRIMES)//2,3*len(engine.PRIMES)//4,len(engine.PRIMES)-2,len(engine.PRIMES)-1],'benchmark_seconds_per_curve':20,'seconds_per_curve':5,'rss_bytes':536870912,'maximum_workers':2,'outer_seconds':1200,'checkpoint_block_size':16,'gate':'The exact equal-area coefficient-bound gate selects two skew families. All four signed rectangular scans, complete small-shape scores and2048 retained562-prime scores replay. Extend this saved population, reusing every matching raw trace in both earlier higher-height pools. Before bulk work, the first uncached row per family must pass sixteen independent character sums and the unchanged two-second cost gate. No point result or catalogue label enters selection.','selection':'Among retained addresses with absolute numerator greater than32768, choose four per family by unchanged S1 score through32749, combined good count, denominator and signed numerator. This uniformly limits point attempts to the portion outside the old square. Higher validation primes through65521 never break ties. No result-dependent refill.','future_point_scope':'At most eight generic17-only point attempts after every2048 score and fixed selection replays, under a separate point protocol. No additional parameters or automatic point-height increase.','boundaries':'Finite selected denominator slices, with a short-prime retention bottleneck. No full32768 coverage, density estimate, rank classifier, upper bound or universal novelty. All failures/censored raw outputs retained; no automatic retry.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['parent_protocol_sha256']!=cert.hashed(parent.D/'protocol.json') or p['parent_result_sha256']!=cert.hashed(parent.D/'result.json') or p['gp_sha256']!=cert.hashed(engine.GP) or p['cache_result_hashes']!=cache_bindings():raise ArithmeticError('frozen higher trace inputs changed')
    return p

def call(row,p,path,seconds,create):
    program=engine.program(row['model'])
    if create:
        if path.exists():raise FileExistsError('preserve trace call')
        c=capture([str(engine.GP),'-q','-s','256000000'],input_text=program,limits=Limits(seconds,p['rss_bytes']),log_path=path.with_suffix('.log'),separate_stderr=True,check=False);checkpoint(path,{'program':program,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    r=cert.read(path)
    if r['program']!=program or r['stderr'] or r['supervision']['outcome']!='completed' or r['supervision']['returncode']!=0:raise ArithmeticError('trace call failed/censored or differs')
    return engine.parse(r['stdout'],row['model']),r

def benchmark(check=False):
    p=protocol();out=D/'benchmark.json';rows=[]
    if not check and out.exists():raise FileExistsError('preserve higher-height trace benchmark')
    for i in p['benchmark_indices']:
        row=p['rows'][i];path=D/'benchmark'/row['family']/'raw.json';(values,ms),raw=call(row,p,path,20,not check);checks=[]
        for j in p['direct_validation_indices']:
            q,t=values[j];v=engine.direct(row['model'],q)
            if t!=v:raise ArithmeticError('independent character sum differs')
            checks.append([q,v])
        if raw['supervision']['wall_seconds']>2:raise ArithmeticError('higher-height trace cost outside full-batch gate')
        rows.append({'index':i,'family':row['family'],'parameter':row['parameter'],'direct_checks':checks,'cpu_ms':ms,'wall_seconds':raw['supervision']['wall_seconds'],'raw_sha256':cert.hashed(path)})
    result={'status':'PASS','protocol_hash':digest(p),'rows':rows,'exact_character_sums':16}
    if check:
        if cert.read(out)!=result:raise ArithmeticError('benchmark replay differs')
    else:checkpoint(out,result)
    print('SKEW S1 TRACE BENCHMARK16 DIRECT CHECKS PASS',flush=True)
def checked_row(row,p,create):
    cache=row['cached_trace'];path=ROOT/cache['path'] if cache else D/row['family']/f"candidate-{row['retained_index']:04}"/'raw.json'
    if cache and cert.hashed(path)!=cache['sha256']:raise ArithmeticError('cached raw trace differs')
    (values,ms),raw=call(row,p,path,5,create and not cache);s=scoring.sums(values)
    return {**row,**s,'combined_selection_units':row['score_units']+s['extension_selection_units'],'combined_good':row['good_primes']+s['extension_good'],'cpu_ms':ms,'raw_path':str(path.relative_to(ROOT)),'raw_sha256':cert.hashed(path)}

def choose(rows):return [r['retained_index'] for r in sorted([r for r in rows if abs(r['numerator'])>32768],key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))[:4]]
def selection(rows):return {f:choose([r for r in rows if r['family']==f]) for f in sorted({r['family'] for r in rows})}
def run():
    p=protocol();out=D/'result.json';benchmark(True)
    if out.exists():raise FileExistsError('preserve full trace batch')
    d={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,d)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for start in range(0,2048,16):
            d['rows']+=list(pool.map(lambda r:checked_row(r,p,True),p['rows'][start:start+16]));checkpoint(out,d)
            if len(d['rows'])%256==0:print('SKEW S1 TRACE EXTENSION',len(d['rows']),'of2048',flush=True)
    d['selection']=selection(d['rows']);d['status']='COMPLETE_FROZEN_TRACE_EXTENSION';checkpoint(out,d)
def replay():
    p=protocol();d=cert.read(D/'result.json');base=cert.read(parent.D/'result.json');atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};expected=attach_cache([{**r,'model':model_at(atlas[r['family']],r['parameter'])} for r in base['rows']])
    if p['rows']!=expected or d['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or d['protocol_hash']!=digest(p) or len(d['rows'])!=2048:raise ArithmeticError('full2048 roster differs')
    for row,got in zip(p['rows'],d['rows']):
        if checked_row(row,p,False)!=got:raise ArithmeticError('exact trace/score replay differs')
    if d['selection']!=selection(d['rows']):raise ArithmeticError('fixed8 selector differs')
    print('REPLAYED2048 SKEW S1 TRACE ROSTERS AND FIXED8 SELECTION',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()

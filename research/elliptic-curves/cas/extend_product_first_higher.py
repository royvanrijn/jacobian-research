#!/usr/bin/env python3
"""Fixed extended-prime scoring of6144 saved higher-height R17 addresses."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import rescreen_higher_product_first as parent
from math import log
import benchmark_r17_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=ROOT/'artifacts/local/elliptic-curves/higher32768-product-first-extended-v1';OLD=ROOT/'artifacts/local/elliptic-curves/higher32768-r17-extended-v1'
def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(scoring.__file__).resolve(),spec.ATLAS)}}
def model_at(f,t):
    t=cert.F(t);q=t.denominator;model=[cert.F(0)]*3+[spec.polynomial(f['A_coefficients_low_to_high'],t)*q**8,spec.polynomial(f['B_coefficients_low_to_high'],t)*q**12]
    if any(x.denominator!=1 for x in model) or 4*model[3]**3+27*model[4]**2==0:raise ArithmeticError('integral nonsingular specialization required')
    return list(map(str,model))
def attach_cache(rows):
    old=cert.read(OLD/'result.json');cached={(r['family'],r['parameter']):r for r in old['rows']};result=[]
    if old['status']!='COMPLETE_FROZEN_TRACE_EXTENSION':raise ArithmeticError('complete old trace cache required')
    for row in rows:
        r=cached.get((row['family'],row['parameter']));value=None
        if r:
            if r['model']!=row['model']:raise ArithmeticError('cached model differs')
            path=OLD/r['family']/f"candidate-{r['retained_index']:04}"/'raw.json'
            if cert.hashed(path)!=r['raw_sha256']:raise ArithmeticError('cached trace changed')
            value={'path':str(path.relative_to(ROOT)),'sha256':r['raw_sha256']}
        result.append({**row,'cached_trace':value})
    return result

def product_sums(values):
    s=scoring.sums(values)
    for key in ('extension_selection_units','validation_units'):s[key]=0
    for q,t in values:
        if t is not None:s['extension_selection_units' if q<=32749 else 'validation_units']+=round(log((q+1-t)/q)*10**12)
    return s

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher-height trace extension')
    p=parent.protocol();d=cert.read(parent.D/'result.json');replay=cert.read(parent.D/'replay.supervisor.json')
    if d['status']!='COMPLETE_FIXED_PRODUCT_FIRST_POPULATION' or replay['outcome']!='completed' or replay['returncode']!=0 or len(d['rows'])!=6144:raise ArithmeticError('complete6144-address short-score replay required')
    atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};rows=[{**r,'model':model_at(atlas[r['family']],r['parameter'])} for r in d['rows']];rows=attach_cache(rows);samples=[next(i for i,r in enumerate(rows) if r['family']==family and r['cached_trace'] is None) for family in sorted(atlas)]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher-product-first-prime-extension.v1','old_trace_result_sha256':cert.hashed(OLD/'result.json'),'cached_rows':sum(r['cached_trace'] is not None for r in rows),'new_rows':sum(r['cached_trace'] is None for r in rows),'sources':sources(),'parent_protocol_sha256':cert.hashed(parent.D/'protocol.json'),'parent_result_sha256':cert.hashed(parent.D/'result.json'),'parent_replay_sha256':cert.hashed(parent.D/'replay.supervisor.json'),'rows':rows,'prime_roster':engine.PRIMES,'gp_sha256':cert.hashed(engine.GP),'benchmark_indices':samples,'direct_validation_indices':[0,1,2,len(engine.PRIMES)//4,len(engine.PRIMES)//2,3*len(engine.PRIMES)//4,len(engine.PRIMES)-2,len(engine.PRIMES)-1],'benchmark_seconds_per_curve':20,'seconds_per_curve':5,'rss_bytes':536870912,'maximum_workers':2,'outer_seconds':1200,'checkpoint_block_size':16,'gate':'The product-first rescreen uses exactly the earlier twelve slices and122368792 addresses. It retains4338 addresses absent from the old6144 S1 pool;1806 overlaps reuse their exact saved trace outputs. The paired second-stage null cannot test these discarded candidates. Extend only the4338 previously uncached addresses, after the first uncached address per family passes a bounded six-curve cost and48 character-sum benchmark. All6144 rows are scored, but final point selection uniformly excludes the1806 old-pool overlaps before ranking. Public points, measured ranks and catalogue labels do not enter selection.','selection':'Among addresses absent from the entire old S1-retained6144 pool, choose four per family by combined quantized product score through32749, then selection-band good count, denominator and signed numerator. Validation32771..65521 is excluded even from ties. Exclusion uses only old pool membership, not ranks or point outcomes. No result-dependent refill.','future_point_scope':'At most24 equally exposed generic17-only attempts using the complete43/49 exact generic maximum classes, after all6144 trace rosters and selections replay and a separate point protocol freezes height/time limits.','boundaries':'Finite selected denominator slices, with a short-prime retention bottleneck. No full32768 coverage, density estimate, rank classifier, upper bound or universal novelty. All failures/censored raw outputs retained; no automatic retry.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['parent_protocol_sha256']!=cert.hashed(parent.D/'protocol.json') or p['parent_result_sha256']!=cert.hashed(parent.D/'result.json') or p['gp_sha256']!=cert.hashed(engine.GP) or p['old_trace_result_sha256']!=cert.hashed(OLD/'result.json'):raise ArithmeticError('frozen higher trace inputs changed')
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
    result={'status':'PASS','protocol_hash':digest(p),'rows':rows,'exact_character_sums':48}
    if check:
        if cert.read(out)!=result:raise ArithmeticError('benchmark replay differs')
    else:checkpoint(out,result)
    print('PRODUCT-FIRST TRACE BENCHMARK48 DIRECT CHECKS PASS',flush=True)
def checked_row(row,p,create):
    cache=row['cached_trace'];path=ROOT/cache['path'] if cache else D/row['family']/f"candidate-{row['retained_index']:04}"/'raw.json'
    if cache and cert.hashed(path)!=cache['sha256']:raise ArithmeticError('cached raw trace differs')
    (values,ms),raw=call(row,p,path,5,create and not cache);s=product_sums(values)
    return {**row,**s,'combined_selection_units':row['score_units']+s['extension_selection_units'],'combined_good':row['good_primes']+s['extension_good'],'cpu_ms':ms,'raw_path':str(path.relative_to(ROOT)),'raw_sha256':cert.hashed(path)}

def choose(rows):return [r['retained_index'] for r in sorted([r for r in rows if r['cached_trace'] is None],key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))[:4]]
def selection(rows):return {f:choose([r for r in rows if r['family']==f]) for f in sorted({r['family'] for r in rows})}
def run():
    p=protocol();out=D/'result.json';benchmark(True)
    if out.exists():raise FileExistsError('preserve full trace batch')
    d={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,d)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for start in range(0,6144,16):
            d['rows']+=list(pool.map(lambda r:checked_row(r,p,True),p['rows'][start:start+16]));checkpoint(out,d)
            if len(d['rows'])%256==0:print('PRODUCT-FIRST TRACE EXTENSION',len(d['rows']),'of6144',flush=True)
    d['selection']=selection(d['rows']);d['status']='COMPLETE_FROZEN_TRACE_EXTENSION';checkpoint(out,d)
def replay():
    p=protocol();d=cert.read(D/'result.json');base=cert.read(parent.D/'result.json');atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};expected=attach_cache([{**r,'model':model_at(atlas[r['family']],r['parameter'])} for r in base['rows']])
    if p['rows']!=expected or d['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or d['protocol_hash']!=digest(p) or len(d['rows'])!=6144:raise ArithmeticError('full6144 roster differs')
    for row,got in zip(p['rows'],d['rows']):
        if checked_row(row,p,False)!=got:raise ArithmeticError('exact trace/score replay differs')
    if d['selection']!=selection(d['rows']):raise ArithmeticError('fixed24 selector differs')
    print('REPLAYED6144 PRODUCT-FIRST TRACE ROSTERS AND FIXED24 SELECTION',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()

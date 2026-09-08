#!/usr/bin/env python3
"""Cached32768-candidate extension, fresh4096 scalar trial and64 new fibres."""
import argparse,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import scan_11952_new_annulus as scan
import benchmark_retained_extended_cache as cache
import benchmark_r17_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
import extend_outer131072_r17 as models
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=scan.ROOT;CAS=scan.CAS;ART=scan.ART
D=ROOT/'artifacts/local/elliptic-curves/11952-new-annulus-scores-v1'
OUT=ART/'new_annulus_11952_selection_v1.json'
OLD=ART/'outer60_mw16_results_v1.json'
CAT=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1/database.json'
FRESH_PRIMES=[q for q in _primes_up_to(131071) if q>=65537]

def sources():
    paths=[Path(__file__).resolve(),Path(scan.__file__),Path(cache.__file__),cache.OUT,cache.TABLE,cache.BINARY,Path(engine.__file__),Path(scoring.__file__),Path(models.__file__),models.spec.ATLAS,OLD,CAT,scan.D/'protocol.json',scan.D/'result.json',scan.D/'replay.json',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    out=D/'bank-protocol.json'
    if out.exists():raise FileExistsError('preserve new11952 score bank')
    upstream=cert.read(scan.D/'controller/ledger.json');raw=cert.read(scan.D/'result.json');replay=cert.read(scan.D/'replay.json');proof=cert.read(cache.OUT)
    if upstream['status']!='PASS' or replay['status']!='PASS' or replay['result_sha256']!=cert.hashed(scan.D/'result.json') or proof['status']!='PASS' or cert.hashed(cache.TABLE)!=proof['cache_sha256'] or cert.hashed(cache.BINARY)!=proof['binary_sha256']:raise ArithmeticError('annular and exact cached-score gates required')
    f=next(r for r in cert.read(models.spec.ATLAS)['families'] if r['family']=='11952');rows=[]
    for shard in raw['shards']:
        for r in shard['rows']:
            i=len(rows);n=r['numerator'];d=r['denominator']
            if not scan.INNER<max(abs(n),d)<=scan.OUTER or cert.F(n,d)!=cert.F(r['parameter']):raise ArithmeticError('fresh primitive annulus differs')
            rows.append({**r,'id':f'11952-{i:07}','retained_index':i,'family':'11952','slice_id':shard['id'],'model':models.model_at(f,r['parameter'])})
    if len(rows)!=32768 or len({r['parameter'] for r in rows})!=32768:raise ArithmeticError('exact32768 distinct retained addresses required')
    checkpoint(out,{'schema':'elliptic-curves.11952-new-annulus-bank.v1','sources':sources(),'rows':rows,'cached_seconds':120,'cached_primes':2948,'rss_bytes':2147483648,'gate':'All eight disjoint131072<H<=524288 slices and32768 canonical short scores replay. The existing exact11952 projective extension cache is parameter-independent and previously agrees with scalar traces and complete encoded-table replay. Retain every short survivor for this cached extension, then freeze one bounded scalar cohort.','scope':'One cached4099..32749 score call on all32768 addresses, then4096 distinct equations selected by combined S1 through32749, with977 prior and593 catalogue equations used only as exclusions. No public rank, point or validation score enters ordering. A separate immutable scalar protocol binds its complete4096 roster before new scalar calls.'})
    print('NEW11952 FROZEN32768 CACHED SCORE BANK',flush=True)

def bank_protocol():
    p=cert.read(D/'bank-protocol.json')
    if p['sources']!=sources():raise ArithmeticError('immutable new11952 score bank differs')
    return p

def cached(check=False):
    p=bank_protocol();input_path=D/'candidates.txt';rawpath=D/'cache-raw.json';out=D/'cache-result.json'
    text='R17-CANDIDATES-V1 32768\n'+''.join(f"{r['numerator']} {r['denominator']}\n" for r in p['rows'])
    cmd=[str(cache.BINARY),str(cache.TABLE),str(input_path)]
    if not check:
        if any(q.exists() for q in (input_path,rawpath,out)):raise FileExistsError('preserve one cached call')
        input_path.write_text(text)
        c=capture(cmd,limits=Limits(p['cached_seconds'],p['rss_bytes']),log_path=D/'cache.log',separate_stderr=True,check=False)
        checkpoint(rawpath,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    if input_path.read_text()!=text:raise ArithmeticError('fixed candidate bank changed')
    raw=cert.read(rawpath);lines=raw['stdout'].splitlines();rows=[]
    if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or len(lines)!=32769 or lines[-1]!='S 32768 2948':raise ArithmeticError('cached invocation failed/censored or framing differs')
    for i,(r,line) in enumerate(zip(p['rows'],lines[:-1])):
        v=line.split()
        if len(v)!=4 or v[:2]!=['R',str(i)]:raise ArithmeticError('cached row address differs')
        s,g=map(int,v[2:])
        if not 0<=g<=2948 or abs(s)>2948*10**13 or (not g and s):raise ArithmeticError('cache score range differs')
        rows.append({'id':r['id'],'extension_selection_units':s,'extension_good':g,'combined_selection_units':r['score_units']+s,'combined_good':r['good_primes']+g})
    result={'status':'PASS','bank_protocol_sha256':cert.hashed(D/'bank-protocol.json'),'input_sha256':cert.hashed(input_path),'raw_sha256':cert.hashed(rawpath),'rows':rows,'wall_seconds':raw['supervision']['wall_seconds']}
    if check:
        if cert.read(out)!=result:raise ArithmeticError('cached score replay differs')
    else:checkpoint(out,result)
    print('NEW11952 CACHED32768 SCORES PASS',result['wall_seconds'],flush=True)

def scalar_roster():
    p=bank_protocol();data=cert.read(D/'cache-result.json');old=cert.read(OLD)
    if data['status']!='PASS' or data['bank_protocol_sha256']!=cert.hashed(D/'bank-protocol.json'):raise ArithmeticError('complete cache required')
    previous=old['previous_equations']+[{'address':OLD.name+':'+r['id'],'curve':r['curve']} for r in old['curves']]
    catalogue=[{'id':r['id'],'curve':r['ainvs']} for r in cert.read(CAT)['curves']]
    if len(previous)!=977 or len(catalogue)!=593 or cert.hashed(CAT)!=old['catalogue']['raw_sha256']:raise ArithmeticError('fixed equation exclusions differ')
    seen={};rows=[];skips=[];pool=[]
    def jvalue(model):
        v=cert.weierstrass_invariants(model);return v['c4']**3/v['discriminant']
    for r in previous+catalogue:
        model=tuple(map(cert.F,r['curve']));seen.setdefault(jvalue(model),[]).append((model,str(r.get('address',r.get('id')))))
    for r,s in zip(p['rows'],data['rows']):
        if r['id']!=s['id']:raise ArithmeticError('cache roster differs')
        pool.append({**r,**s})
    pool.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
    for r in pool:
        model=tuple(map(cert.F,r['model']));j=jvalue(model);matches=[name for m,name in seen.get(j,[]) if cert.isomorphic(m,model)]
        if matches:skips.append({'id':r['id'],'parameter':r['parameter'],'matches':matches});continue
        rows.append(r);seen.setdefault(j,[]).append((model,r['id']))
        if len(rows)==4096:break
    if len(rows)!=4096:raise ArithmeticError('fixed4096 eligible equations unavailable; no refill')
    return rows,skips,previous,catalogue

def scalar_freeze():
    out=D/'protocol.json'
    if out.exists():raise FileExistsError('preserve scalar4096 roster')
    replay=cert.read(D/'controller/cache-check.supervisor.json')
    if replay['outcome']!='completed' or replay['returncode']!=0:raise ArithmeticError('cache transcript replay required')
    rows,skips,previous,catalogue=scalar_roster()
    checkpoint(out,{'schema':'elliptic-curves.11952-new-annulus-scalars.v1','sources':sources(),'bank_protocol_sha256':cert.hashed(D/'bank-protocol.json'),'cache_result_sha256':cert.hashed(D/'cache-result.json'),'rows':rows,'skipped':skips,'previous_equations':previous,'catalogue_equations':catalogue,'selected_curves':64,'seconds_per_curve':20,'score_wall_seconds':1800,'validation_wall_seconds':300,'maximum_workers':2,'checkpoint_block':16,'rss_bytes':1073741824,'gp_sha256':cert.hashed(engine.GP),'fresh_validation_primes':FRESH_PRIMES,'cost_gate_rows':8,'cost_gate_maximum_projected_serial_seconds':2400,'gate':'The full32768 new-annulus retained bank passes cached extension and canonical short-score replay. Completed11952 generic17 searches have found new rank27 fibres inside the prior square. This prospective4096-curve truncation tests fresh addresses beyond131072, not the rank density of the full annulus.','ordering':'First4096 exact Q-isomorphism-distinct nonsingular models after977 prior and593 catalogue equation exclusions, under combined S1 through32749, total good count, denominator and signed numerator. Reorder only this frozen4096 by combined S1 through65521 and freeze64 finalists. All65537..131071 validation is obtained afterward and does not change selection.','scope':'Fresh scalar4099..65521 calls on exactly4096 new-annulus fibres; require cached4099..32749 agreement and8192 independent character sums. First8 runtime gate reused once, two scalar workers, no retries. After full scalar replay select64, replay exact equation exclusions and ordering, and obtain wholly disjoint validation. No new point search, rank prediction, point absence or universal novelty.'})
    print('NEW11952 FROZEN4096 DISTINCT SCALAR FIBRES',len(skips),'EXCLUSIONS',flush=True)

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP) or p['bank_protocol_sha256']!=cert.hashed(D/'bank-protocol.json') or p['cache_result_sha256']!=cert.hashed(D/'cache-result.json'):raise ArithmeticError('frozen scalar inputs differ')
    return p

def trace(row,p,create):
    folder=D/row['id'];path=folder/'raw.json';code=engine.program(row['model'])
    if create:
        if path.exists():raise FileExistsError('preserve one scalar trace call')
        c=capture([str(engine.GP),'-q','-s','256000000'],input_text=code,limits=Limits(p['seconds_per_curve'],p['rss_bytes']),log_path=folder/'gp.log',separate_stderr=True,check=False)
        checkpoint(path,{'program':code,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['program']!=code or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('scalar trace failed or censored')
    values,ms=engine.parse(raw['stdout'],row['model']);s=scoring.sums(values)
    if (s['extension_selection_units'],s['extension_good'])!=(row['extension_selection_units'],row['extension_good']):raise ArithmeticError('fresh scalar versus cached first extension differs')
    checks=[]
    for prime in (4099,32771):
        value=next(a for q,a in values if q==prime)
        if engine.direct(row['model'],prime)!=value:raise ArithmeticError('independent scalar character sum differs')
        checks.append([prime,value])
    return {'id':row['id'],'scores':s,'combined_late_units':row['combined_selection_units']+s['validation_units'],
        'combined_late_good':row['combined_good']+s['validation_good'],'direct_checks':checks,
        'cpu_ms':ms,'wall_seconds':raw['supervision']['wall_seconds'],'raw_sha256':cert.hashed(path)}

def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve one late-band score execution')
    result={'status':'RUNNING_COST_GATE','protocol_sha256':cert.hashed(D/'protocol.json'),'rows':[]};checkpoint(out,result)
    try:
        for row in p['rows'][:p['cost_gate_rows']]:result['rows'].append(trace(row,p,True));checkpoint(out,result)
        projected=sum(r['wall_seconds'] for r in result['rows'])*len(p['rows'])/p['cost_gate_rows']
        result['projected_serial_seconds']=projected
        if projected>p['cost_gate_maximum_projected_serial_seconds']:raise ArithmeticError('fixed first8 runtime gate failed')
        result['status']='RUNNING_FROZEN4096';checkpoint(out,result)
        with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
            for start in range(p['cost_gate_rows'],len(p['rows']),p['checkpoint_block']):
                result['rows']+=list(pool.map(lambda row:trace(row,p,True),p['rows'][start:start+p['checkpoint_block']]))
                checkpoint(out,result)
                if len(result['rows'])%256==8:print('NEW11952 TRACE PROGRESS',len(result['rows']),flush=True)
        result['status']='COMPLETE_FROZEN4096';checkpoint(out,result)
    except Exception as exc:
        result.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(out,result);raise

def check():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FROZEN4096' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['rows'])!=4096:raise ArithmeticError('complete score roster required')
    for row,actual in zip(p['rows'],d['rows']):
        if trace(row,p,False)!=actual:raise ArithmeticError('exact score replay differs')
    print('REPLAYED ALL4096 NEW11952 SCORES',flush=True)

def select():
    p=protocol();d=cert.read(D/'result.json');check=cert.read(D/'controller/check.supervisor.json')
    if OUT.exists():raise FileExistsError('preserve frozen64 late-band finalists')
    if d['status']!='COMPLETE_FROZEN4096' or len(d['rows'])!=4096 or check['outcome']!='completed' or check['returncode']!=0:raise ArithmeticError('complete independent score rerun required')
    rows=[]
    for row,score in zip(p['rows'],d['rows']):
        if row['id']!=score['id']:raise ArithmeticError('score roster differs')
        rows.append({**row,**score})
    rows.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
    checkpoint(OUT,{'schema':'elliptic-curves.11952-new-annulus-finalists.v1','status':'PASS_FROZEN64_SELECTION',
        'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'scores_sha256':cert.hashed(D/'result.json'),
        'score_replay_sha256':cert.hashed(D/'controller/check.supervisor.json'),'selected':rows[:64],
        'claim_boundary':p['scope']})
    print('FROZEN64 NEW11952 FINALISTS',flush=True)

def fresh_program(model):
    return engine.program(model).replace('forprime(p=4099,65521','forprime(p=65537,131071')

def validate(check=False):
    p=protocol();selection=cert.read(OUT);out=D/'fresh-validation.json'
    if not check and out.exists():raise FileExistsError('preserve disjoint validation')
    def one(row):
        path=D/row['id']/'fresh-raw.json';code=fresh_program(row['model'])
        if not check:
            if path.exists():raise FileExistsError('preserve one disjoint trace call')
            c=capture([str(engine.GP),'-q','-s','256000000'],input_text=code,limits=Limits(p['seconds_per_curve'],p['rss_bytes']),log_path=path.with_suffix('.log'),separate_stderr=True,check=False)
            checkpoint(path,{'program':code,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
        raw=cert.read(path);lines=raw['stdout'].splitlines();primes=p['fresh_validation_primes'];a,b=map(int,row['model'][3:]);units=good=0;values=[]
        if raw['program']!=code or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or len(lines)!=len(primes)+2 or lines[-1]!='DONE' or not lines[-2].startswith('MS|'):raise ArithmeticError('disjoint trace execution/framing differs')
        for q,line in zip(primes,lines):
            v=line.split('|');bad=(4*a*a*a+27*b*b)%q==0
            if len(v)!=(2 if bad else 3) or v[0]!=('B' if bad else 'T') or int(v[1])!=q:raise ArithmeticError('disjoint reduction framing differs')
            t=None if bad else int(v[2]);values.append([q,t])
            if t is not None:
                if t*t>4*q:raise ArithmeticError('disjoint Hasse bound differs')
                units+=round((2-t)*math.log(q)/(q+1-t)*10**12);good+=1
        checks=[]
        for i in (0,len(primes)//2,len(primes)-1):
            q,t=values[i]
            if engine.direct(row['model'],q)!=t:raise ArithmeticError('fresh direct character sum differs')
            checks.append([q,t])
        return {'id':row['id'],'status':'PASS','validation_units':units,'validation_good':good,'direct_checks':checks,'raw_sha256':cert.hashed(path)}
    with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:rows=list(pool.map(one,selection['selected']))
    result={'status':'PASS','selection_sha256':cert.hashed(OUT),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,
        'claim_boundary':'Wholly disjoint65537..131071 validation recorded after all64 finalists are frozen. These values do not change selection or map geometry, and are not rank bounds.'}
    if check:
        if cert.read(out)!=result:raise ArithmeticError('disjoint validation replay differs')
    else:checkpoint(out,result)
    print('NEW11952 FINALISTS FRESH VALIDATION64 PASS',flush=True)


def selection_check():
    p=protocol();rows,skips,previous,catalogue=scalar_roster()
    if (p['rows'],p['skipped'],p['previous_equations'],p['catalogue_equations'])!=(rows,skips,previous,catalogue):raise ArithmeticError('exact scalar roster and equation exclusions differ')
    data=cert.read(D/'result.json');saved=cert.read(OUT);merged=[]
    if data['status']!='COMPLETE_FROZEN4096' or len(data['rows'])!=4096:raise ArithmeticError('complete frozen scores required')
    for r,t in zip(rows,data['rows']):
        if r['id']!=t['id']:raise ArithmeticError('scalar score roster differs')
        merged.append({**r,**t})
    merged.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
    if saved['selected']!=merged[:64] or saved['status']!='PASS_FROZEN64_SELECTION' or saved['sources']!=sources() or saved['protocol_sha256']!=cert.hashed(D/'protocol.json') or saved['scores_sha256']!=cert.hashed(D/'result.json') or saved['score_replay_sha256']!=cert.hashed(D/'controller/check.supervisor.json'):raise ArithmeticError('exact64 selection differs')
    print('NEW11952 EXACT4096 EXCLUSIONS AND64 ORDERING PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','cache','cache-check','scalar-freeze','run','check','select','selection-check','validate','validation-check']);a=p.parse_args()
    if a.stage in ('cache','cache-check'):cached(a.stage=='cache-check')
    elif a.stage=='scalar-freeze':scalar_freeze()
    elif a.stage=='selection-check':selection_check()
    elif a.stage=='validation-check':validate(True)
    else:globals()[a.stage]()

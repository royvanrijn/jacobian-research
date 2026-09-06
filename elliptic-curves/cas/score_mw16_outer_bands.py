#!/usr/bin/env python3
"""Target-free extended scores for15360 fresh outer MW16 parameters."""
import argparse,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_mw16_outer_bands as scan
import benchmark_mw16_extended_prime_traces as engine
import extend_mw16_retained_prime_scores as models
import extend_retained_r17_prime_scores as scoring
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=scan.ROOT;CAS=scan.CAS;ART=scan.ART
D=ROOT/'artifacts/local/elliptic-curves/mw16-outer-band-scores-v1'
OUT=ART/'mw16_outer_band_selection_v1.json'
FRESH_PRIMES=[q for q in _primes_up_to(131071) if q>=65537]

def sources():
    paths=[Path(__file__).resolve(),Path(scan.__file__),Path(engine.__file__),Path(models.__file__),Path(scoring.__file__),Path(cert.__file__),scan.spec.ATLAS,scan.D/'protocol.json',scan.D/'result.json',scan.D/'replay.json',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py',CAS/'mod2_reduction_independence.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve outer score protocol')
    upstream=cert.read(scan.D/'controller/ledger.json');replay=cert.read(scan.D/'replay.json');raw=cert.read(scan.D/'result.json')
    if upstream['status']!='PASS' or replay['status']!='PASS' or replay['result_sha256']!=cert.hashed(scan.D/'result.json'):raise ArithmeticError('all outer population gates required')
    families={r['fibration_id']:r for r in cert.read(scan.spec.ATLAS)['families']};pools=[]
    for band in (1,2,3):
        for family in sorted(families):
            pool=[{**r,'band':band,'family':family,'slice_id':s['id']} for s in raw['shards'] if s['band']==band and s['family']==family for r in s['rows']]
            pool.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
            if len(pool)!=1024 or len({r['parameter'] for r in pool})!=1024:raise ArithmeticError('fixed band-family1024 differs')
            pools.append([{**r,'id':f'b{band}-{family}-{i:04}','retained_index':i,'model':models.model_at(families[family],r['parameter'])} for i,r in enumerate(pool)])
    rows=[pool[i] for i in range(1024) for pool in pools]
    if len(rows)!=15360:raise ArithmeticError('fixed15360 roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-outer-band-scores.v1','sources':sources(),'rows':rows,'selected_curves':60,'seconds_per_curve':20,'score_wall_seconds':3600,'validation_wall_seconds':300,'maximum_workers':2,'checkpoint_block':30,'rss_bytes':1073741824,'gp_sha256':cert.hashed(engine.GP),'fresh_validation_primes':FRESH_PRIMES,'cost_gate_rows':15,'cost_gate_maximum_projected_serial_seconds':4800,'gate':'All three outer bands and every signed MW16 frame passed small exact annular scanner comparisons, finite runtime gates, independent primitive counts and15360 canonical562-prime score replays. Extend every retained address before selecting four per family per band. This is a new-parameter incidence experiment with generic16-only point work planned after complete score and fresh-prime validation replay.','ordering':'Four per family per band, ordered by combined quantized S1 through65521, total good-prime count, denominator, signed numerator. Deduplicate selected equations over Q across band/family groups in fixed band then family order, filling each group from its same frozen1024 pool. No catalogue, public points, previous ranks or disjoint validation enters selection.','scope':'Finite15360 survivors from286812899 primitive addresses in disjoint stratified bands beyond4096 and through262144. First15 rows are one per band/family and form a runtime gate, reused once. Two direct character sums on every scalar trace. Wholly disjoint65537..131071 validation on frozen60 finalists. No retries, broadening, rank prediction, universal novelty or claim of complete parameter coverage.'})
    print('FROZEN15360 OUTER MW16 TRACE ROSTER',flush=True)
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP):raise ArithmeticError('frozen late-band inputs changed')
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
    checks=[]
    for prime in (4099,32771):
        value=next(a for q,a in values if q==prime)
        if engine.direct(row['model'],prime)!=value:raise ArithmeticError('independent scalar character sum differs')
        checks.append([prime,value])
    return {'id':row['id'],'scores':s,'combined_late_units':row['score_units']+s['extension_selection_units']+s['validation_units'],
        'combined_late_good':row['good_primes']+s['extension_good']+s['validation_good'],'direct_checks':checks,
        'cpu_ms':ms,'wall_seconds':raw['supervision']['wall_seconds'],'raw_sha256':cert.hashed(path)}

def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve one late-band score execution')
    result={'status':'RUNNING_COST_GATE','protocol_sha256':cert.hashed(D/'protocol.json'),'rows':[]};checkpoint(out,result)
    try:
        for row in p['rows'][:p['cost_gate_rows']]:result['rows'].append(trace(row,p,True));checkpoint(out,result)
        projected=sum(r['wall_seconds'] for r in result['rows'])*len(p['rows'])/p['cost_gate_rows']
        result['projected_serial_seconds']=projected
        if projected>p['cost_gate_maximum_projected_serial_seconds']:raise ArithmeticError('fixed first15 runtime gate failed')
        result['status']='RUNNING_FROZEN15360';checkpoint(out,result)
        with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
            for start in range(p['cost_gate_rows'],len(p['rows']),p['checkpoint_block']):
                result['rows']+=list(pool.map(lambda row:trace(row,p,True),p['rows'][start:start+p['checkpoint_block']]))
                checkpoint(out,result)
                if len(result['rows'])%300==15:print('MW16 OUTER TRACE PROGRESS',len(result['rows']),flush=True)
        result['status']='COMPLETE_FROZEN15360';checkpoint(out,result)
    except Exception as exc:
        result.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(out,result);raise

def check():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FROZEN15360' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['rows'])!=15360:raise ArithmeticError('complete score roster required')
    for row,actual in zip(p['rows'],d['rows']):
        if trace(row,p,False)!=actual:raise ArithmeticError('exact score replay differs')
    print('REPLAYED ALL15360 MW16 OUTER SCORES',flush=True)

def selected_rows(p,d):
    pools={};seen={};selected=[];skips=[]
    for row,score in zip(p['rows'],d['rows']):
        if row['id']!=score['id']:raise ArithmeticError('score roster differs')
        pools.setdefault((row['band'],row['family']),[]).append({**row,**score})
    for key,pool in sorted(pools.items()):
        pool.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
        count=0
        for row in pool:
            model=tuple(map(cert.F,row['model']));v=cert.weierstrass_invariants(model);j=v['c4']**3/v['discriminant']
            matches=[name for m,name in seen.get(j,[]) if cert.isomorphic(m,model)]
            if matches:skips.append({'id':row['id'],'matches':matches});continue
            selected.append(row);seen.setdefault(j,[]).append((model,row['id']));count+=1
            if count==4:break
        if count!=4:raise ArithmeticError('insufficient distinct equations; no population expansion')
    if len(selected)!=60:raise ArithmeticError('fixed60 differs')
    return selected,skips

def select():
    p=protocol();d=cert.read(D/'result.json');gate=cert.read(D/'controller/check.supervisor.json')
    if OUT.exists():raise FileExistsError('preserve frozen60')
    if d['status']!='COMPLETE_FROZEN15360' or len(d['rows'])!=15360 or gate['outcome']!='completed' or gate['returncode']!=0:raise ArithmeticError('complete exact score replay required')
    selected,skips=selected_rows(p,d)
    checkpoint(OUT,{'schema':'elliptic-curves.mw16-outer-band-finalists.v1','status':'PASS_FROZEN60_SELECTION','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'scores_sha256':cert.hashed(D/'result.json'),'score_replay_sha256':cert.hashed(D/'controller/check.supervisor.json'),'selected':selected,'skips':skips,'claim_boundary':p['scope']})
    print('FROZEN60 DISTINCT OUTER MW16 FINALISTS',flush=True)

def selection_check():
    p=protocol();d=cert.read(D/'result.json');saved=cert.read(OUT);selected,skips=selected_rows(p,d)
    if saved['selected']!=selected or saved['skips']!=skips or saved['protocol_sha256']!=cert.hashed(D/'protocol.json') or saved['scores_sha256']!=cert.hashed(D/'result.json') or saved['sources']!=sources():raise ArithmeticError('selection replay differs')
    print('REPLAYED60 EXACT OUTER EQUATION SELECTION',flush=True)

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
        'claim_boundary':'Wholly disjoint65537..131071 validation recorded after all60 finalists are frozen. These values do not change selection or map geometry, and are not rank bounds.'}
    if check:
        if cert.read(out)!=result:raise ArithmeticError('disjoint validation replay differs')
    else:checkpoint(out,result)
    print('MW16 OUTER FINALISTS FRESH VALIDATION60 PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','check','select','validate','validation-check','selection_check']);a=p.parse_args()
    if a.stage=='validation-check':validate(True)
    else:globals()[a.stage]()

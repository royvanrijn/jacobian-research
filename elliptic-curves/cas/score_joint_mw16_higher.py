#!/usr/bin/env python3
"""Target-free scalar selection on10240 new higher MW16 equations."""
import argparse,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_joint_mw16_higher_annuli as scan
import benchmark_mw16_extended_prime_traces as engine
import extend_mw16_retained_prime_scores as models
import extend_retained_r17_prime_scores as scoring
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=scan.ROOT;CAS=scan.CAS;ART=scan.ART
D=ROOT/'artifacts/local/elliptic-curves/joint-mw16-higher-scores-v1'
OUT=ART/'joint_mw16_higher_selection_v1.json'
FRESH_PRIMES=[q for q in _primes_up_to(131071) if q>=65537]

def sources():
    paths=[Path(__file__).resolve(),Path(scan.__file__),Path(engine.__file__),Path(models.__file__),Path(scoring.__file__),Path(cert.__file__),scan.old.spec.ATLAS,scan.D/'protocol.json',scan.D/'result.json',scan.D/'replay.json',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py',CAS/'mod2_reduction_independence.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prospective_rows(raw):
    families={r['fibration_id']:r for r in cert.read(scan.old.spec.ATLAS)['families']};pools=[];seen={};skips=[]
    for band in (2,3):
        for family in sorted(families):
            candidates=[{**r,'band':band,'family':family,'slice_id':s['id']} for s in raw['shards'] if s['band']==band and s['family']==family for r in s['rows']]
            candidates.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
            if len(candidates)!=8192 or len({r['parameter'] for r in candidates})!=8192:raise ArithmeticError('fixed8192 prospective band/family pool required')
            pool=[]
            for i,r in enumerate(candidates):
                row={**r,'id':f'higher-b{band}-{family}-{i:04}','retained_index':i,'model':models.model_at(families[family],r['parameter'])}
                model=tuple(map(cert.F,row['model']));inv=cert.weierstrass_invariants(model);j=inv['c4']**3/inv['discriminant']
                matches=[name for other,name in seen.get(j,[]) if cert.isomorphic(model,other)]
                if matches:skips.append({'id':row['id'],'matches':matches});continue
                pool.append(row);seen.setdefault(j,[]).append((model,row['id']))
                if len(pool)==1024:break
            if len(pool)!=1024:raise ArithmeticError('insufficient within-roster distinct equations; no expansion')
            pools.append(pool)
    return [pool[i] for i in range(1024) for pool in pools],skips

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher scalar roster')
    replay=cert.read(scan.D/'replay.json');raw=cert.read(scan.D/'result.json')
    if cert.read(scan.D/'controller/ledger.json')['status']!='PASS' or replay['status']!='PASS' or replay['result_sha256']!=cert.hashed(scan.D/'result.json'):raise ArithmeticError('all full-score higher-slice proofs required')
    rows,skips=prospective_rows(raw)
    if len(rows)!=10240:raise ArithmeticError('fixed10240 higher scalar cases required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-joint-higher-scores.v1','sources':sources(),'rows':rows,'prospective_skips':skips,
        'selected_curves':60,'seconds_per_curve':20,'score_wall_seconds':2400,'validation_wall_seconds':300,
        'maximum_workers':4,'checkpoint_block':80,'rss_bytes':1073741824,'gp_sha256':cert.hashed(engine.GP),
        'fresh_validation_primes':FRESH_PRIMES,'cost_gate_rows':20,'cost_gate_maximum_projected_serial_seconds':6400,
        'gate':'All191215782 primitive addresses in twenty untouched higher slices have been scored through32749 before retention, with all81920 retained scores replayed. The complete five-family cache matches all40960 earlier scalar cases. Choose1024 within-roster distinct equations per higher band/family from the fixed8192-address signed pool, giving equal scalar exposure without public-record targets.',
        'ordering':'Candidate1024 pools use combined S1 through32749,good count,denominator,signed numerator in band/family order; only prospective Q-isomorphism deduplication is permitted. Final six per band/family use S1 through65521 with the same tie rules. Disjoint65537..131071 validation cannot alter selection.',
        'scope':'Exactly10240 distinct prospective equations beyond16384 and through262144,1024 per band/family. The first20 rows include two per band/family and form a fixed serial-cost gate reused once. Four workers,20seconds per curve,2400seconds total,80-row checkpoints. Require scalar agreement with the first cached extension and20480 independent character sums, then complete score replay before60 finalists and disjoint validation. No catalogue,known record equations,parameters,points,ranks,j-invariants or jump labels enter selection/execution. No prior-point-outcome filtering, retry, compact rescan, adaptive refill or rank inference.'})
    print('FROZEN10240 DISTINCT HIGHER MW16 SCALAR CASES',flush=True)

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP):raise ArithmeticError('frozen higher scalar sources changed')
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
    if (s['extension_selection_units'],s['extension_good'])!=(row['extension_selection_units'],row['extension_good']):raise ArithmeticError('cached first extension differs from independent scalar traces')
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
        if projected>p['cost_gate_maximum_projected_serial_seconds']:raise ArithmeticError('fixed first20 runtime gate failed')
        result['status']='RUNNING_FROZEN10240';checkpoint(out,result)
        with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
            for start in range(p['cost_gate_rows'],len(p['rows']),p['checkpoint_block']):
                result['rows']+=list(pool.map(lambda row:trace(row,p,True),p['rows'][start:start+p['checkpoint_block']]))
                checkpoint(out,result)
                if len(result['rows'])%800==20:print('MW16 HIGHER JOINT TRACE PROGRESS',len(result['rows']),flush=True)
        result['status']='COMPLETE_FROZEN10240';checkpoint(out,result)
    except Exception as exc:
        result.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(out,result);raise

def check():
    p=protocol();d=cert.read(D/'result.json')
    rows,skips=prospective_rows(cert.read(scan.D/'result.json'))
    if p['rows']!=rows or p['prospective_skips']!=skips:raise ArithmeticError('prospective full-score selection replay differs')
    if d['status']!='COMPLETE_FROZEN10240' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['rows'])!=10240:raise ArithmeticError('complete score roster required')
    for row,actual in zip(p['rows'],d['rows']):
        if trace(row,p,False)!=actual:raise ArithmeticError('exact higher score replay differs')
    print('REPLAYED ALL10240 MW16 HIGHER JOINT SCORES',flush=True)

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
            if count==6:break
        if count!=6:raise ArithmeticError('insufficient distinct equations; no population expansion')
    if len(selected)!=60:raise ArithmeticError('fixed60 differs')
    return selected,skips

def select():
    p=protocol();d=cert.read(D/'result.json');gate=cert.read(D/'controller/check.supervisor.json')
    if OUT.exists():raise FileExistsError('preserve frozen60')
    if d['status']!='COMPLETE_FROZEN10240' or len(d['rows'])!=10240 or gate['outcome']!='completed' or gate['returncode']!=0:raise ArithmeticError('complete exact higher score replay required')
    selected,skips=selected_rows(p,d)
    checkpoint(OUT,{'schema':'elliptic-curves.mw16-joint-higher-finalists.v1','status':'PASS_FROZEN60_SELECTION','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'scores_sha256':cert.hashed(D/'result.json'),'score_replay_sha256':cert.hashed(D/'controller/check.supervisor.json'),'selected':selected,'skips':skips,'claim_boundary':p['scope']})
    print('FROZEN60 DISTINCT FRESH MW16 FINALISTS',flush=True)

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
    print('MW16 HIGHER JOINT FINALISTS FRESH VALIDATION60 PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','check','select','validate','validation-check','selection_check']);a=p.parse_args()
    if a.stage=='validation-check':validate(True)
    else:globals()[a.stage]()

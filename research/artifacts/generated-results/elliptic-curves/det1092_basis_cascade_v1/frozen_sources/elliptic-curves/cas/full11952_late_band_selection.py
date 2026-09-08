#!/usr/bin/env python3
"""Fixed4096 fresh equations, second prime-band selection, then64 finalists."""
import argparse,json,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import score_full11952_retained as parent
import benchmark_r17_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART
D=ROOT/'artifacts/local/elliptic-curves/full11952-late-band-selection-v1'
OUT=ART/'full11952_late_band_selection_v1.json'
OLD=ART/'full11952_64_r17_results_v1.json'
CAT=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1/database.json'
FRESH_PRIMES=[q for q in _primes_up_to(131071) if q>=65537]

def sources():
    paths=[Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(engine.__file__).resolve(),Path(scoring.__file__).resolve(),OLD,CAT,parent.OUT,parent.META,parent.EXT,parent.D/'selected-validation/result.json',ART/'full11952_64_point_portable_replay_v1.json',parent.extension.spec.ATLAS]
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve late-band selection protocol')
    old=cert.read(OLD);proof=cert.read(ART/'full11952_64_point_portable_replay_v1.json')
    if proof['status']!='PASS' or proof['logical_stages']!=194 or len(old['curves'])!=64:raise ArithmeticError('previous full cohort must be terminal and replayed')
    previous=old['previous_equations']+[{'address':OLD.name+':'+r['id'],'curve':r['curve']} for r in old['curves']]
    if len(previous)!=853:raise ArithmeticError('fixed853 previous equations required')
    catalogue=[{'id':r['id'],'curve':r['ainvs']} for r in cert.read(CAT)['curves']]
    if len(catalogue)!=593:raise ArithmeticError('fixed593 catalogue equations required')
    roster=[];skips=[];seen={}
    def jvalue(model):
        v=cert.weierstrass_invariants(model);return v['c4']**3/v['discriminant']
    for r in previous+catalogue:
        model=tuple(map(cert.F,r['curve']));seen.setdefault(jvalue(model),[]).append((model,str(r.get('address',r.get('id')))))
    p=parent.protocol();prefix=parent.ranked_prefix(p)
    family=next(f for f in cert.read(parent.extension.spec.ATLAS)['families'] if f['family']=='11952')
    for position,(u,g,nd,nn,i,s,sg,e,eg) in enumerate(prefix):
        n,d=-nn,-nd;t=str(cert.F(n,d));model=tuple(map(cert.F,parent.extension.model_at(family,t)));j=jvalue(model)
        matches=[name for m,name in seen.get(j,[]) if cert.isomorphic(m,model)]
        if matches:skips.append({'retained_index':i,'parameter':t,'matches':matches});continue
        row={'id':f'11952-{i:07}','family':'11952','retained_index':i,'parameter':t,'numerator':n,'denominator':d,
             'score_units':s,'good_primes':sg,'extension_selection_units':e,'extension_good':eg,
             'combined_selection_units':u,'combined_good':g,'model':list(map(str,model)),
             'old_prefix_position':position+1}
        roster.append(row);seen.setdefault(j,[]).append((model,row['id']))
        if len(roster)==4096:break
    if len(roster)!=4096:raise ArithmeticError('insufficient fixed prefix; no automatic expansion')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.full11952-late-band-selection.v1','sources':sources(),
        'rows':roster,'previous_equations':previous,'catalogue_equations':catalogue,'skipped':skips,
        'old_prefix_rows':32768,'old_prefix_sha256':digest(prefix),'selected_curves':64,
        'seconds_per_curve':20,'score_wall_seconds':1800,'validation_wall_seconds':300,
        'maximum_workers':2,'checkpoint_block':16,'rss_bytes':1073741824,
        'gp_sha256':cert.hashed(engine.GP),'fresh_validation_primes':FRESH_PRIMES,
        'cost_gate_rows':8,'cost_gate_maximum_projected_serial_seconds':2400,
        'gate':'The completed64 cohort has unused32771..65521 traces. Combining these with selection through32749 moves the known29 control from29th to1st, and puts five of the six measured bounds>=19 in the top6. This retrospective within-selected-cohort observation motivates one fresh selector trial; it is not a calibrated rank predictor or sensitivity estimate. Existing scalar benchmarks also support the finite runtime gate.',
        'ordering':'Take the first4096 nonsingular Q-isomorphism-distinct equations in the existing32768-row extended-score prefix, excluding all853 previous cohort equations and593 pinned catalogue equations. For these frozen4096, select64 by combined quantized S1 through65521, total good-prime count descending, denominator then signed numerator ascending. No point, public rank or new validation value enters selection.',
        'scope':'No new parameter scan or full-million second-band cache. The4096-curve roster is a declared further truncation of the retained population. Fresh scalar traces4099..65521 must match the earlier cached4099..32749 sums, with two direct character sums per curve. The first8 frozen rows provide a runtime gate and are reused once; no retries. After score replay freeze64, then obtain wholly disjoint65537..131071 validation with three direct sums per finalist. A separate point protocol is required before any point search; no automatic refill or adaptive wave.'})
    print('FROZEN LATE-BAND4096',len(skips),'EXCLUSIONS',roster[-1]['old_prefix_position'],flush=True)

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
                if len(result['rows'])%256==8:print('LATE-BAND TRACE PROGRESS',len(result['rows']),flush=True)
        result['status']='COMPLETE_FROZEN4096';checkpoint(out,result)
    except Exception as exc:
        result.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(out,result);raise

def check():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FROZEN4096' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['rows'])!=4096:raise ArithmeticError('complete score roster required')
    for row,actual in zip(p['rows'],d['rows']):
        if trace(row,p,False)!=actual:raise ArithmeticError('exact score replay differs')
    print('REPLAYED ALL4096 LATE-BAND SCORES',flush=True)

def select():
    p=protocol();d=cert.read(D/'result.json');check=cert.read(D/'check.supervisor.json')
    if OUT.exists():raise FileExistsError('preserve frozen64 late-band finalists')
    if d['status']!='COMPLETE_FROZEN4096' or len(d['rows'])!=4096 or check['outcome']!='completed' or check['returncode']!=0:raise ArithmeticError('complete independent score rerun required')
    rows=[]
    for row,score in zip(p['rows'],d['rows']):
        if row['id']!=score['id']:raise ArithmeticError('score roster differs')
        rows.append({**row,**score})
    rows.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
    checkpoint(OUT,{'schema':'elliptic-curves.full11952-late-band-finalists.v1','status':'PASS_FROZEN64_SELECTION',
        'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'scores_sha256':cert.hashed(D/'result.json'),
        'score_replay_sha256':cert.hashed(D/'check.supervisor.json'),'selected':rows[:64],
        'claim_boundary':p['scope']})
    print('FROZEN64 LATE-BAND FINALISTS',flush=True)

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
    print('LATE-BAND FINALISTS FRESH VALIDATION64 PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','check','select','validate','validation-check']);a=p.parse_args()
    if a.stage=='validation-check':validate(True)
    else:globals()[a.stage]()

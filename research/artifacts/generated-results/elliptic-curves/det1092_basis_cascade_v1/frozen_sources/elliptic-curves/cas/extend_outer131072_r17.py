#!/usr/bin/env python3
"""Fixed extended-prime scoring of the saved outer32768-to131072 candidates."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import scan_outer131072_r17_v2 as parent
import benchmark_r17_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=parent.LOCAL/'outer131072-r17-extended-v1';OLD=ART/'compact192_r17_results_v1.json'
def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(scoring.__file__).resolve(),spec.ATLAS,OLD)}}
def model_at(f,t):
    t=cert.F(t);q=t.denominator;model=[cert.F(0)]*3+[spec.polynomial(f['A_coefficients_low_to_high'],t)*q**8,spec.polynomial(f['B_coefficients_low_to_high'],t)*q**12]
    if any(x.denominator!=1 for x in model) or 4*model[3]**3+27*model[4]**2==0:raise ArithmeticError('integral nonsingular specialization required')
    return list(map(str,model))
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher-height trace extension')
    p=parent.protocol();d=cert.read(parent.D/'result.json');replay=cert.read(parent.D/'replay.supervisor.json')
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or replay['outcome']!='completed' or replay['returncode']!=0 or len(d['rows'])!=6144:raise ArithmeticError('complete6144-address short-score replay required')
    atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};rows=[{**r,'model':model_at(atlas[r['family']],r['parameter'])} for r in d['rows'] if max(abs(r['numerator']),r['denominator'])>32768];samples=[next(i for i,r in enumerate(rows) if r['family']==family and r['retained_index']>=4) for family in sorted(atlas)]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.outer131072-r17-prime-extension.v1','original_retained_rows':6144,'outer_rows':len(rows),'excluded_inner_rows':6144-len(rows),'maximum_outer_rows':6144,'sources':sources(),'parent_protocol_sha256':cert.hashed(parent.D/'protocol.json'),'parent_result_sha256':cert.hashed(parent.D/'result.json'),'parent_replay_sha256':cert.hashed(parent.D/'replay.supervisor.json'),'rows':rows,'prime_roster':engine.PRIMES,'gp_sha256':cert.hashed(engine.GP),'benchmark_indices':samples,'direct_validation_indices':[0,1,2,len(engine.PRIMES)//4,len(engine.PRIMES)//2,3*len(engine.PRIMES)//4,len(engine.PRIMES)-2,len(engine.PRIMES)-1],'benchmark_seconds_per_curve':20,'seconds_per_curve':5,'rss_bytes':536870912,'maximum_workers':2,'outer_seconds':1800,'checkpoint_block_size':16,'gate':'All twelve fixed H131072 modulo1024 slices and6144 returned562-prime scores have exact replay. Exclude every address inside the earlier32768 square before extending traces. Extend every remaining outer address using the unchanged trace policy, after six fixed cost/character-sum checks pass: the first outer retained index at least4 in each family. Public record parameters, points, measured ranks and catalogue labels do not enter models, traces or ordering.','selection':'In lexicographic family order select eight per family by combined quantized S1 through32749, then combined selection-band good count, denominator and signed numerator. Exclude rational isomorphs of all741 previously measured equations in the completed compact192 proof chain and of earlier chosen rows, using equations only. Validation32771..65521 is excluded even from ties. No point-attempt failures or subsequent catalogue matches are refilled.','future_point_scope':'At most48 outer-only generic17-point attempts, eight per family, using all43/49 exact generic maximum classes after every outer trace and the fixed selection replay. A separate point protocol must freeze all maps before any points, height125000 and ten seconds per chart, with a provisional28 target pending exact proof. No point search occurs during this extension.','boundaries':'Finite selected denominator slices, with a short-prime retention bottleneck. No full131072 coverage, density estimate, rank classifier, upper bound or universal novelty. All failures/censored raw outputs retained; no automatic retry.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['parent_protocol_sha256']!=cert.hashed(parent.D/'protocol.json') or p['parent_result_sha256']!=cert.hashed(parent.D/'result.json') or p['gp_sha256']!=cert.hashed(engine.GP):raise ArithmeticError('frozen higher trace inputs changed')
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
    print('OUTER131072 TRACE BENCHMARK48 DIRECT CHECKS PASS',flush=True)
def checked_row(row,p,create):
    path=D/row['family']/f"candidate-{row['retained_index']:04}"/'raw.json';(values,ms),raw=call(row,p,path,5,create);s=scoring.sums(values)
    return {**row,**s,'combined_selection_units':row['score_units']+s['extension_selection_units'],'combined_good':row['good_primes']+s['extension_good'],'cpu_ms':ms,'raw_sha256':cert.hashed(path)}
def selection(rows):
    old=cert.read(OLD);previous=old['previous_equations']+[{'curve':r['curve']} for r in old['curves']]
    if len(previous)!=741:raise ArithmeticError('fixed741 previous equations required')
    def j_value(model):
        inv=cert.weierstrass_invariants(model)
        return inv['c4']**3/inv['discriminant']
    seen={}
    for row in previous:
        model=tuple(map(cert.F,row['curve']));seen.setdefault(j_value(model),[]).append(model)
    selected={}
    for family in sorted({r['family'] for r in rows}):
        chosen=[]
        candidates=sorted((r for r in rows if r['family']==family),key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
        for row in candidates:
            if max(abs(row['numerator']),row['denominator'])<=32768:raise ArithmeticError('inner parameter leaked into outer trace cohort')
            model=tuple(map(cert.F,row['model']));j=j_value(model)
            if any(cert.isomorphic(model,q) for q in seen.get(j,[])):continue
            chosen.append(row['retained_index']);seen.setdefault(j,[]).append(model)
            if len(chosen)==8:break
        if len(chosen)!=8:raise ArithmeticError('eight eligible distinct outer curves per family required')
        selected[family]=chosen
    return selected

def run():
    p=protocol();out=D/'result.json';benchmark(True)
    if out.exists():raise FileExistsError('preserve full trace batch')
    d={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,d)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for start in range(0,len(p['rows']),16):
            d['rows']+=list(pool.map(lambda r:checked_row(r,p,True),p['rows'][start:start+16]));checkpoint(out,d)
            if len(d['rows'])%256==0:print('OUTER131072 TRACE EXTENSION',len(d['rows']),'of',len(p['rows']),flush=True)
    d['selection']=selection(d['rows']);d['status']='COMPLETE_FROZEN_TRACE_EXTENSION';checkpoint(out,d)
def replay():
    p=protocol();d=cert.read(D/'result.json');base=cert.read(parent.D/'result.json');atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};expected=[{**r,'model':model_at(atlas[r['family']],r['parameter'])} for r in base['rows'] if max(abs(r['numerator']),r['denominator'])>32768]
    if p['rows']!=expected or d['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or d['protocol_hash']!=digest(p) or len(d['rows'])!=len(p['rows']):raise ArithmeticError('full6144 roster differs')
    for row,got in zip(p['rows'],d['rows']):
        if checked_row(row,p,False)!=got:raise ArithmeticError('exact trace/score replay differs')
    if d['selection']!=selection(d['rows']):raise ArithmeticError('fixed24 selector differs')
    print('REPLAYED OUTER TRACE ROSTERS AND FIXED48 SELECTION',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()
